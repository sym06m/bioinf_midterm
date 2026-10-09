# Project 18: short vs long reads, E. coli K-12 MG1655.
# NOTE: this workflow mirrors the Colab notebooks (notebooks/) but has NOT been executed end to end yet.
# Run:  snakemake --cores 2 -n        (dry run)     |     snakemake --cores 2     (everything)
configfile: "config.yaml"

S, L = config["short_acc"], config["long_acc"]
GS, SEED, T = config["genome_size"], config["seed"], config["threads"]
COVS, REP, REPCOVS = config["covs"], config["rep_seed"], config["rep_covs"]
HS, HL = config["hybrid"]["short"], config["hybrid"]["long"]
FLYE = config["flye_mode"]

# plain-string helpers (no f-strings: older Snakemake parsers fail on them with newer Python)
SR1 = "data/raw/" + S + "_1.fastq.gz"
SR2 = "data/raw/" + S + "_2.fastq.gz"
LR = "data/raw/" + L + ".fastq.gz"
HSR1 = "sub/short_" + str(HS) + "x_R1.fq.gz"
HSR2 = "sub/short_" + str(HS) + "x_R2.fq.gz"
HLR = "sub/long_" + str(HL) + "x_seed" + str(SEED) + ".fq.gz"
HYB = "asm/hybridspades_s" + str(HS) + "_l" + str(HL)
POL = "asm/polish_s" + str(HS) + "_l" + str(HL)

wildcard_constraints:
    cov=r"\d+", seed=r"\d+"

ASSEMBLIES = (
    ["asm/short_{}x/contigs.fasta".format(c) for c in COVS]
    + ["asm/long_{}x/assembly.fasta".format(c) for c in COVS]
    + ["asm/long_{}x-rep{}/assembly.fasta".format(c, REP) for c in REPCOVS]
    + [HYB + "/contigs.fasta", POL + "/polished.fasta"]
)

rule all:
    input: "results/quast/transposed_report.tsv", "results/indel_context.csv", "figures/fig1_assembly_metrics.png"

rule fetch_short:
    output: SR1, SR2
    shell: "python scripts/fetch_ena.py {S} data/raw"

rule fetch_long:
    output: LR
    shell: "python scripts/fetch_ena.py {L} data/raw"

rule fetch_ref:
    output: "data/ref/ref.fasta"
    params: url=config["reference_url"]
    shell: "mkdir -p data/ref && wget -q -O data/ref/ref.fna.gz {params.url} && gunzip -c data/ref/ref.fna.gz > {output}"

rule qc_short:
    input: r1=SR1, r2=SR2
    output: r1="data/qc/short_R1.trim.fq.gz", r2="data/qc/short_R2.trim.fq.gz", j="results/qc/fastp.json"
    threads: T
    shell:
        """
        mkdir -p results/qc && fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe -q 20 -l 100 --cut_tail --cut_tail_mean_quality 20 -w {threads} -j {output.j} -h results/qc/fastp.html
        """

rule qc_long:
    input: LR
    output: "results/qc/nanoplot/NanoStats.txt"
    threads: T
    shell: "NanoPlot --fastq {input} -t {threads} -o results/qc/nanoplot --loglength"

rule sub_short:
    input: r1="data/qc/short_R1.trim.fq.gz", r2="data/qc/short_R2.trim.fq.gz"
    output: r1="sub/short_{cov}x_R1.fq.gz", r2="sub/short_{cov}x_R2.fq.gz"
    shell: "rasusa reads -c {wildcards.cov} -g {GS} -s {SEED} -o {output.r1} -o {output.r2} {input.r1} {input.r2}"

rule sub_long:
    input: LR
    output: "sub/long_{cov}x_seed{seed}.fq.gz"
    shell: "rasusa reads -c {wildcards.cov} -g {GS} -s {wildcards.seed} -o {output} {input}"

rule asm_short:
    input: r1="sub/short_{cov}x_R1.fq.gz", r2="sub/short_{cov}x_R2.fq.gz"
    output: "asm/short_{cov}x/contigs.fasta"
    threads: T
    shell: "spades.py --isolate -1 {input.r1} -2 {input.r2} -t {threads} -m 10 -o asm/short_{wildcards.cov}x"

rule asm_long:
    input: "sub/long_{cov}x_seed" + str(SEED) + ".fq.gz"
    output: "asm/long_{cov}x/assembly.fasta"
    threads: T
    shell: "flye {FLYE} {input} --out-dir asm/long_{wildcards.cov}x --threads {threads}"

rule asm_long_rep:
    input: "sub/long_{cov}x_seed{seed}.fq.gz"
    output: "asm/long_{cov}x-rep{seed}/assembly.fasta"
    threads: T
    shell: "flye {FLYE} {input} --out-dir asm/long_{wildcards.cov}x-rep{wildcards.seed} --threads {threads}"

rule asm_hybrid_spades:
    input: r1=HSR1, r2=HSR2, l=HLR
    output: HYB + "/contigs.fasta"
    threads: T
    shell: "spades.py -1 {input.r1} -2 {input.r2} --nanopore {input.l} -t {threads} -m 10 -o " + HYB

rule polish:
    input: draft="asm/long_" + str(HL) + "x/assembly.fasta", r1=HSR1, r2=HSR2
    output: POL + "/polished.fasta"
    threads: T
    params: d=POL
    shell:
        """
        cp {input.draft} {params.d}/draft.fasta && bwa index {params.d}/draft.fasta && bwa mem -t {threads} -a {params.d}/draft.fasta {input.r1} > {params.d}/a1.sam && bwa mem -t {threads} -a {params.d}/draft.fasta {input.r2} > {params.d}/a2.sam && polypolish filter --in1 {params.d}/a1.sam --in2 {params.d}/a2.sam --out1 {params.d}/f1.sam --out2 {params.d}/f2.sam && polypolish polish {params.d}/draft.fasta {params.d}/f1.sam {params.d}/f2.sam > {output} && rm -f {params.d}/*.sam
        """

rule quast:
    input: ref="data/ref/ref.fasta", asms=ASSEMBLIES
    output: "results/quast/transposed_report.tsv"
    threads: T
    params: labels=lambda wc, input: ",".join(p.split("/")[1] for p in input.asms)
    shell: "quast.py -r {input.ref} -t {threads} --min-contig 500 --no-icarus -o results/quast {input.asms} --labels {params.labels}"

rule analysis:
    input: ASSEMBLIES, "data/ref/ref.fasta"
    output: "results/indel_context.csv", "results/missing_regions.csv"
    shell: "python scripts/analysis.py"

rule figures:
    input: "results/indel_context.csv", "results/missing_regions.csv", "results/quast/transposed_report.tsv"
    output: "figures/fig1_assembly_metrics.png"
    shell: "python scripts/make_figures.py"
