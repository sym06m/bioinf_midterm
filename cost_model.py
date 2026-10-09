"""Marginal sequencing cost of each assembly (USD) from the bases it used.  Usage: python scripts/cost_model.py

Model: cost = Gb used x price per Gb. Library preparation, instrument purchase, labour and compute are NOT included.
A whole flow cell / kit is bought at once, so this per-Gb model only holds when a run is shared by many samples."""
import re
import pandas as pd

GS = 4_641_652
PRICE = {
    "short": {"usd_per_gb": round(1608 / 8.5, 2),
              "source": "Univ. of Kentucky Genomics Core Illumina handout 2022-23: MiSeq v2 500-cycle (2x250) USD 1,608, max 8.5 Gb. "
                        "https://www.ukhealthcare.uky.edu/sites/default/files/2022-11/OG%20SRF_NGS_Illumina%20Sequencing%202022-2023%20Handout.pdf"},
    "long": {"usd_per_gb": round(700 / 30, 2),
             "source": "ONT store list price, MinION flow cell R10.4.1 (FLO-MIN114) USD 700, https://store.nanoporetech.com/priceList.html ; "
                       "up to 30 Gb per flow cell (FEMS Microbiol Ecol 2021;97(3):fiab001, PMC8068755)"},
}
LOW_YIELD_LONG = round(700 / 10, 2)   # sensitivity case: flow cell yields only 10 Gb (an assumption, not from a source)

def gb_used(label):
    m = re.match(r"(short|long)_(\d+)x", label)
    if m:
        return {m.group(1): int(m.group(2)) * GS / 1e9}
    m = re.match(r"(hybridspades|polish)_s(\d+)_l(\d+)", label)
    return {"short": int(m.group(2)) * GS / 1e9, "long": int(m.group(3)) * GS / 1e9}

q = pd.read_csv("results/assembly_summary_v2.csv")
q = q[q["label"] != "long_30x"].copy()   # failed run
rows = []
for _, r in q.iterrows():
    g = gb_used(r["label"])
    base = sum(PRICE[k]["usd_per_gb"] * v for k, v in g.items())
    low = sum((LOW_YIELD_LONG if k == "long" else PRICE[k]["usd_per_gb"]) * v for k, v in g.items())
    rows.append({"label": r["label"], "gb_short": round(g.get("short", 0), 3), "gb_long": round(g.get("long", 0), 3),
                 "cost_usd": round(base, 2), "cost_usd_if_long_10Gb": round(low, 2),
                 "N50": r["N50"], "NGA50": r["NGA50"], "misassemblies": r["# misassemblies"],
                 "genome_fraction": r["Genome fraction (%)"],
                 "mismatches_100kb": r["# mismatches per 100 kbp"], "indels_100kb": r["# indels per 100 kbp"]})
out = pd.DataFrame(rows)
out.to_csv("results/cost_vs_quality.csv", index=False)
print(f"short {PRICE['short']['usd_per_gb']} USD/Gb | long {PRICE['long']['usd_per_gb']} USD/Gb")
print(out.to_string(index=False))
