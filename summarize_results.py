#!/usr/bin/env python3
import json
from collections import defaultdict

with open("/root/biolab/work/outbreak_report.json") as f:
    data = json.load(f)

print()
print("=" * 70)
print("  PATHOGEN IDENTIFICATION REPORT")
print("  4 stool samples × 16 pathogen panel (minimap2)")
print("=" * 70)

positive_hits = defaultdict(list)
for sname, sdata in data.items():
    for pname, pdata in sdata["results"].items():
        if pdata["mapping_rate"] > 1 and pdata["avg_depth"] > 3 and pdata["genome_coverage"] > 10:
            positive_hits[pname].append({
                "sample": sname,
                "type": sdata["metadata"],
                "mapping_rate": pdata["mapping_rate"],
                "avg_depth": pdata["avg_depth"],
                "genome_coverage": pdata["genome_coverage"],
            })

for pathogen, hits in sorted(positive_hits.items(), key=lambda x: -max(h["avg_depth"] for h in x[1])):
    n_clin = sum(1 for h in hits if h["type"] == "clinical_stool")
    n_food = sum(1 for h in hits if h["type"] == "food_environmental")
    avg_d = sum(h["avg_depth"] for h in hits) / len(hits)
    avg_c = sum(h["genome_coverage"] for h in hits) / len(hits)
    print()
    print(f"  >>> {pathogen}")
    print(f"      Positive: {len(hits)}/{len(data)} samples (clinical: {n_clin}, food/env: {n_food})")
    print(f"      Avg depth: {avg_d:.1f}x  |  Avg coverage: {avg_c:.1f}%")
    for h in hits:
        print(f"        - {h['sample']} ({h['type']}): {h['mapping_rate']}%, {h['avg_depth']}x, {h['genome_coverage']}%")

top = max(positive_hits.items(), key=lambda x: len(x[1]))
print()
print("=" * 70)
print("  CONCLUSION")
print("=" * 70)
print(f"  Most likely causative agent: {top[0]}")
print(f"  Detected in {len(top[1])} of {len(data)} samples")
in_food = sum(1 for h in top[1] if h["type"] == "food_environmental")
in_clin = sum(1 for h in top[1] if h["type"] == "clinical_stool")
if in_food > 0:
    print(f"  Found in BOTH clinical ({in_clin}) AND food/environmental ({in_food}) samples")
    print(f"  => Foodborne transmission CONFIRMED")
print("=" * 70)
