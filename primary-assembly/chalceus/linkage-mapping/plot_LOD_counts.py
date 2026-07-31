#!/usr/bin/env python3

from collections import Counter
from pathlib import Path
import csv
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

LOD_VALUES = [10, 11, 12, 15]
MAP_PATTERN = "map{lod}.txt"

# Plot the expected major linkage groups
MAX_LG = 12

OUTPUT_PLOT = "lod_marker_counts.png"
OUTPUT_TABLE = "lod_marker_counts.tsv"


# ------------------------------------------------------------
# Read a Lep-MAP3 map file
# ------------------------------------------------------------

def count_markers(map_file):
    """
    Count markers assigned to each linkage group.

    Lep-MAP3 SeparateChromosomes2 output has one linkage-group
    assignment per non-comment line. LG 0 represents unassigned
    markers/singletons.
    """
    counts = Counter()

    with open(map_file, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            first_field = line.split()[0]

            try:
                lg = int(first_field)
            except ValueError:
                raise ValueError(
                    f"Cannot interpret line {line_number} of "
                    f"{map_file}: {line!r}"
                )

            counts[lg] += 1

    return counts


# ------------------------------------------------------------
# Load all LOD results
# ------------------------------------------------------------

all_counts = {}

for lod in LOD_VALUES:
    map_file = Path(MAP_PATTERN.format(lod=lod))

    if not map_file.is_file():
        raise FileNotFoundError(f"Missing input file: {map_file}")

    all_counts[lod] = count_markers(map_file)

    assigned = sum(
        count for lg, count in all_counts[lod].items() if lg != 0
    )
    unassigned = all_counts[lod].get(0, 0)
    number_of_lgs = sum(
        1 for lg, count in all_counts[lod].items()
        if lg != 0 and count > 0
    )

    print(
        f"LOD {lod}: "
        f"{number_of_lgs} LGs, "
        f"{assigned:,} assigned markers, "
        f"{unassigned:,} unassigned markers"
    )


# ------------------------------------------------------------
# Write count table
# ------------------------------------------------------------

with open(OUTPUT_TABLE, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow(
        ["LOD", "LG", "marker_count"]
    )

    for lod in LOD_VALUES:
        for lg in sorted(all_counts[lod]):
            writer.writerow(
                [lod, lg, all_counts[lod][lg]]
            )


# ------------------------------------------------------------
# Plot major linkage groups
# ------------------------------------------------------------

linkage_groups = list(range(1, MAX_LG + 1))

plt.figure(figsize=(14, 8))

for lod in LOD_VALUES:
    marker_counts = [
        all_counts[lod].get(lg, 0)
        for lg in linkage_groups
    ]

    plt.plot(
        linkage_groups,
        marker_counts,
        marker="o",
        linewidth=2,
        markersize=7,
        label=f"LOD {lod}",
    )

plt.xlabel("Linkage group (LG)", fontsize=13)
plt.ylabel("Number of markers", fontsize=13)
plt.title(
    f"LG1–LG{MAX_LG} marker counts at different LOD limits",
    fontsize=15,
)
plt.xticks(linkage_groups)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(OUTPUT_PLOT, dpi=300)
plt.close()

print()
print(f"Plot written to: {OUTPUT_PLOT}")
print(f"Counts written to: {OUTPUT_TABLE}")
