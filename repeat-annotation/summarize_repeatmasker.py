#!/usr/bin/env python3

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fasta_lengths(path):
    lengths = {}
    current = None

    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())

    return lengths


def parse_repeatmasker_out(path):
    records = []

    with open_text(path) as handle:
        for line in handle:
            fields = line.split()

            # Skip the two header lines and malformed records.
            if len(fields) < 14 or not fields[0].isdigit():
                continue

            try:
                divergence = float(fields[1])
                deletions = float(fields[2])
                insertions = float(fields[3])
                sequence = fields[4]
                start = int(fields[5])
                end = int(fields[6])
                repeat_name = fields[9]
                class_family = fields[10]
            except (ValueError, IndexError):
                continue

            length = end - start + 1

            if "/" in class_family:
                repeat_class, repeat_family = class_family.split("/", 1)
            else:
                repeat_class = class_family
                repeat_family = "Unclassified"

            records.append({
                "sequence": sequence,
                "start": start,
                "end": end,
                "length": length,
                "divergence": divergence,
                "deletions": deletions,
                "insertions": insertions,
                "repeat_name": repeat_name,
                "class_family": class_family,
                "repeat_class": repeat_class,
                "repeat_family": repeat_family,
            })

    return records


def merge_intervals(records):
    """Calculate non-overlapping masked bases for each sequence."""
    intervals = defaultdict(list)

    for record in records:
        intervals[record["sequence"]].append(
            (record["start"], record["end"])
        )

    masked = {}

    for sequence, values in intervals.items():
        values.sort()
        total = 0
        current_start, current_end = values[0]

        for start, end in values[1:]:
            if start <= current_end + 1:
                current_end = max(current_end, end)
            else:
                total += current_end - current_start + 1
                current_start, current_end = start, end

        total += current_end - current_start + 1
        masked[sequence] = total

    return masked


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize RepeatMasker .out results."
    )
    parser.add_argument(
        "-o", "--repeat-out",
        required=True,
        help="RepeatMasker .out or .out.gz file"
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="Original genome FASTA or FASTA.gz"
    )
    parser.add_argument(
        "-p", "--prefix",
        default="repeat_summary",
        help="Output prefix (default: repeat_summary)"
    )
    parser.add_argument(
        "--top",
        type=int,
        default=25,
        help="Number of top repeat names to report (default: 25)"
    )
    args = parser.parse_args()

    records = parse_repeatmasker_out(args.repeat_out)
    genome_lengths = read_fasta_lengths(args.fasta)

    if not records:
        raise SystemExit(
            "No RepeatMasker records were found. Check the input .out file."
        )

    genome_size = sum(genome_lengths.values())
    merged_masked = merge_intervals(records)
    masked_bases = sum(merged_masked.values())

    class_stats = defaultdict(lambda: {"count": 0, "bases": 0})
    class_family_stats = defaultdict(lambda: {"count": 0, "bases": 0})
    name_stats = defaultdict(lambda: {"count": 0, "bases": 0})
    scaffold_stats = defaultdict(lambda: {"count": 0, "bases": 0})
    divergence_stats = defaultdict(lambda: {"count": 0, "bases": 0})

    for record in records:
        length = record["length"]

        class_stats[record["repeat_class"]]["count"] += 1
        class_stats[record["repeat_class"]]["bases"] += length

        class_family_stats[record["class_family"]]["count"] += 1
        class_family_stats[record["class_family"]]["bases"] += length

        name_stats[record["repeat_name"]]["count"] += 1
        name_stats[record["repeat_name"]]["bases"] += length

        scaffold_stats[record["sequence"]]["count"] += 1
        scaffold_stats[record["sequence"]]["bases"] += length

        lower = int(record["divergence"] // 2) * 2
        label = f"{lower:02d}-{lower + 2:02d}"
        divergence_stats[label]["count"] += 1
        divergence_stats[label]["bases"] += length

    # Overall human-readable report
    with open(f"{args.prefix}.txt", "w") as out:
        out.write("RepeatMasker summary\n")
        out.write("====================\n\n")
        out.write(f"Genome size:              {genome_size:,} bp\n")
        out.write(f"Repeat annotations:       {len(records):,}\n")
        out.write(f"Sum of repeat fragments:  {sum(r['length'] for r in records):,} bp\n")
        out.write(f"Non-overlapping masked:   {masked_bases:,} bp\n")
        out.write(
            f"Genome masked:            "
            f"{100 * masked_bases / genome_size:.2f}%\n\n"
        )

        out.write("Repeat classes\n")
        out.write("--------------------\n")
        out.write(f"{'Class':30s} {'Count':>12s} {'Bases':>15s} {'Genome %':>10s}\n")

        for name, stats in sorted(
            class_stats.items(),
            key=lambda item: item[1]["bases"],
            reverse=True
        ):
            percentage = 100 * stats["bases"] / genome_size
            out.write(
                f"{name:30s} {stats['count']:12,d} "
                f"{stats['bases']:15,d} {percentage:9.3f}%\n"
            )

    def write_stats(path, stats, first_column):
        with open(path, "w") as out:
            out.write(f"{first_column}\tcount\tbases\tgenome_percent\n")
            for name, values in sorted(
                stats.items(),
                key=lambda item: item[1]["bases"],
                reverse=True
            ):
                percentage = 100 * values["bases"] / genome_size
                out.write(
                    f"{name}\t{values['count']}\t{values['bases']}\t"
                    f"{percentage:.6f}\n"
                )

    write_stats(
        f"{args.prefix}.by_class.tsv",
        class_stats,
        "repeat_class"
    )
    write_stats(
        f"{args.prefix}.by_class_family.tsv",
        class_family_stats,
        "class_family"
    )
    write_stats(
        f"{args.prefix}.by_repeat_name.tsv",
        name_stats,
        "repeat_name"
    )

    # Per-sequence table uses merged intervals to avoid double counting.
    with open(f"{args.prefix}.by_sequence.tsv", "w") as out:
        out.write(
            "sequence\tsequence_length\trepeat_hits\t"
            "masked_bases\tmasked_percent\n"
        )

        for sequence, sequence_length in sorted(
            genome_lengths.items(),
            key=lambda item: item[1],
            reverse=True
        ):
            count = scaffold_stats[sequence]["count"]
            bases = merged_masked.get(sequence, 0)
            percentage = (
                100 * bases / sequence_length if sequence_length else 0
            )
            out.write(
                f"{sequence}\t{sequence_length}\t{count}\t"
                f"{bases}\t{percentage:.6f}\n"
            )

    with open(f"{args.prefix}.divergence.tsv", "w") as out:
        out.write("divergence_bin\tcount\tbases\tgenome_percent\n")

        for label in sorted(
            divergence_stats,
            key=lambda value: int(value.split("-")[0])
        ):
            stats = divergence_stats[label]
            percentage = 100 * stats["bases"] / genome_size
            out.write(
                f"{label}\t{stats['count']}\t{stats['bases']}\t"
                f"{percentage:.6f}\n"
            )

    print(f"Parsed {len(records):,} RepeatMasker records")
    print(f"Genome size: {genome_size:,} bp")
    print(f"Non-overlapping masked bases: {masked_bases:,} bp")
    print(f"Genome masked: {100 * masked_bases / genome_size:.2f}%")
    print(f"Results written with prefix: {args.prefix}")


if __name__ == "__main__":
    main()
