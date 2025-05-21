import gzip
import matplotlib.pyplot as plt
import seaborn as sns

def parse_fastq_lengths_and_qualities(fastq_path):
    lengths = []
    qualities = []

    with gzip.open(fastq_path, "rt") as f:
        while True:
            try:
                _ = next(f)  # read name
                seq = next(f).strip()  # sequence
                _ = next(f)  # '+'
                qual = next(f).strip()  # quality
                lengths.append(len(seq))
                avg_qual = sum(ord(c) - 33 for c in qual) / len(qual)
                qualities.append(avg_qual)
            except StopIteration:
                break

    return lengths, qualities

# Input files
files = {
    "bc2041": "bc2041.fastq.gz",
    "bc2042": "bc2042.fastq.gz"
}

for label, filepath in files.items():
    print(f"Processing {label}...")

    lengths, qualities = parse_fastq_lengths_and_qualities(filepath)

    # Create joint plot
    sns.set(style="white")
    g = sns.jointplot(
        x=lengths,
        y=qualities,
        kind="scatter",
        s=5,
        alpha=0.5,
        color="steelblue",
        marginal_kws=dict(bins=100, fill=True),
    )
    g.set_axis_labels("Read length (bp)", "Average quality score")
    g.fig.suptitle(f"Read Length vs Quality - {label}", y=1.03)
    plt.tight_layout()

    # Save the figure
    g.savefig(f"{label}_read_quality_plot.png", dpi=300)
    print(f"Saved plot: {label}_read_quality_plot.png")
