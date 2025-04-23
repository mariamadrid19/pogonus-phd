import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import numpy as np

# === Load contig lengths from FASTA ===
fasta_file = "Pogonus_gilv.asm.bp.p_ctg.fa"
# Change the name of the fasta file accordingly
lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

# === Plot 1: Contig Length Histogram ===
plt.figure(figsize=(10, 6))
sns.histplot(lengths, bins=100, log_scale=(True, True), color="skyblue")
plt.xlabel("Contig Length (bp, log scale)")
plt.ylabel("Count (log scale)")
plt.title("Contig Length Distribution")
plt.tight_layout()
plt.savefig("contig_length_histogram.png")
plt.close()

# === Plot 2: Cumulative Assembly Size ===
sorted_lengths = sorted(lengths, reverse=True)
cum_lengths = np.cumsum(sorted_lengths)
total_length = cum_lengths[-1]

plt.figure(figsize=(10, 6))
plt.plot(cum_lengths / 1e6, label="Cumulative length")
plt.axhline(y=total_length * 0.5 / 1e6, color='red', linestyle='--', label="N50 threshold")
plt.xlabel("Number of Contigs")
plt.ylabel("Cumulative Length (Mb)")
plt.title("Cumulative Assembly Size")
plt.legend()
plt.tight_layout()
plt.savefig("cumulative_assembly_plot.png")
plt.close()

# === Plot 3: Nx Statistics Bar Plot ===
# Adjust these values accordingly 
nx_values = {
    "N50": 12774182,
    "N60": 9968333,
    "N70": 5987516,
    "N80": 4244452,
    "N90": 363732,
    "N100": 1482
}

plt.figure(figsize=(8, 5))
plt.bar(nx_values.keys(), [v / 1e6 for v in nx_values.values()], color="teal")
plt.ylabel("Length (Mb)")
plt.title("Nx Statistics for Assembly")
plt.tight_layout()
plt.savefig("nx_bar_plot.png")
plt.close()

print("Plots saved as PNGs:")
print(" - contig_length_histogram.png")
print(" - cumulative_assembly_plot.png")
print(" - nx_bar_plot.png")
