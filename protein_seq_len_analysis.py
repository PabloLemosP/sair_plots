import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})

# Read the data
df = pd.read_csv("data/correct_unique_sequences.csv")
df["seq_len"] = df["input_receptor"].str.len()
df.to_csv("data/with_count_unique_seqs.csv", index=False)

# Frequency of sequence lengths
freq = df["seq_len"].value_counts().sort_index()
freq_df = freq.reset_index()
freq_df.columns = ["seq_len", "frequency"]
freq_df.to_csv("data/seq_len_frequency.csv", index=False)

bins = 50
hist, bin_edges = np.histogram(df["seq_len"], bins=bins)
# Convert frequency to percentage
hist_percentage = hist / hist.sum()  # * 100

dark_blue = "#4882b4"

plt.figure(figsize=(6, 4))
for i in range(len(hist)):
    plt.bar(
        (bin_edges[i] + bin_edges[i + 1]) / 2,
        hist_percentage[i],
        width=bin_edges[i + 1] - bin_edges[i],
        color=dark_blue,
        align="center",
        edgecolor="black",  # Add black outline
        linewidth=1,  # Set outline thickness
    )

plt.gca().yaxis.set_major_formatter(PercentFormatter(1, decimals=0))

plt.xlabel(r"Sequence Length")
plt.ylabel(r"Frequency")

# Set serif font for all tick labels and axes labels
plt.xticks(fontsize=11, fontfamily="serif")
plt.yticks(fontsize=11, fontfamily="serif")

plt.tight_layout()
plt.savefig("figs/seq_len_histogram_percentage.png", dpi=600, bbox_inches="tight")
plt.close()
