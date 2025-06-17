import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

pocket_clusters_file = (
    "/scratch/buckets/sandboxaq-sno-scratch-dev/maarten/pocket_clusters.txt"
)

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})

fig, axis = plt.subplots(ncols=2, nrows=1, figsize=(6, 3), sharey=False)


with open(pocket_clusters_file, "r") as f:
    lines = f.readlines()
pocket_clusters = [line.strip().split("\t") for line in lines]
proteins = [cluster[0] for cluster in pocket_clusters]
counts = [int(cluster[1]) for cluster in pocket_clusters]
# index of the protein with the most clusters
max_index = np.argmax(counts)
print(f"Max clusters: {proteins[max_index]} with {counts[max_index]} clusters")

counts_blub = [c if c <= 10 else 11 for c in counts]
bins = list(range(1, 13))
axis[1].hist(counts_blub, bins=bins, edgecolor="black", rwidth=0.8, density=True)
axis[1].set_xticks(
    ticks=[i + 0.5 for i in range(1, 12)],
    labels=[str(i) for i in range(1, 11)] + ["$>10$"],
)

# Set the font size for the x-ticks
axis[1].tick_params(axis="x", labelsize=8)
axis[1].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))


# Add labels and title
axis[1].set_title("Distinct pockets per protein")
axis[1].set_xlabel("Number of pockets")
axis[1].set_ylabel("Fraction of proteins")

# plt.tight_layout()
# plt.savefig("./figs/pocket_per_protein.png", bbox_inches="tight", dpi=600)
# plt.close()

pocket_clusters_file = (
    "/scratch/buckets/sandboxaq-sno-scratch-dev/maarten/subclust_pockets.csv"
)


df = pd.read_csv(pocket_clusters_file, sep=",")
bins = np.arange(0.5, 5.6, 1)  # bins centered on 1, 2, 3, 4, 5
# fig, axis = plt.subplots(ncols=1, nrows=1, figsize=(6, 4))
axis[0].hist(df["numpockets"], bins=bins, edgecolor="black", rwidth=0.8, density=True)
axis[0].set_xticks([1, 2, 3, 4, 5])
axis[0].tick_params(axis="x", labelsize=8)
axis[0].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))

# Add labels and titlex
axis[0].set_title("Distinct pockets per batch")
axis[0].set_xlabel("Number of pockets")
axis[0].set_ylabel("Fraction of predictions")

plt.tight_layout()
# plt.savefig("./figs/pockets_per_minibatch.png", bbox_inches="tight", dpi=600)
plt.savefig(
    "./figs/pockets_per_protein_and_minibatch.png", bbox_inches="tight", dpi=600
)
