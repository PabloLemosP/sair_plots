import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

pocket_clusters_file = "/scratch/buckets/sandboxaq-sno-scratch-dev/maarten/combined_pocket_lddt_vs_iptm.csv"

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})


combined_df = pd.read_csv(pocket_clusters_file)

spearman_corr = spearmanr(combined_df["lddt"], combined_df["iptm"])
print(f"Spearman correlation between lddt and iptm: {spearman_corr}")
fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(
    combined_df["lddt"],
    bins=50,
    edgecolor="black",
    weights=np.ones_like(combined_df["lddt"]) / len(combined_df["lddt"]),
)
ax.set_title("Distribution of pocket LDDT scores")
ax.set_xlabel("Pocket LDDT score")
ax.set_ylabel("Frequency")
from matplotlib.ticker import PercentFormatter

plt.gca().yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
plt.savefig("figs/pocket_lddt_distribution.png", bbox_inches="tight", dpi=600)
plt.close()

spearman_corr = spearmanr(combined_df["confidence_score"], combined_df["qtmscore"])
print(f"Spearman correlation between confidence_score and qtmscore: {spearman_corr}")
