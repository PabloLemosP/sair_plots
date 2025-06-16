import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.stats import spearmanr


pocket_clusters_file = '/scratch/buckets/sandboxaq-sno-scratch-dev/maarten/combined_pocket_lddt_vs_iptm.csv'

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 8})


combined_df = pd.read_csv(pocket_clusters_file)

spearman_corr = spearmanr(combined_df['lddt'], combined_df['iptm'])
print(f'Spearman correlation between lddt and iptm: {spearman_corr}')

plt.hist(combined_df['lddt'], bins=50, edgecolor='black')
plt.title("Distribution of pocket lddt scores")
plt.xlabel("Pocket lddt score")
plt.ylabel("Frequency")
plt.savefig("figs/pocket_lddt_distribution.png", bbox_inches="tight", dpi=600)
plt.close()

spearman_corr = spearmanr(combined_df['confidence_score'], combined_df['qtmscore'])
print(f'Spearman correlation between confidence_score and qtmscore: {spearman_corr}')
