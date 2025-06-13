import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pocket_clusters_file = '/scratch/buckets/sandboxaq-sno-scratch-dev/maarten/pocket_clusters.txt'

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 8})

fig, axis = plt.subplots(ncols=2, nrows=1, figsize=(4, 2), sharex=True, sharey=True)


with open(pocket_clusters_file, 'r') as f:
    lines = f.readlines()
pocket_clusters = [line.strip().split('\t') for line in lines]
proteins = [cluster[0] for cluster in pocket_clusters]
counts = [int(cluster[1]) for cluster in pocket_clusters]
# index of the protein with the most clusters
max_index = np.argmax(counts)
print(f'Max clusters: {proteins[max_index]} with {counts[max_index]} clusters')

counts_blub = [c if c <= 10 else 11 for c in counts]
bins = list(range(1, 13))
plt.hist(counts_blub, bins=bins, edgecolor='black', rwidth=0.8)
plt.xticks(ticks=[i+0.5 for i in range(1, 12)], labels=[str(i) for i in range(1, 11)] + ['$>10$'])


# Add labels and title
plt.title("Distinct pockets per protein")
plt.xlabel("Number of pockets")
plt.ylabel("Number of proteins")

plt.tight_layout()
plt.savefig("./figs/pocket_per_protein.png", bbox_inches="tight", dpi=600)
plt.close()

pocket_clusters_file = '/scratch/buckets/sandboxaq-sno-scratch-dev/maarten/subclust_pockets.csv'


df = pd.read_csv(pocket_clusters_file, sep=',')
bins = np.arange(0.5, 5.6, 1)  # bins centered on 1, 2, 3, 4, 5
plt.hist(df['numpockets'], bins=bins, edgecolor='black', rwidth=0.8)
plt.xticks([1, 2, 3, 4, 5])

# Add labels and titlex
plt.title("Distinct pockets per batch")
plt.xlabel("Number of pockets")
plt.ylabel("Number of predictions")

plt.tight_layout()
plt.savefig("./figs/pockets_per_minibatch.png", bbox_inches="tight", dpi=600)
