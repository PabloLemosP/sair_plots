import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 8})

PATH_TO_DF = "/scratch/buckets/sb-alg-dgx-2q25/pablo/final_dfs/sair_v0.parquet"

metrics = [
    "ptm",
    "iptm",
    "complex_ipde",
    "complex_pde",
    "complex_iplddt",
    "complex_plddt",
    "confidence_score",
    "chains_ptm",
    "interaction_ptm",
]


def get_spearman_results(dfs):
    spearman_results = []
    for metric in metrics:
        corrs = []
        for df in dfs:
            # The omit NaN part of spearmanr is broken, so we need to handle NaNs manually
            metric_values = df[metric].values
            potency_values = df["potency"].values
            metric_nans = np.isnan(metric_values)
            metric_values = metric_values[~metric_nans]
            potency_values = potency_values[~metric_nans]

            corr, _ = spearmanr(metric_values, potency_values)
            # print(
            #     f"Spearman correlation for {metric} in {df['assay'].iloc[0]}: {corr:.4f}"
            # )
            corrs.append(corr)
        spearman_results.append(corrs)
    return spearman_results


if __name__ == "__main__":
    df = pd.read_parquet(PATH_TO_DF)
    assert df["source"].unique().tolist() == [
        "ChEMBL",
        "BindingDB",
    ], "Unexpected sources in the dataframe"

    # Drop duplicates based on 'entry_id' keeping ChEMBL entries first, as they are the ones with assay info
    priority = {"ChEMBL": 0, "BindingDB": 1}
    df["source_priority"] = df["source"].map(priority)
    df = (
        df.sort_values("source_priority")
        .drop_duplicates(["entry_id", "index"], keep="first")
        .drop(columns="source_priority")
    )

    df_biochem = df[df["assay"] == "biochem"]
    df_cell = df[df["assay"] == "cell"]
    # df_homogenate = df[df["assay"] == "homogenate"]

    spearman_results = get_spearman_results([df, df_biochem, df_cell])
    assay_labels = ["All Sources", "Bioch", "Cell"]
    title_labels = [
        "PTM (+)",
        "iPTM (+)",
        "Complex iPDE (-)",
        "Complex PDE (-)",
        "Complex iPLDDT (+)",
        "Complex pLDDT (+)",
        "Confidence Score (+)",
        "Chains PTM (+)",
        "Interaction PTM (+)",
    ]

    fig, axes = plt.subplots(3, 3, figsize=(6, 5), sharey=True)
    for i, (ax, metric, title) in enumerate(zip(axes.flat, metrics, title_labels)):
        ax.bar(
            assay_labels,
            spearman_results[i],
            color=["tab:blue", "tab:orange", "tab:green", "tab:red"],
        )
        ax.set_title(title, fontsize=9)
        ax.set_ylim(-0.4, 0.4)
        if i % 3 == 0:
            ax.set_ylabel("Spearman Correlation (IC 50)", fontsize=9)
        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        for j, v in enumerate(spearman_results[i]):
            ax.text(
                j,
                v,
                f"{v:.2f}",
                ha="center",
                va="bottom" if v >= 0 else "top",
                fontsize=7,
                color="black",
            )
    plt.tight_layout()
    plt.savefig("./figs/spearman_correlations.png", dpi=600, bbox_inches="tight")
