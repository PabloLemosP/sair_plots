import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})

PATH_TO_DF = "/scratch/buckets/sb-alg-dgx-2q25/pablo/final_dfs/sair_v0.parquet"


def analyse_df_and_return_scores(df_subset):
    """
    Analyzes a DataFrame subset and returns a dictionary of calculated scores.
    """
    # Initialize a dictionary to store scores for each method
    method_scores = {}

    # Define methods and their respective score columns and whether to negate them
    methods_info = {
        "Vina": {"score_col": "vina_score", "negate": True},
        "Vina minimized": {"score_col": "vina_score_min", "negate": True},
        "Vinardo": {"score_col": "vinardo_score", "negate": True},
        "OnionNet": {"score_col": "onionnet_score", "negate": False},
        "AevPlig": {"score_col": "aevplig_score", "negate": False},
        "iPTM": {"score_col": "iptm", "negate": False},
    }

    for method_name, info in methods_info.items():
        scores_raw = df_subset[info["score_col"]].values
        nan_mask = np.isnan(scores_raw)
        scores = scores_raw[~nan_mask]
        if info["negate"]:
            scores = -scores  # Apply negation if required

        potency = df_subset["potency"].values[~nan_mask]

        # Ensure there's enough data to calculate correlations and AUC
        if len(scores) > 1 and len(potency) > 1 and len(np.unique(potency > 7)) > 1:
            spearman_corr = spearmanr(scores, potency).correlation
            pearson_corr = pearsonr(scores, potency)[0]
            kendall_corr = kendalltau(scores, potency).correlation
            auc = roc_auc_score(potency > 7, scores)

            method_scores[method_name] = {
                "Spearman": spearman_corr,
                "Pearson": pearson_corr,
                "Kendall": kendall_corr,
                "AUC": auc,
            }
        else:
            # Handle cases with insufficient data by assigning NaN or 0
            method_scores[method_name] = {
                "Spearman": np.nan,
                "Pearson": np.nan,
                "Kendall": np.nan,
                "AUC": np.nan,
            }

    return method_scores


if __name__ == "__main__":
    df = pd.read_parquet(PATH_TO_DF)
    df = df[df["source"] == "ChEMBL"]  # Filter for ChEMBL data
    df_highconf = df[df["confidence_score"] > 0.8]  # Filter for high confidence scores
    # Prepare data for plotting
    all_dfs = []
    all_dfs_highconf = []

    # Overall Scores
    overall_scores_dict = analyse_df_and_return_scores(df)
    for method, scores in overall_scores_dict.items():
        all_dfs.append({"Assay Type": "All Sources", "Method": method, **scores})

    # High Confidence Scores
    overall_highconf_scores_dict = analyse_df_and_return_scores(df_highconf)
    for method, scores in overall_highconf_scores_dict.items():
        all_dfs_highconf.append(
            {"Assay Type": "All Sources", "Method": method, **scores}
        )

    # Biochemical Assay Scores
    df_biochem = df[df["assay"] == "biochem"]
    biochem_scores_dict = analyse_df_and_return_scores(df_biochem)
    for method, scores in biochem_scores_dict.items():
        all_dfs.append({"Assay Type": "Biochemical", "Method": method, **scores})

    # High Confidence Biochemical Assay Scores
    df_biochem_highconf = df_highconf[df_highconf["assay"] == "biochem"]
    biochem_highconf_scores_dict = analyse_df_and_return_scores(df_biochem_highconf)
    for method, scores in biochem_highconf_scores_dict.items():
        all_dfs_highconf.append(
            {"Assay Type": "Biochemical", "Method": method, **scores}
        )
    # High Confidence Cell Assay Scores

    # Cell Assay Scores
    df_cell = df[df["assay"] == "cell"]
    cell_scores_dict = analyse_df_and_return_scores(df_cell)
    for method, scores in cell_scores_dict.items():
        all_dfs.append({"Assay Type": "Cell", "Method": method, **scores})

    # High Confidence Cell Assay Scores
    df_cell_highconf = df_highconf[df_highconf["assay"] == "cell"]
    cell_highconf_scores_dict = analyse_df_and_return_scores(df_cell_highconf)
    for method, scores in cell_highconf_scores_dict.items():
        all_dfs_highconf.append({"Assay Type": "Cell", "Method": method, **scores})

    # Create the final DataFrame for plotting
    df_plot = pd.DataFrame(all_dfs)
    df_plot_highconf = pd.DataFrame(all_dfs_highconf)

    # Set up the plots
    metrics = ["Spearman", "Pearson", "Kendall", "AUC"]
    methods = df_plot["Method"].unique()
    assay_types = df_plot["Assay Type"].unique()

    n_metrics = len(metrics)
    n_methods = len(methods)
    n_assay_types = len(assay_types)

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(6, 6), sharey=False)
    plot_titles = [
        "Spearman Correlation",
        "Pearson Correlation",
        "Kendall Correlation",
        "AUC",
    ]
    axes = axes.flatten()  # Flatten axes array for easy iteration

    bar_width = 0.15  # Adjust bar width for better spacing
    group_spacing = 0.8  # Space between groups of methods

    for i, metric in enumerate(metrics):
        ax = axes[i]  # Get the current subplot axis
        metric_df = df_plot[["Assay Type", "Method", metric]]
        metric_df_highconf = df_plot_highconf[["Assay Type", "Method", metric]]

        # Calculate x-positions for each group of bars (per method)
        indices = np.arange(n_methods) * (n_assay_types * bar_width + group_spacing)
        colors = ["tab:blue", "tab:orange", "tab:green"]

        for j, assay_type in enumerate(assay_types):
            # Extract scores for the current metric and assay type for all methods
            # Use .loc to ensure correct indexing and avoid SettingWithCopyWarning
            scores_for_plot = [
                metric_df.loc[
                    (metric_df["Method"] == method)
                    & (metric_df["Assay Type"] == assay_type),
                    metric,
                ].values[0]
                for method in methods
            ]
            ax.bar(
                indices + j * bar_width,
                scores_for_plot,
                bar_width,
                label=assay_type,
                color=colors[j],
                # alpha=0.5,
            )

            # Add high confidence scores as with alpha blending
            scores_for_plot_highconf = [
                metric_df_highconf.loc[
                    (metric_df_highconf["Method"] == method)
                    & (metric_df_highconf["Assay Type"] == assay_type),
                    metric,
                ].values[0]
                for method in methods
            ]
            if assay_type == "All Sources":
                ax.bar(
                    indices + j * bar_width,
                    scores_for_plot_highconf,
                    bar_width,
                    alpha=0.5,
                    color=colors[j],
                    label="High Conf",
                )
            else:
                ax.bar(
                    indices + j * bar_width,
                    scores_for_plot_highconf,
                    bar_width,
                    alpha=0.5,
                    color=colors[j],
                )

        # ax.set_xlabel("Method")
        ax.set_ylabel(plot_titles[i])
        # ax.set_title(f'Scores for {metric} Metric')
        # Center x-ticks under the method groups
        ax.set_xticks(
            indices + (n_assay_types - 1) * bar_width / 2,
            labels=methods,
            rotation=45,
            ha="right",
            fontsize=8,
        )
        # ax.set_xticklabels(methods, fontsize=9)
        # ax.legend(title='Assay Type')
        ax.grid(axis="y", linestyle="--", alpha=0.7)
        ax.axhline(
            0, color="grey", linewidth=0.8
        )  # Add a horizontal line at 0 for reference

    # axes[0].axhline(
    #     0, color="black", linewidth=0.8, ls="--"
    # )  # Add a horizontal line at 0 for reference
    # axes[1].axhline(
    #     0, color="black", linewidth=0.8, ls="--"
    # )  # Add a horizontal line at 0 for reference
    # axes[2].axhline(
    #     0, color="black", linewidth=0.8, ls="--"
    # )  # Add a horizontal line at 0 for reference
    # axes[3].axhline(
    #     0.5, color="black", linewidth=0.8, ls="--"
    # )  # Add a horizontal line at 0 for reference
    axes[3].set_ylim(0.5, 0.75)

    axes[3].legend(fontsize=8, loc="upper left")

    plt.tight_layout()  # Adjust layout to prevent labels from overlapping
    plt.savefig(
        "./figs/scores_by_method_and_assay_type.png", dpi=600, bbox_inches="tight"
    )
