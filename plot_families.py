import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import PATH_TO_DF
from matplotlib.ticker import PercentFormatter

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})


def prepare_plot_data(df):
    df_biochemical = df[df["assay"] == "biochem"]
    df_cell = df[df["assay"] == "cell"]
    # Create a dictionary to hold the data for the DataFrame
    plot_data = {"Family": df["family"].unique()}
    plot_data["All Sources"] = [
        df[df["family"] == family].shape[0] / len(df) for family in plot_data["Family"]
    ]
    plot_data["Biochemical"] = [
        df_biochemical[df_biochemical["family"] == family].shape[0]
        / len(df_biochemical)
        for family in plot_data["Family"]
    ]
    plot_data["Cell"] = [
        df_cell[df_cell["family"] == family].shape[0] / len(df_cell)
        for family in plot_data["Family"]
    ]
    return pd.DataFrame(plot_data)


if __name__ == "__main__":
    df = pd.read_parquet(PATH_TO_DF)
    df_plot = prepare_plot_data(df)
    df_plot = df_plot.set_index("Family")
    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(6, 4))  # Keeping your desired figsize
    df_plot.plot(kind="bar", ax=ax, width=0.8)

    # Add labels and title
    ax.set_title("Distribution of Protein Families in SAIR Dataset", fontsize=10)
    ax.set_xlabel("Protein Family")
    ax.set_ylabel("Frequency")

    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))

    # Adjust tick parameters for readability on a smaller plot
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)

    plt.legend(title="Analysis Type")
    plt.tight_layout()

    plt.savefig("./figs/protein_families.png", bbox_inches="tight", dpi=600)
