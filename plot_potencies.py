import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})

PATH_TO_DF = "/scratch/buckets/sb-alg-dgx-2q25/pablo/final_dfs/sair_v0.parquet"

if __name__ == "__main__":
    df = pd.read_parquet(PATH_TO_DF)
    assert df["source"].unique().tolist() == [
        "ChEMBL",
        "BindingDB",
    ], "Unexpected sources in the dataframe"

    df_chembl = df[df["source"] == "ChEMBL"]
    df_bindingdb = df[df["source"] == "BindingDB"]

    df_cell = df[df["assay"] == "cell"]
    df_biochem = df[df["assay"] == "biochem"]
    df_homogenate = df[df["assay"] == "homogenate"]

    fig, axis = plt.subplots(ncols=3, nrows=2, figsize=(6, 4), sharex=True, sharey=True)
    axis[0, 0].hist(
        df["pIC50"],
        bins=20,
        alpha=0.5,
        label="All",
        color="blue",
        density=True,
    )
    axis[0, 1].hist(
        df_chembl["pIC50"],
        bins=20,
        alpha=0.5,
        label="ChEMBL",
        color="red",
        density=True,
    )
    axis[0, 2].hist(
        df_bindingdb["pIC50"],
        bins=20,
        alpha=0.5,
        label="BindingDB",
        color="grey",
        density=True,
    )
    # axis[0, 0].set_xlabel(r'$- \log_{10} \left( \mathrm{IC50 [nM]} \right) $')
    axis[0, 0].set_ylabel("Frequency")
    axis[0, 0].set_title("All Sources")
    axis[0, 1].set_title("ChEMBL")
    axis[0, 2].set_title("BindingDB")
    # axis[0, 1].set_xlabel(r'$- \log_{10} \left( \mathrm{IC50 [nM]} \right) $')
    # axis[0, 2].set_xlabel(r'$- \log_{10} \left( \mathrm{IC50 [nM]} \right) $')
    # axis[0].legend()
    axis[0, 0].set_xlim(3.5, 12)
    axis[0, 1].set_xlim(3.5, 12)
    axis[0, 2].set_xlim(3.5, 12)

    axis[1, 0].hist(
        df_cell["pIC50"], bins=20, alpha=0.5, label="Cell", color="green", density=True
    )
    axis[1, 1].hist(
        df_biochem["pIC50"],
        bins=20,
        alpha=0.5,
        label="Biochem",
        color="orange",
        density=True,
    )
    axis[1, 2].hist(
        df_homogenate["pIC50"],
        bins=20,
        alpha=0.5,
        label="Homogenate",
        color="pink",
        density=True,
    )
    axis[1, 0].set_xlabel(r"$- \log_{10} \left( \mathrm{IC50 [nM]} \right) $")
    axis[1, 0].set_ylabel("Frequency")
    axis[1, 0].set_title("Cell")
    axis[1, 1].set_title("Biochemical")
    axis[1, 2].set_title("Homogenate")
    axis[1, 1].set_xlabel(r"$- \log_{10} \left( \mathrm{IC50 [nM]} \right) $")
    axis[1, 2].set_xlabel(r"$- \log_{10} \left( \mathrm{IC50 [nM]} \right) $")
    axis[1, 0].set_xlim(3.5, 12)
    axis[1, 1].set_xlim(3.5, 12)
    axis[1, 2].set_xlim(3.5, 12)

    plt.tight_layout()
    plt.savefig(
        "./figs/ic50_distribution_by_source_and_assay.png",
        dpi=600,
        bbox_inches="tight",
    )
