import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams.update({"font.size": 11})

PATH_TO_DF = "/scratch/buckets/sb-alg-dgx-2q25/pablo/final_dfs/sair_v0.parquet"


def analyze_failures(df, analysis_name):
    print(f"{analysis_name} Failure Analysis:")
    print("=" * 30)

    df_failed = df[df["all_passed"] == False]

    pb_columns = [
        "mol_pred_loaded",
        "sanitization",
        "inchi_convertible",
        "all_atoms_connected",
        "bond_lengths",
        "bond_angles",
        "internal_steric_clash",
        "aromatic_ring_flatness",
        "double_bond_flatness",
        "internal_energy",
        "mol_cond_loaded",
        "passes_valence_checks",
        "passes_kekulization",
        "number_clashes",
        "number_short_outlier_bonds",
        "number_long_outlier_bonds",
        "number_valid_bonds",
        "number_valid_angles",
        "number_valid_noncov_pairs",
        "number_aromatic_rings",
        "number_double_bonds",
    ]

    len_df_pre_nan = len(df_failed)
    df_failed_cleaned = df_failed.dropna(
        subset=pb_columns
    )  # Only consider pb_columns for dropping NaNs
    num_nan_rows_removed = len_df_pre_nan - len(df_failed_cleaned)
    nan_failure_rate = (
        (num_nan_rows_removed / len_df_pre_nan * 100) if len_df_pre_nan > 0 else 0.0
    )

    print(
        f"Number of NaN rows removed: {num_nan_rows_removed}. Failure rate: {nan_failure_rate:.2f}%"
    )

    if len(df_failed_cleaned) == 0:
        print("No failures after NaN removal for this analysis.")
        # Return a dictionary with 0.0 for all relevant failure types
        return {
            pb: 0.0
            for pb in pb_columns
            if pb
            not in [
                "mol_pred_loaded",
                "sanitization",
                "inchi_convertible",
                "all_atoms_connected",
                "mol_cond_loaded",
                "passes_valence_checks",
                "passes_kekulization",
            ]
        }

    fail_counts = (df_failed_cleaned[pb_columns] == False).sum(axis=0)

    denominator = len(df_failed_cleaned)
    if denominator == 0:
        failure_rates = {pb: 0.0 for pb in pb_columns}
    else:
        failure_rates = {
            pb: (count / denominator) * 100 for pb, count in fail_counts.items()
        }

    # print("Failures per PoseBuster:")
    # for pb, rate in failure_rates.items():
    #     print(f"{pb}: {int(fail_counts[pb])} failures. Failure rate: {rate:.2f}%")

    # For these, there is no failure rate, so we dont show them
    filtered_rates = {
        k: v
        for k, v in failure_rates.items()
        if k
        not in [
            "mol_pred_loaded",
            "sanitization",
            "inchi_convertible",
            "all_atoms_connected",
            "mol_cond_loaded",
            "passes_valence_checks",
            "passes_kekulization",
        ]
    }
    return filtered_rates


def prepare_plot_data(df):
    df_biochemical = df[df["assay"] == "biochem"]
    df_cell = df[df["assay"] == "cell"]
    df_homogenate = df[df["assay"] == "homogenate"]

    # Collect failure rates
    overall_failures = analyze_failures(df, "Overall")
    biochemical_failures = analyze_failures(df_biochemical, "Biochemical Assays")
    cell_failures = analyze_failures(df_cell, "Cell Assays")

    all_failure_types = sorted(
        list(
            set(
                list(overall_failures.keys())
                + list(biochemical_failures.keys())
                + list(cell_failures.keys())
            )
        )
    )

    # Create a dictionary to hold the data for the DataFrame
    plot_data = {"Failure Type": all_failure_types}
    plot_data["Overall"] = [overall_failures.get(ft, 0.0) for ft in all_failure_types]
    plot_data["Biochemical"] = [
        biochemical_failures.get(ft, 0.0) for ft in all_failure_types
    ]
    plot_data["Cell"] = [cell_failures.get(ft, 0.0) for ft in all_failure_types]
    return pd.DataFrame(plot_data)


if __name__ == "__main__":
    df_plot = prepare_plot_data(pd.read_parquet(PATH_TO_DF))
    df_plot = df_plot.set_index("Failure Type")

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(6, 4))  # Keeping your desired figsize
    df_plot.plot(kind="bar", ax=ax, width=0.8)

    # Add labels and title
    ax.set_title("Failure Rates per PoseBuster Category")
    ax.set_xlabel("PoseBuster Failure Type")
    ax.set_ylabel("Failure Rate (\%)")

    # Adjust tick parameters for readability on a smaller plot
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)

    plt.legend(title="Analysis Type")
    plt.tight_layout()
    plt.savefig("./figs/posebusters_failure_rates.png", bbox_inches="tight", dpi=600)
