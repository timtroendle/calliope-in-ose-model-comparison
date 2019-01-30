"""Compare installed capacities from the model runs."""
import pandas as pd


def compare_installed_capacities(path_to_baseline, path_to_low_cost, path_to_output):
    """Compare installed capacities from the model runs."""
    baseline = pd.read_csv(path_to_baseline, index_col=0)
    low_cost = pd.read_csv(path_to_low_cost, index_col=0)

    (low_cost - baseline).loc["EUR", :].to_csv(
        path_to_output,
        index=True,
        header=["installed capacity diff to baseline"],
        float_format="%.1f"
    )


if __name__ == "__main__":
    compare_installed_capacities(
        path_to_baseline=snakemake.input.baseline,
        path_to_low_cost=snakemake.input.low_cost,
        path_to_output=snakemake.output[0]
    )
