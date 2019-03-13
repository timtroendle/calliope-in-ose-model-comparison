"""Compare installed storage capacities from the model runs."""
import pandas as pd


def compare_installed_storage_capacities(path_to_baseline, path_to_other, path_to_output):
    """Compare installed cstorage apacities from the model runs."""
    baseline = pd.read_csv(path_to_baseline, index_col=0)
    other = pd.read_csv(path_to_other, index_col=0)

    (other - baseline).loc["DEU", :].to_csv(
        path_to_output,
        index=True,
        header=["installed storage capacity diff in DEU to baseline"],
        float_format="%.1f"
    )


if __name__ == "__main__":
    compare_installed_storage_capacities(
        path_to_baseline=snakemake.input.baseline,
        path_to_other=snakemake.input.other,
        path_to_output=snakemake.output[0]
    )
