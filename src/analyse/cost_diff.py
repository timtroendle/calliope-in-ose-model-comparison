"""Excavate cost diff from model runs."""
import calliope
import pandas as pd


def compare_cost(path_to_baseline, path_to_low_cost, path_to_output):
    """Excavate cost diff from model runs."""
    baseline_cost = _read_cost(path_to_baseline)
    low_cost_cost = _read_cost(path_to_low_cost)

    pd.Series({
        "baseline": baseline_cost,
        "low-cost": low_cost_cost,
        "diff to baseline": low_cost_cost - baseline_cost
    }).rename(index="levelised cost of electricity").to_csv(
        path_to_output,
        index=True,
        header=True,
        float_format="%.1f"
    )


def _read_cost(path_to_model):
    model = calliope.read_netcdf(path_to_model)
    return (model.get_formatted_array("total_levelised_cost")
                 .to_dataframe()
                 .loc[("electricity", "monetary"), "total_levelised_cost"])


if __name__ == "__main__":
    compare_cost(
        path_to_baseline=snakemake.input.baseline,
        path_to_low_cost=snakemake.input.low_cost,
        path_to_output=snakemake.output[0]
    )
