"""Excavate levelised costs from model runs."""
from pathlib import Path

import calliope
import pandas as pd


def compare_cost(paths_to_results, path_to_output):
    """Excavate levelised costs from model runs."""
    costs = {
        Path(path_to_result).parent.stem: _read_cost(path_to_result)
        for path_to_result in paths_to_results
    }

    pd.Series(costs).rename(index="levelised cost of electricity").to_csv(
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
        paths_to_results=snakemake.input.results,
        path_to_output=snakemake.output[0]
    )
