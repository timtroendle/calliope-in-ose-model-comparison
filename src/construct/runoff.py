"""Create run of river capacity factors from Swiss generation time series."""
from pathlib import Path
import zipfile

import pandas as pd


def run_of_river_capacity_factors(path_to_generation, path_to_result):
    """Create run of river capacity factors from Swiss generation time series."""
    data = pd.read_csv(path_to_generation, index_col=0)
    max_power = data["MaxPower_kW"]
    data.drop(
        columns=["Nr", "XCoord_LV03", "YCoord_LV03", "Altitude_m",
                 "WatershedArea_km2", "AnnualProduction_MWh", "MaxPower_kW"],
        inplace=True
    )
    data = data.sum(axis="index") / max_power.sum(axis="index")
    data.index = pd.date_range(start="2015-01-01 00:00", end="2015-12-31 23:00", freq="1H")
    data.name = "capacityfactor"
    data.to_csv(
        path_to_result,
        index=True,
        header=True
    )


if __name__ == "__main__":
    with zipfile.ZipFile(snakemake.input.zip[0], 'r') as zipped:
        zipped.extractall("./build")
    path_to_data = Path(".") / "build" / "opsd-runofriver_production-2017-10-16" / "RunOfRiverPlants_PowerProd_1h_2015.csv"
    assert path_to_data.exists()
    run_of_river_capacity_factors(
        path_to_generation=path_to_data,
        path_to_result=snakemake.output[0]
    )
