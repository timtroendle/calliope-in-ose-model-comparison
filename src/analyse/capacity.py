"""Excavate installed capacities from the model results."""
import calliope
import pandas as pd


def excavate_installed_capacities(path_to_results, path_to_capacities):
    """Excavate installed capacities from the model results."""
    model = calliope.read_netcdf(path_to_results)

    capacities = _read_capacities(model)

    capacities.loc[:, "locs"] = capacities.loc[:, "locs"].str[:3]
    capacities.dropna(subset=["energy_cap"], inplace=True)
    capacities = capacities.groupby(["locs", "techs"]).sum().reset_index()
    capacities = capacities.pivot(index="locs", columns="techs", values="energy_cap") / 1e6 # from kW to GW
    del capacities["demand_elec"]
    capacities.loc["EUR"] = capacities.sum()
    capacities.to_csv(
        path_to_capacities,
        index=True,
        header=True,
        float_format="%.1f"
    )


def _read_capacities(model):
    data = model.get_formatted_array("energy_cap").to_dataframe().reset_index()
    data.loc[data.techs.map(lambda tech: "transmission" in tech), "techs"] = pd.np.nan
    return data.dropna(subset=["techs"])


if __name__ == "__main__":
    excavate_installed_capacities(
        path_to_results=snakemake.input.results[0],
        path_to_capacities=snakemake.output[0]
    )
