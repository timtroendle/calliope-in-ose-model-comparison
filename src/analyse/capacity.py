"""Excavate installed capacities from the model results."""
import calliope
import pandas as pd


def excavate_installed_capacities(path_to_results, path_to_capacities):
    """Excavate installed capacities from the model results."""
    model = calliope.read_netcdf(path_to_results)

    capacities = _read_capacities(model)

    capacities.loc[:, "location"] = capacities.loc[:, "location"].str[:3]
    capacities.dropna(subset=["energy_cap"], inplace=True)
    capacities = capacities.groupby(["location", "tech"]).sum().reset_index()
    capacities = capacities.pivot(index="location", columns="tech", values="energy_cap") / 1e6 # from kW to GW
    del capacities["demand_elec"]
    capacities.loc["EUR"] = capacities.sum()
    capacities.to_csv(
        path_to_capacities,
        index=True,
        header=True,
        float_format="%.1f"
    )


def _read_capacities(model):
    data = model.results["energy_cap"].to_dataframe().reset_index()

    data["location"] = data["loc_techs"].map(lambda x: x.split("::")[0])
    data["tech"] = data["loc_techs"].map(lambda x: x.split("::")[1])

    data.loc[data.tech.map(lambda tech: "transmission" in tech), "tech"] = pd.np.nan
    return data.dropna(subset=["tech"])


if __name__ == "__main__":
    excavate_installed_capacities(
        path_to_results=snakemake.input.results[0],
        path_to_capacities=snakemake.output[0]
    )
