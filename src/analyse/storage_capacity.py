"""Excavate installed storage capacities from the model results."""
import calliope


def excavate_installed_storage_capacities(path_to_results, scaling_factor, path_to_capacities_raw,
                                          path_to_capacities_publish):
    """Excavate installed storage capacities from the model results."""
    model = calliope.read_netcdf(path_to_results)

    capacities = _read_capacities(model)

    capacities = capacities.groupby(["locs", "techs"]).sum().reset_index()
    capacities = (capacities.pivot(index="locs", columns="techs", values="storage_cap")
                            .div(scaling_factor * 1e3)) # from MWh to GWh
    capacities.loc["EUR"] = capacities.sum()
    capacities.to_csv(
        path_to_capacities_raw,
        index=True,
        header=True,
    )
    capacities.to_csv(
        path_to_capacities_publish,
        index=True,
        header=True,
        float_format="%.1f"
    )


def _read_capacities(model):
    return model.get_formatted_array("storage_cap").to_dataframe().reset_index()


if __name__ == "__main__":
    excavate_installed_storage_capacities(
        path_to_results=snakemake.input.results[0],
        path_to_capacities_raw=snakemake.output.raw,
        path_to_capacities_publish=snakemake.output.publish,
        scaling_factor=snakemake.params.scaling_factor
    )
