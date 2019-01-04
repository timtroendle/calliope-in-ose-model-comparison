"""This module contains the visualisation of results."""
import matplotlib.pyplot as plt
import seaborn as sns
import calliope


def visualise_model_results(path_to_results, path_to_figure):
    """Plot the results."""
    model = calliope.read_netcdf(path_to_results)

    generation = _create_generation_timeseries(model)
    consumption = _create_consumption_timeseries(model)

    consumption = consumption.reindex(columns=generation.columns, fill_value=0)

    sns.set_context('paper')
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    ax.plot(generation + consumption)
    ax.legend(generation.columns)
    ax.set_xlabel("time")
    ax.set_ylabel("electricity generation")
    fig.savefig(path_to_figure, dpi=300)


def _create_generation_timeseries(model):
    data = model.results["carrier_prod"].to_dataframe().reset_index()

    data["location"] = data["loc_tech_carriers_prod"].map(lambda x: x.split("::")[0])
    data["tech"] = data["loc_tech_carriers_prod"].map(lambda x: x.split("::")[1])
    data["carrier"] = data["loc_tech_carriers_prod"].map(lambda x: x.split("::")[2])

    ts = data.groupby(["tech", "timesteps"]).sum().reset_index().pivot(
        index="timesteps",
        columns="tech",
        values="carrier_prod"
    )
    return ts.drop(columns=[column for column in ts.columns if "transmission" in column])


def _create_consumption_timeseries(model):
    data = model.results["carrier_con"].to_dataframe().reset_index()

    data["location"] = data["loc_tech_carriers_con"].map(lambda x: x.split("::")[0])
    data["tech"] = data["loc_tech_carriers_con"].map(lambda x: x.split("::")[1])
    data["carrier"] = data["loc_tech_carriers_con"].map(lambda x: x.split("::")[2])

    ts = data.groupby(["tech", "timesteps"]).sum().reset_index().pivot(
        index="timesteps",
        columns="tech",
        values="carrier_con"
    )
    return ts.drop(columns=[column for column in ts.columns if "transmission" in column])


if __name__ == "__main__":
    visualise_model_results(
        path_to_results=snakemake.input.results[0],
        path_to_figure=snakemake.output[0]
    )
