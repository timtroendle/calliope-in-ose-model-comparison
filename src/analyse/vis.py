"""This module contains the visualisation of results."""
import matplotlib.pyplot as plt
import seaborn as sns
import calliope


def visualise_model_results(path_to_results, path_to_figure):
    """Plot the results."""
    model = calliope.read_netcdf(path_to_results)

    generation = _create_generation_timeseries(model) / 1e6 # from kW to GW
    consumption = _create_consumption_timeseries(model) / 1e6 # from kW to GW

    generation = generation.reindex(columns=generation.columns.union(consumption.columns), fill_value=0)
    consumption = consumption.reindex(columns=generation.columns.union(consumption.columns), fill_value=0)

    sns.set_context('paper')
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    ax.plot(generation + consumption)
    ax.legend(generation.columns)
    ax.set_xlabel("time")
    ax.set_ylabel("electricity generation [GW]")
    fig.savefig(path_to_figure, dpi=300)


def _create_generation_timeseries(model):
    ts = model.get_formatted_array("carrier_prod").to_dataframe(
        name="carrier_prod"
    ).reset_index().groupby(["techs", "timesteps"]).sum().reset_index().pivot(
        index="timesteps",
        columns="techs",
        values="carrier_prod"
    ).rename(
        columns=lambda tech: tech if tech[-4] != "_" else tech[:-4]
    )
    return ts.groupby(
        ts.columns,
        axis='columns'
    ).sum().drop(columns=[column for column in ts.columns if "transmission" in column])


def _create_consumption_timeseries(model):
    ts = model.get_formatted_array("carrier_con").to_dataframe(
        name="carrier_con"
    ).reset_index().groupby(["techs", "timesteps"]).sum().reset_index().pivot(
        index="timesteps",
        columns="techs",
        values="carrier_con"
    ).rename(
        columns=lambda tech: tech if tech[-4] != "_" else tech[:-4]
    )
    return ts.groupby(
        ts.columns,
        axis='columns'
    ).sum().drop(columns=[column for column in ts.columns if "transmission" in column])


if __name__ == "__main__":
    visualise_model_results(
        path_to_results=snakemake.input.results[0],
        path_to_figure=snakemake.output[0]
    )
