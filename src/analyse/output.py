"""Collect all results in standardised format."""
import functools
from pathlib import Path

import calliope
import pandas as pd


TOTAL_COVERED_REGION = "total covered region"
GERMANY = "DEU"
LOCATIONS = [TOTAL_COVERED_REGION, GERMANY]
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly",
            "wind_onshore_competing", "wind_offshore", "hydro_run_of_river", "biomass"]


def excavate_all_results(paths_to_scenarios, scaling_factors, path_to_output):
    """Collect all results in standardised format"""
    scenarios = {
        Path(path_to_scenario).parent.stem: calliope.read_netcdf(path_to_scenario)
        for path_to_scenario in paths_to_scenarios
    }
    variables = set_up_variables()
    scaling_factors = prepare_scaling_factors(scaling_factors)
    excavate_data(scenarios, variables, scaling_factors).to_csv(
        path_to_output,
        header=True,
        index=True
    )


def excavate_data(scenarios, variables, scaling_factors):
    index = pd.MultiIndex.from_product(
        [scenarios.keys(), LOCATIONS, variables.keys()],
        names=["scenario", "region", "variable"]
    )
    data = pd.Series(index=index, name="value")
    for variable_name, variable_function in variables.items():
        if variable_function:
            excavate_variable(data, scenarios, variable_name, variable_function, scaling_factors)
    return data


def excavate_variable(data, scenarios, variable_name, variable_function, scaling_factors):
    values = {
        scenario_name: variable_function(scenario_data, scaling_factors)
        for scenario_name, scenario_data in scenarios.items()
    }
    for scenario, scenario_values in values.items():
        data.loc[(scenario, GERMANY, variable_name)] = scenario_values.loc[GERMANY]
        data.loc[(scenario, TOTAL_COVERED_REGION, variable_name)] = scenario_values.sum()
        if variable_name == "Energy|Electricity|Renewable curtailment|Relative":
            data.loc[(scenario, TOTAL_COVERED_REGION, variable_name)] = scenario_values.mean()


def prepare_scaling_factors(scaling_factors):
    scaling_factors["energy"] = scaling_factors["power"] / 1e-6 # from MWh to TWh
    scaling_factors["stored_energy"] = scaling_factors["power"] / 1e-3 # from MWh to GWh
    scaling_factors["specific_cost"] = scaling_factors["monetary"] / scaling_factors["power"]
    scaling_factors["power"] = scaling_factors["power"] / 1e-3 # from MW to GW
    return scaling_factors


def set_up_variables():
    return {
        "Cost|Total system": excavate_cost_total_system,
        "Cost|Capacity": excavate_cost_capacity,
        "Cost|Average total system": excavate_levelised_total_cost,
        "Cost|Variable": excavate_cost_variable,
        "Cost|Network": None,
        "Capacity|Electricity|Nuclear": lambda data, sf: excavate_capacity_for_tech(data, ["nuclear"], sf),
        "Capacity|Electricity|Lignite": lambda data, sf: excavate_capacity_for_tech(data, ["lignite"], sf),
        "Capacity|Electricity|Hard coal": lambda data, sf: excavate_capacity_for_tech(data, ["coal"], sf),
        "Capacity|Electricity|Gas CCGT": lambda data, sf: excavate_capacity_for_tech(data, ["ccgt"], sf),
        "Capacity|Electricity|Gas OCGT": None,
        "Capacity|Electricity|Oil": None,
        "Capacity|Electricity|Other nonrenewable": None,
        "Capacity|Electricity|Wind onshore": lambda data, sf: excavate_capacity_for_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf),
        "Capacity|Electricity|Wind offshore": lambda data, sf: excavate_capacity_for_tech(data, ["wind_offshore"], sf),
        "Capacity|Electricity|Solar PV": lambda data, sf: excavate_capacity_for_tech(data, ["open_field_pv", "roof_mounted_pv"], sf),
        "Capacity|Electricity|Solar CSP": None,
        "Capacity|Electricity|Hydro ROR": lambda data, sf: excavate_capacity_for_tech(data, ["hydro_run_of_river"], sf),
        "Capacity|Electricity|Hydro Reservoir": None,
        "Capacity|Electricity|Bioenergy": lambda data, sf: excavate_capacity_for_tech(data, ["biomass"], sf),
        "Capacity|Electricity|Other renewable": None,
        "Capacity|Electricity|Storage|Liion|Power": lambda data, sf: excavate_capacity_for_tech(data, ["battery"], sf),
        "Capacity|Electricity|Storage|Liion|Energy": lambda data, sf: excavate_storage_capacity_for_tech(data, ["battery"], sf),
        "Capacity|Electricity|Storage|Pumped hydro|Power": lambda data, sf: excavate_capacity_for_tech(data, ["pumped_hydro"], sf),
        "Capacity|Electricity|Storage|Pumped hydro|Energy": lambda data, sf: excavate_storage_capacity_for_tech(data, ["pumped_hydro"], sf),
        "Capacity|Electricity|Storage|Other|Power": None,
        "Capacity|Electricity|Storage|Other|Energy": None,
        "Capacity|Electricity|Gross import": excavate_transmission_capacity,
        "Capacity|Electricity|Gross export": excavate_transmission_capacity,
        "Capacity|Sector coupling|Electric vehicles": None,
        "Capacity|Sector coupling|Powertoheat": None,
        "Capacity|Sector coupling|Electrolysis": None,
        "Energy|Electricity|Nuclear": lambda data, sf: excavate_generation_for_tech(data, ["nuclear"], sf),
        "Energy|Electricity|Lignite": lambda data, sf: excavate_generation_for_tech(data, ["lignite"], sf),
        "Energy|Electricity|Hard coal": lambda data, sf: excavate_generation_for_tech(data, ["coal"], sf),
        "Energy|Electricity|Gas CCGT": lambda data, sf: excavate_generation_for_tech(data, ["ccgt"], sf),
        "Energy|Electricity|Gas OCGT": None,
        "Energy|Electricity|Oil": None,
        "Energy|Electricity|Other nonrenewable": None,
        "Energy|Electricity|Wind onshore": lambda data, sf: excavate_generation_for_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf),
        "Energy|Electricity|Wind offshore": lambda data, sf: excavate_generation_for_tech(data, ["wind_offshore"], sf),
        "Energy|Electricity|Solar PV": lambda data, sf: excavate_generation_for_tech(data, ["open_field_pv", "roof_mounted_pv"], sf),
        "Energy|Electricity|Solar CSP": None,
        "Energy|Electricity|Hydro ROR": lambda data, sf: excavate_generation_for_tech(data, ["hydro_run_of_river"], sf),
        "Energy|Electricity|Hydro Reservoir": None,
        "Energy|Electricity|Bioenergy": lambda data, sf: excavate_generation_for_tech(data, ["biomass"], sf),
        "Energy|Electricity|Other renewable": None,
        "Energy|Electricity|Renewable curtailment|Absolute": None, # TODO add
        "Energy|Electricity|Renewable curtailment|Relative": excavate_relative_curtailment,
        "Energy|Electricity|Storage|Liion": lambda data, sf: excavate_generation_for_tech(data, ["battery"], sf),
        "Energy|Electricity|Storage|Pumped hydro": lambda data, sf: excavate_generation_for_tech(data, ["pumped_hydro"], sf),
        "Energy|Electricity|Storage|Other": None,
        "Energy|Electricity|Gross import": excavate_transmission_generation,
        "Energy|Electricity|Gross export": excavate_transmission_consumption,
        "Energy|Sector coupling|Electric vehicles|G2V": None,
        "Energy|Sector coupling|Powertoheat": None,
        "Energy|Sector coupling|Electrolysis": None,
        "Cycles|Electricity|Storage|Liion": None,
        "Cycles|Electricity|Storage|Pumped hydro": None,
        "Cycles|Electricity|Storage|Other": None,
        "Carbon emissions": excavate_co2_total_system,
        "Price|Electricity|Weighted average": None,
        "Price|Electricity|Storage|Liion|Weighted average charging": None,
        "Price|Electricity|Storage|Liion|Weighted average discharging": None,
        "Startups|Electricity|Total number": None,
        "Startups|Electricity|Total cost": None,
        "Energy|Electricity|Peak demand": excavate_peak_demand,
        "Energy|Electricity|Total demand": excavate_demand,
    }


def excavate_cost_total_system(data, scaling_factors):
    return (data.get_formatted_array("cost")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost
                .sum()
                .loc["monetary"]
                .div(scaling_factors["monetary"]))


def excavate_cost_capacity(data, scaling_factors):
    return (data.get_formatted_array("cost_investment")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost_investment
                .sum()
                .loc["monetary"]
                .div(scaling_factors["monetary"]))


def excavate_levelised_total_cost(data, scaling_factors):
    lcoe = (data.get_formatted_array("total_levelised_cost")
                .to_dataframe()
                .loc[("electricity", "monetary"), "total_levelised_cost"]) / scaling_factors["specific_cost"]
    return pd.Series([lcoe], index=[GERMANY])


def excavate_cost_variable(data, scaling_factors):
    return (data.get_formatted_array("cost_var")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost_var
                .sum()
                .loc["monetary"]
                .div(scaling_factors["monetary"]))


def excavate_capacity_for_tech(data, techs, scaling_factors):
    return (excavate_capacity(data, scaling_factors["power"]).loc[techs]
                                                             .reset_index()
                                                             .groupby("locs")
                                                             .energy_cap
                                                             .sum())


def excavate_transmission_capacity(data, scaling_factors):
    capacity = excavate_capacity(data, scaling_factors["power"]).reset_index()
    capacity = (capacity.drop(index=capacity[~capacity.techs.str.contains("ac_transmission")].index)
                        .reset_index()
                        .groupby("locs")
                        .energy_cap
                        .sum())
    if len(capacity.index) == 0: # happens for Germany only scenarios
        return pd.Series([0], index=[GERMANY])
    else:
        return capacity


@functools.lru_cache(maxsize=10, typed=False)
def excavate_capacity(data, scaling_factor):
    return (data.get_formatted_array("energy_cap")
                .to_dataframe()
                .reset_index()
                .set_index(["techs", "locs"])
                .div(scaling_factor))


def excavate_storage_capacity_for_tech(data, techs, scaling_factors):
    return (excavate_storage_capacity(data, scaling_factors["stored_energy"]).loc[techs]
                                                                             .reset_index()
                                                                             .groupby("locs")
                                                                             .storage_cap
                                                                             .sum())


@functools.lru_cache(maxsize=10, typed=False)
def excavate_storage_capacity(data, scaling_factor):
    return (data.get_formatted_array("storage_cap")
                .to_dataframe()
                .reset_index()
                .set_index(["techs", "locs"])
                .div(scaling_factor))


def excavate_generation_for_tech(data, techs, scaling_factors):
    return (excavate_generation(data, scaling_factors["energy"]).loc[techs]
                                                                .reset_index()
                                                                .groupby("locs")
                                                                .carrier_prod
                                                                .sum())


def excavate_transmission_generation(data, scaling_factors):
    generation = excavate_generation(data, scaling_factors["energy"]).reset_index()
    generation = (generation.drop(index=generation[~generation.techs.str.contains("ac_transmission")].index)
                            .reset_index()
                            .groupby("locs")
                            .carrier_prod
                            .sum())
    if len(generation.index) == 0: # happens for Germany only scenarios
        return pd.Series([0], index=[GERMANY])
    else:
        return generation


@functools.lru_cache(maxsize=10, typed=False)
def excavate_generation(data, scaling_factor):
    return (data.get_formatted_array("carrier_prod")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .carrier_prod
                .sum()
                .div(scaling_factor))


@functools.lru_cache(maxsize=10, typed=False)
def excavate_resource(data):
    return (data.get_formatted_array("resource")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .resource
                .mean())


@functools.lru_cache(maxsize=10, typed=False)
def excavate_capacity_factor(data):
    return (data.get_formatted_array("capacity_factor")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .capacity_factor
                .mean())


def excavate_relative_curtailment(data, scaling_factors):
    resource = (excavate_resource(data).loc[RE_TECHS]
                                       .reset_index()
                                       .groupby("locs")
                                       .resource
                                       .mean())
    capacity_factor = (excavate_capacity_factor(data).loc[RE_TECHS]
                                                     .reset_index()
                                                     .groupby("locs")
                                                     .capacity_factor
                                                     .mean())
    return capacity_factor / resource


def excavate_transmission_consumption(data, scaling_factors):
    consumption = excavate_consumption(data, scaling_factors["energy"]).reset_index()
    consumption = (consumption.drop(index=consumption[~consumption.techs.str.contains("ac_transmission")].index)
                              .reset_index()
                              .groupby("locs")
                              .carrier_con
                              .sum())
    if len(consumption.index) == 0: # happens for Germany only scenarios
        return pd.Series([0], index=[GERMANY])
    else:
        return consumption


def excavate_demand(data, scaling_factors):
    consumption = excavate_consumption(data, scaling_factors["energy"]).reset_index()
    consumption = (consumption.drop(index=consumption[consumption.techs.str.contains("ac_transmission")].index)
                              .reset_index()
                              .groupby("locs")
                              .carrier_con
                              .sum())
    return consumption


@functools.lru_cache(maxsize=10, typed=False)
def excavate_consumption(data, scaling_factor):
    return (data.get_formatted_array("carrier_con")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .carrier_con
                .sum()
                .mul(-1)
                .div(scaling_factor))


def excavate_peak_demand(data, scaling_factors):
    return (data.get_formatted_array("carrier_con")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .carrier_con
                .max()
                .loc["demand_elec"]
                .mul(-1)
                .div(scaling_factors["power"]))


def excavate_co2_total_system(data, scaling_factors):
    return (data.get_formatted_array("cost")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost
                .sum()
                .loc["co2"]
                .div(scaling_factors["co2"]))


if __name__ == "__main__":
    excavate_all_results(
        paths_to_scenarios=snakemake.input.scenarios,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
