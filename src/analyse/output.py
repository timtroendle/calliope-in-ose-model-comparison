"""Collect all results in standardised format."""
import functools
from pathlib import Path

import calliope
import pandas as pd


MODEL_NAME = "euro-calliope"
TOTAL_COVERED_REGION = "Total covered region"
GERMANY_IN = "DEU"
GERMANY_OUT = "DE"
LOCATIONS = [TOTAL_COVERED_REGION, GERMANY_OUT]
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly",
            "wind_onshore_competing", "wind_offshore", "hydro_run_of_river"]
SCENARIO_NAME_MAP = {
    "baseline": "Baseline battery costs|Full geographic coverage|Without sector coupling",
    "low-cost": "50percent battery costs|Full geographic coverage|Without sector coupling",
    "lowest-cost": "25percent battery costs|Full geographic coverage|Without sector coupling",
    "baseline-germany": "Baseline battery costs|Germany only|Without sector coupling",
    "low-cost-germany": "50percent battery costs|Germany only|Without sector coupling",
    "lowest-cost-germany": "25percent battery costs|Germany only|Without sector coupling"
}


def excavate_all_results(paths_to_scenarios, scaling_factors, path_to_output):
    """Collect all results in standardised format"""
    scenarios = {
        Path(path_to_scenario).parent.stem: calliope.read_netcdf(path_to_scenario)
        for path_to_scenario in paths_to_scenarios
    }
    variables = _set_up_variables()
    scaling_factors = _prepare_scaling_factors(scaling_factors)
    data = _excavate_data(scenarios, variables, scaling_factors).reset_index()
    data["model"] = MODEL_NAME
    data["scenario"] = data["scenario"].map(SCENARIO_NAME_MAP)
    data[["model", "scenario", "region", "variable", "value"]].to_csv(
        path_to_output,
        header=True,
        index=False
    )


def _excavate_data(scenarios, variables, scaling_factors):
    index = pd.MultiIndex.from_product(
        [scenarios.keys(), LOCATIONS, variables.keys()],
        names=["scenario", "region", "variable"]
    )
    data = pd.Series(index=index, name="value")
    for variable_name, variable_function in variables.items():
        if variable_function:
            _excavate_variable(data, scenarios, variable_name, variable_function, scaling_factors)
    return data


def _excavate_variable(data, scenarios, variable_name, variable_function, scaling_factors):
    values = {
        scenario_name: variable_function(scenario_data, scaling_factors)
        for scenario_name, scenario_data in scenarios.items()
    }
    for scenario, scenario_values in values.items():
        data.loc[(scenario, GERMANY_OUT, variable_name)] = scenario_values.loc[GERMANY_IN]
        data.loc[(scenario, TOTAL_COVERED_REGION, variable_name)] = scenario_values.sum()
        if variable_name == "Energy|Electricity|Renewable curtailment|Relative":
            data.loc[(scenario, TOTAL_COVERED_REGION, variable_name)] = scenario_values.mean()


def _prepare_scaling_factors(scaling_factors):
    scaling_factors["energy"] = scaling_factors["power"] / 1e-6 # from MWh to TWh
    scaling_factors["stored_energy"] = scaling_factors["power"] / 1e-3 # from MWh to GWh
    scaling_factors["specific_cost"] = scaling_factors["monetary"] / scaling_factors["power"]
    scaling_factors["power"] = scaling_factors["power"] / 1e-3 # from MW to GW
    return scaling_factors


def _set_up_variables():
    return {
        "Cost|Total system": _excavate_cost_total_system,
        "Cost|Capacity": _excavate_cost_capacity,
        "Cost|Average total system": _excavate_levelised_total_cost,
        "Cost|Variable": _excavate_cost_variable,
        "Cost|Network": None,
        "Capacity|Electricity|Nuclear": lambda data, sf: _excavate_capacity_for_tech(data, ["nuclear"], sf),
        "Capacity|Electricity|Lignite": lambda data, sf: _excavate_capacity_for_tech(data, ["lignite"], sf),
        "Capacity|Electricity|Hard coal": lambda data, sf: _excavate_capacity_for_tech(data, ["coal"], sf),
        "Capacity|Electricity|Gas CCGT": lambda data, sf: _excavate_capacity_for_tech(data, ["ccgt"], sf),
        "Capacity|Electricity|Gas OCGT": None,
        "Capacity|Electricity|Oil": None,
        "Capacity|Electricity|Other nonrenewable": None,
        "Capacity|Electricity|Wind onshore": lambda data, sf: _excavate_capacity_for_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf),
        "Capacity|Electricity|Wind offshore": lambda data, sf: _excavate_capacity_for_tech(data, ["wind_offshore"], sf),
        "Capacity|Electricity|Solar PV": lambda data, sf: _excavate_capacity_for_tech(data, ["open_field_pv", "roof_mounted_pv"], sf),
        "Capacity|Electricity|Solar CSP": None,
        "Capacity|Electricity|Hydro ROR": lambda data, sf: _excavate_capacity_for_tech(data, ["hydro_run_of_river"], sf),
        "Capacity|Electricity|Hydro Reservoir": None,
        "Capacity|Electricity|Bioenergy": lambda data, sf: _excavate_capacity_for_tech(data, ["biomass"], sf),
        "Capacity|Electricity|Other renewable": None,
        "Capacity|Electricity|Storage|Liion|Power": lambda data, sf: _excavate_capacity_for_tech(data, ["battery"], sf),
        "Capacity|Electricity|Storage|Liion|Energy": lambda data, sf: _excavate_storage_capacity_for_tech(data, ["battery"], sf),
        "Capacity|Electricity|Storage|Pumped hydro|Power": lambda data, sf: _excavate_capacity_for_tech(data, ["pumped_hydro"], sf),
        "Capacity|Electricity|Storage|Pumped hydro|Energy": lambda data, sf: _excavate_storage_capacity_for_tech(data, ["pumped_hydro"], sf),
        "Capacity|Electricity|Storage|Other|Power": None,
        "Capacity|Electricity|Storage|Other|Energy": None,
        "Capacity|Electricity|Gross import": _excavate_transmission_capacity,
        "Capacity|Electricity|Gross export": _excavate_transmission_capacity,
        "Capacity|Sector coupling|Electric vehicles": None,
        "Capacity|Sector coupling|Powertoheat": None,
        "Capacity|Sector coupling|Electrolysis": None,
        "Energy|Electricity|Nuclear": lambda data, sf: _excavate_generation_for_tech(data, ["nuclear"], sf),
        "Energy|Electricity|Lignite": lambda data, sf: _excavate_generation_for_tech(data, ["lignite"], sf),
        "Energy|Electricity|Hard coal": lambda data, sf: _excavate_generation_for_tech(data, ["coal"], sf),
        "Energy|Electricity|Gas CCGT": lambda data, sf: _excavate_generation_for_tech(data, ["ccgt"], sf),
        "Energy|Electricity|Gas OCGT": None,
        "Energy|Electricity|Oil": None,
        "Energy|Electricity|Other nonrenewable": None,
        "Energy|Electricity|Wind onshore": lambda data, sf: _excavate_generation_for_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf),
        "Energy|Electricity|Wind offshore": lambda data, sf: _excavate_generation_for_tech(data, ["wind_offshore"], sf),
        "Energy|Electricity|Solar PV": lambda data, sf: _excavate_generation_for_tech(data, ["open_field_pv", "roof_mounted_pv"], sf),
        "Energy|Electricity|Solar CSP": None,
        "Energy|Electricity|Hydro ROR": lambda data, sf: _excavate_generation_for_tech(data, ["hydro_run_of_river"], sf),
        "Energy|Electricity|Hydro Reservoir": None,
        "Energy|Electricity|Bioenergy": lambda data, sf: _excavate_generation_for_tech(data, ["biomass"], sf),
        "Energy|Electricity|Other renewable": None,
        "Energy|Electricity|Renewable curtailment|Absolute": _excavate_absolute_curtailment,
        "Energy|Electricity|Renewable curtailment|Relative": _excavate_relative_curtailment,
        "Energy|Electricity|Storage|Liion": lambda data, sf: _excavate_generation_for_tech(data, ["battery"], sf),
        "Energy|Electricity|Storage|Pumped hydro": lambda data, sf: _excavate_generation_for_tech(data, ["pumped_hydro"], sf),
        "Energy|Electricity|Storage|Other": None,
        "Energy|Electricity|Gross import": _excavate_transmission_generation,
        "Energy|Electricity|Gross export": _excavate_transmission_consumption,
        "Energy|Sector coupling|Electric vehicles|G2V": None,
        "Energy|Sector coupling|Powertoheat": None,
        "Energy|Sector coupling|Electrolysis": None,
        "Cycles|Electricity|Storage|Liion": lambda data, sf: _excavate_cycles_for_tech(data, ["battery"], sf),
        "Cycles|Electricity|Storage|Pumped hydro": lambda data, sf: _excavate_cycles_for_tech(data, ["pumped_hydro"], sf),
        "Cycles|Electricity|Storage|Other": None,
        "Carbon emissions": _excavate_co2_total_system,
        "Price|Electricity|Weighted average": None,
        "Price|Electricity|Storage|Liion|Weighted average charging": None,
        "Price|Electricity|Storage|Liion|Weighted average discharging": None,
        "Startups|Electricity|Total number": None,
        "Startups|Electricity|Total cost": None,
        "Energy|Electricity|Peak demand": _excavate_peak_demand,
        "Energy|Electricity|Total demand": _excavate_demand,
    }


def _excavate_cost_total_system(data, scaling_factors):
    return (data.get_formatted_array("cost")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost
                .sum()
                .loc["monetary"]
                .div(scaling_factors["monetary"]))


def _excavate_cost_capacity(data, scaling_factors):
    return (data.get_formatted_array("cost_investment")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost_investment
                .sum()
                .loc["monetary"]
                .div(scaling_factors["monetary"]))


def _excavate_levelised_total_cost(data, scaling_factors):
    lcoe = (data.get_formatted_array("total_levelised_cost")
                .to_dataframe()
                .loc[("electricity", "monetary"), "total_levelised_cost"]) / scaling_factors["specific_cost"]
    return pd.Series([lcoe], index=[GERMANY_IN])


def _excavate_cost_variable(data, scaling_factors):
    return (data.get_formatted_array("cost_var")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost_var
                .sum()
                .loc["monetary"]
                .div(scaling_factors["monetary"]))


def _excavate_capacity_for_tech(data, techs, scaling_factors):
    return (_excavate_capacity(data, scaling_factors["power"]).loc[techs]
                                                              .reset_index()
                                                              .groupby("locs")
                                                              .energy_cap
                                                              .sum())


def _excavate_transmission_capacity(data, scaling_factors):
    capacity = _excavate_capacity(data, scaling_factors["power"]).reset_index()
    capacity = (capacity.drop(index=capacity[~capacity.techs.str.contains("ac_transmission")].index)
                        .reset_index()
                        .groupby("locs")
                        .energy_cap
                        .sum())
    if len(capacity.index) == 0: # happens for Germany only scenarios
        return pd.Series([0], index=[GERMANY_IN])
    else:
        return capacity


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_capacity(data, scaling_factor):
    return (data.get_formatted_array("energy_cap")
                .to_dataframe()
                .reset_index()
                .set_index(["techs", "locs"])
                .div(scaling_factor))


def _excavate_storage_capacity_for_tech(data, techs, scaling_factors):
    return (_excavate_storage_capacity(data, scaling_factors["stored_energy"]).loc[techs]
                                                                              .reset_index()
                                                                              .groupby("locs")
                                                                              .storage_cap
                                                                              .sum())


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_storage_capacity(data, scaling_factor):
    return (data.get_formatted_array("storage_cap")
                .to_dataframe()
                .reset_index()
                .set_index(["techs", "locs"])
                .div(scaling_factor))


def _excavate_generation_for_tech(data, techs, scaling_factors):
    return (_excavate_generation(data, scaling_factors["energy"]).loc[techs]
                                                                 .reset_index()
                                                                 .groupby("locs")
                                                                 .carrier_prod
                                                                 .sum())


def _excavate_transmission_generation(data, scaling_factors):
    generation = _excavate_generation(data, scaling_factors["energy"]).reset_index()
    generation = (generation.drop(index=generation[~generation.techs.str.contains("ac_transmission")].index)
                            .reset_index()
                            .groupby("locs")
                            .carrier_prod
                            .sum())
    if len(generation.index) == 0: # happens for Germany only scenarios
        return pd.Series([0], index=[GERMANY_IN])
    else:
        return generation


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_generation(data, scaling_factor):
    return (data.get_formatted_array("carrier_prod")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .carrier_prod
                .sum()
                .div(scaling_factor))


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_renewable_generation_potential(data, scaling_factor):
    resource = data.get_formatted_array("resource").sel(techs=RE_TECHS)
    capacity = data.get_formatted_array("energy_cap").sel(techs=RE_TECHS) / scaling_factor
    return (resource * capacity).sum(dim=["timesteps", "techs"])


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_renewable_generation(data, scaling_factor):
    return (data.get_formatted_array("carrier_prod")
                .sel(techs=RE_TECHS)
                .sum(dim=["timesteps", "techs"])) / scaling_factor


def _excavate_relative_curtailment(data, scaling_factors):
    potential = _excavate_renewable_generation_potential(data, scaling_factors["energy"])
    generated = _excavate_renewable_generation(data, scaling_factors["energy"])
    return ((potential - generated) / potential).to_pandas().electricity


def _excavate_absolute_curtailment(data, scaling_factors):
    potential = _excavate_renewable_generation_potential(data, scaling_factors["energy"])
    generated = _excavate_renewable_generation(data, scaling_factors["energy"])
    return (potential - generated).to_pandas().electricity


def _excavate_transmission_consumption(data, scaling_factors):
    consumption = _excavate_consumption(data, scaling_factors["energy"]).reset_index()
    consumption = (consumption.drop(index=consumption[~consumption.techs.str.contains("ac_transmission")].index)
                              .reset_index()
                              .groupby("locs")
                              .carrier_con
                              .sum())
    if len(consumption.index) == 0: # happens for Germany only scenarios
        return pd.Series([0], index=[GERMANY_IN])
    else:
        return consumption


def _excavate_demand(data, scaling_factors):
    consumption = _excavate_consumption(data, scaling_factors["energy"]).reset_index()
    consumption = (consumption.drop(index=consumption[consumption.techs.str.contains("ac_transmission")].index)
                              .reset_index()
                              .groupby("locs")
                              .carrier_con
                              .sum())
    return consumption


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_consumption(data, scaling_factor):
    return (data.get_formatted_array("carrier_con")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .carrier_con
                .sum()
                .mul(-1)
                .div(scaling_factor))


def _excavate_peak_demand(data, scaling_factors):
    return (data.get_formatted_array("carrier_con")
                .to_dataframe()
                .reset_index()
                .groupby(["techs", "locs"])
                .carrier_con
                .min()
                .loc["demand_elec"]
                .mul(-1)
                .div(scaling_factors["power"]))


def _excavate_co2_total_system(data, scaling_factors):
    return (data.get_formatted_array("cost")
                .to_dataframe()
                .reset_index()
                .groupby(["costs", "locs"])
                .cost
                .sum()
                .loc["co2"]
                .div(scaling_factors["co2"]))


def _excavate_cycles_for_tech(data, tech, scaling_factors):
    charge = data.get_formatted_array("carrier_con").sel(techs=tech)
    storage_cap = data.get_formatted_array("storage_cap").sel(techs=tech)
    return (charge.sum(dim="timesteps") * (-1) / storage_cap)


if __name__ == "__main__":
    excavate_all_results(
        paths_to_scenarios=snakemake.input.scenarios,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
