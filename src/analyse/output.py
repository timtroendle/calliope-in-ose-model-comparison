"""Collect all results in standardised format."""
import functools
from pathlib import Path
from dataclasses import dataclass

import calliope
import numpy as np
import pandas as pd
import xarray as xr
import pycountry


MODEL_NAME = "euro-calliope"
TOTAL_COVERED_REGION = "Total covered region"
GERMANY = "DEU"
ALL_COUNTRIES = ['ALB', 'AUT', 'BEL', 'BGR', 'BIH', 'CHE', 'CYP', 'CZE',
                 'DEU', 'DNK', 'ESP', 'EST', 'FIN', 'FRA', 'GBR', 'GRC', 'HRV', 'HUN',
                 'IRL', 'ITA', 'LTU', 'LUX', 'LVA', 'MKD', 'MNE', 'NLD', 'NOR', 'POL',
                 'PRT', 'ROU', 'SRB', 'SVK', 'SVN', 'SWE']
LOCATIONS = [TOTAL_COVERED_REGION, GERMANY]
CARRIER = "electricity"
NUMBER_HOURS = 8760
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly",
            "wind_onshore_competing", "wind_offshore", "hydro_run_of_river"]
SCENARIO_NAME_MAP = {
    "baseline": "Baseline battery costs|Full geographic coverage",
    "low-cost": "50percent battery costs|Full geographic coverage",
    "lowest-cost": "25percent battery costs|Full geographic coverage",
    "baseline-germany": "Baseline battery costs|Germany only",
    "low-cost-germany": "50percent battery costs|Germany only",
    "lowest-cost-germany": "25percent battery costs|Germany only"
}


@dataclass
class Variable:
    name: str
    value_function: callable = None
    aggregation_function: callable = np.sum


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
    data["region"] = data["region"].map(_rename_region)
    data[["model", "scenario", "region", "variable", "value"]].to_csv(
        path_to_output,
        header=True,
        index=False
    )


def _rename_region(region):
    if region == TOTAL_COVERED_REGION:
        return region
    else:
        return pycountry.countries.lookup(region).alpha_2


def _excavate_data(scenarios, variables, scaling_factors):
    european_scenarios = [scenario for scenario in scenarios.keys() if "germany" not in scenario]
    germany_scenarios = [scenario for scenario in scenarios.keys() if "germany" in scenario]
    index_europe = pd.MultiIndex.from_product(
        [european_scenarios, ALL_COUNTRIES + [TOTAL_COVERED_REGION], [variable.name for variable in variables]],
        names=["scenario", "region", "variable"]
    )
    index_germany = pd.MultiIndex.from_product(
        [germany_scenarios, [TOTAL_COVERED_REGION], [variable.name for variable in variables]],
        names=["scenario", "region", "variable"]
    )
    data = pd.concat([
        pd.Series(index=index_europe, name="value"),
        pd.Series(index=index_germany, name="value")
    ])
    for variable in variables:
        if variable.value_function:
            _excavate_variable(data, scenarios, variable, scaling_factors)
    return data


def _excavate_variable(data, scenarios, variable, scaling_factors):
    values = {
        scenario_name: variable.value_function(scenario_data, scaling_factors)
        for scenario_name, scenario_data in scenarios.items()
    }
    for scenario, scenario_values in values.items():
        if "germany" not in scenario:
            for country in ALL_COUNTRIES:
                data.loc[(scenario, country, variable.name)] = scenario_values.loc[country]
        try:
            total_value = scenario_values.agg(variable.aggregation_function)
        except AttributeError: # scenario_values are xarray DataArray
            assert isinstance(scenario_values, xr.DataArray)
            total_value = scenario_values.to_pandas().agg(variable.aggregation_function)
        data.loc[(scenario, TOTAL_COVERED_REGION, variable.name)] = total_value


def _prepare_scaling_factors(scaling_factors):
    scaling_factors["energy"] = scaling_factors["power"] / 1e-6 # from MWh to TWh
    scaling_factors["stored_energy"] = scaling_factors["power"] / 1e-3 # from MWh to GWh
    scaling_factors["specific_cost"] = scaling_factors["monetary"] / scaling_factors["power"]
    scaling_factors["power"] = scaling_factors["power"] / 1e-3 # from MW to GW
    return scaling_factors


def _set_up_variables():
    return [
        Variable("Cost|Total system", _excavate_cost_total_system),
        Variable("Cost|Capacity", _excavate_cost_capacity),
        Variable("Cost|Average total system", _excavate_levelised_total_cost),
        Variable("Cost|Variable", _excavate_cost_variable),
        Variable("Cost|Network"),
        Variable("Capacity|Electricity|Nuclear", lambda data, sf: _excavate_capacity_for_tech(data, ["nuclear"], sf)),
        Variable("Capacity|Electricity|Lignite", lambda data, sf: _excavate_capacity_for_tech(data, ["lignite"], sf)),
        Variable("Capacity|Electricity|Hard coal", lambda data, sf: _excavate_capacity_for_tech(data, ["coal"], sf)),
        Variable("Capacity|Electricity|Gas CCGT", lambda data, sf: _excavate_capacity_for_tech(data, ["ccgt"], sf)),
        Variable("Capacity|Electricity|Gas OCGT"),
        Variable("Capacity|Electricity|Oil", lambda data, sf: _excavate_capacity_for_tech(data, ["oils"], sf)),
        Variable("Capacity|Electricity|Other nonrenewable"),
        Variable("Capacity|Electricity|Wind onshore", lambda data, sf: _excavate_capacity_for_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf)),
        Variable("Capacity|Electricity|Wind offshore", lambda data, sf: _excavate_capacity_for_tech(data, ["wind_offshore"], sf)),
        Variable("Capacity|Electricity|Solar PV", lambda data, sf: _excavate_capacity_for_tech(data, ["open_field_pv", "roof_mounted_pv"], sf)),
        Variable("Capacity|Electricity|Solar CSP"),
        Variable("Capacity|Electricity|Hydro ROR", lambda data, sf: _excavate_capacity_for_tech(data, ["hydro_run_of_river"], sf)),
        Variable("Capacity|Electricity|Hydro Reservoir", lambda data, sf: _excavate_capacity_for_tech(data, ["hydro_reservoir"], sf)),
        Variable("Capacity|Electricity|Bioenergy", lambda data, sf: _excavate_capacity_for_tech(data, ["biomass"], sf)),
        Variable("Capacity|Electricity|Other renewable"),
        Variable("Capacity|Electricity|Storage|Liion|Power", lambda data, sf: _excavate_capacity_for_tech(data, ["battery"], sf)),
        Variable("Capacity|Electricity|Storage|Liion|Energy", lambda data, sf: _excavate_storage_capacity_for_tech(data, ["battery"], sf)),
        Variable("Capacity|Electricity|Storage|Pumped hydro|Power", lambda data, sf: _excavate_capacity_for_tech(data, ["pumped_hydro"], sf)),
        Variable("Capacity|Electricity|Storage|Pumped hydro|Energy", lambda data, sf: _excavate_storage_capacity_for_tech(data, ["pumped_hydro"], sf)),
        Variable("Capacity|Electricity|Storage|Other|Power"),
        Variable("Capacity|Electricity|Storage|Other|Energy"),
        Variable("Capacity|Electricity|Gross import", _excavate_transmission_capacity),
        Variable("Capacity|Electricity|Gross export", _excavate_transmission_capacity),
        Variable("Capacity|Sector coupling|Electric vehicles"),
        Variable("Capacity|Sector coupling|Powertoheat"),
        Variable("Capacity|Sector coupling|Electrolysis"),
        Variable("Energy|Electricity|Nuclear", lambda data, sf: _excavate_generation_for_tech(data, ["nuclear"], sf)),
        Variable("Energy|Electricity|Lignite", lambda data, sf: _excavate_generation_for_tech(data, ["lignite"], sf)),
        Variable("Energy|Electricity|Hard coal", lambda data, sf: _excavate_generation_for_tech(data, ["coal"], sf)),
        Variable("Energy|Electricity|Gas CCGT", lambda data, sf: _excavate_generation_for_tech(data, ["ccgt"], sf)),
        Variable("Energy|Electricity|Gas OCGT"),
        Variable("Energy|Electricity|Oil", lambda data, sf: _excavate_generation_for_tech(data, ["oils"], sf)),
        Variable("Energy|Electricity|Other nonrenewable"),
        Variable("Energy|Electricity|Wind onshore", lambda data, sf: _excavate_generation_for_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf)),
        Variable("Energy|Electricity|Wind offshore", lambda data, sf: _excavate_generation_for_tech(data, ["wind_offshore"], sf)),
        Variable("Energy|Electricity|Solar PV", lambda data, sf: _excavate_generation_for_tech(data, ["open_field_pv", "roof_mounted_pv"], sf)),
        Variable("Energy|Electricity|Solar CSP"),
        Variable("Energy|Electricity|Hydro ROR", lambda data, sf: _excavate_generation_for_tech(data, ["hydro_run_of_river"], sf)),
        Variable("Energy|Electricity|Hydro Reservoir", lambda data, sf: _excavate_generation_for_tech(data, ["hydro_reservoir"], sf)),
        Variable("Energy|Electricity|Bioenergy", lambda data, sf: _excavate_generation_for_tech(data, ["biomass"], sf)),
        Variable("Energy|Electricity|Other renewable"),
        Variable("Energy|Electricity|Renewable curtailment|Absolute", _excavate_absolute_curtailment),
        Variable("Energy|Electricity|Renewable curtailment|Relative", _excavate_relative_curtailment, np.mean),
        Variable("Energy|Electricity|Storage|Liion", lambda data, sf: _excavate_generation_for_tech(data, ["battery"], sf)),
        Variable("Energy|Electricity|Storage|Pumped hydro", lambda data, sf: _excavate_generation_for_tech(data, ["pumped_hydro"], sf)),
        Variable("Energy|Electricity|Storage|Other"),
        Variable("Energy|Electricity|Gross import", _excavate_transmission_generation),
        Variable("Energy|Electricity|Gross export", _excavate_transmission_consumption),
        Variable("Full load hours|Electricity|Nuclear", lambda data, sf: _full_load_hours_tech(data, ["nuclear"], sf), np.mean),
        Variable("Full load hours|Electricity|Lignite", lambda data, sf: _full_load_hours_tech(data, ["lignite"], sf), np.mean),
        Variable("Full load hours|Electricity|Hard coal", lambda data, sf: _full_load_hours_tech(data, ["coal"], sf), np.mean),
        Variable("Full load hours|Electricity|Gas CCGT", lambda data, sf: _full_load_hours_tech(data, ["ccgt"], sf), np.mean),
        Variable("Full load hours|Electricity|Gas OCGT"),
        Variable("Full load hours|Electricity|Oil", lambda data, sf: _full_load_hours_tech(data, ["oils"], sf), np.mean),
        Variable("Full load hours|Electricity|Other nonrenewable"),
        Variable("Full load hours|Electricity|Wind onshore", lambda data, sf: _full_load_hours_tech(data, ["wind_onshore_monopoly", "wind_onshore_competing"], sf), np.mean),
        Variable("Full load hours|Electricity|Wind offshore", lambda data, sf: _full_load_hours_tech(data, ["wind_offshore"], sf), np.mean),
        Variable("Full load hours|Electricity|Solar PV", lambda data, sf: _full_load_hours_tech(data, ["open_field_pv", "roof_mounted_pv"], sf), np.mean),
        Variable("Full load hours|Electricity|Solar CSP"),
        Variable("Full load hours|Electricity|Hydro ROR", lambda data, sf: _full_load_hours_tech(data, ["hydro_run_of_river"], sf), np.mean),
        Variable("Full load hours|Electricity|Hydro Reservoir", lambda data, sf: _full_load_hours_tech(data, ["hydro_reservoir"], sf), np.mean),
        Variable("Full load hours|Electricity|Bioenergy", lambda data, sf: _full_load_hours_tech(data, ["biomass"], sf), np.mean),
        Variable("Full load hours|Electricity|Other renewable"),
        Variable("Full load hours|Electricity|Storage|Liion", lambda data, sf: _full_load_hours_tech(data, ["battery"], sf), np.mean),
        Variable("Full load hours|Electricity|Storage|Pumped hydro", lambda data, sf: _full_load_hours_tech(data, ["pumped_hydro"], sf), np.mean),
        Variable("Full load hours|Electricity|Storage|Other"),
        Variable("Cycles|Electricity|Storage|Liion", lambda data, sf: _excavate_cycles_for_tech(data, ["battery"], sf)),
        Variable("Cycles|Electricity|Storage|Pumped hydro", lambda data, sf: _excavate_cycles_for_tech(data, ["pumped_hydro"], sf)),
        Variable("Cycles|Electricity|Storage|Other"),
        Variable("Carbon emissions", _excavate_co2_total_system),
        Variable("Price|Electricity|Weighted average"),
        Variable("Price|Electricity|Storage|Liion|Weighted average charging"),
        Variable("Price|Electricity|Storage|Liion|Weighted average discharging"),
        Variable("Market value|Electricity|Wind onshore"),
        Variable("Market value|Electricity|Wind offshore"),
        Variable("Market value|Electricity|Solar PV"),
        Variable("Market value|Electricity|Solar CSP"),
        Variable("Startups|Electricity|Total number"),
        Variable("Startups|Electricity|Total cost"),
        Variable("Energy|Electricity|Peak demand", _excavate_peak_demand),
        Variable("Energy|Electricity|Total demand", _excavate_demand),
    ]


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
    series = pd.Series(index=ALL_COUNTRIES, dtype=np.float32)
    series.loc[GERMANY] = lcoe
    return series


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
        return pd.Series([0], index=[GERMANY])
    else:
        return capacity


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_capacity(data, scaling_factor):
    return (data.get_formatted_array("energy_cap")
                .to_dataframe()
                .reset_index()
                .set_index(["techs", "locs"])
                .div(scaling_factor))


@functools.lru_cache(maxsize=10, typed=False)
def _excavate_capacity_xarray(data, scaling_factor):
    return data.get_formatted_array("energy_cap") / scaling_factor


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
        return pd.Series([0], index=[GERMANY])
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
def _excavate_generation_xarray(data, scaling_factor):
    return (data.get_formatted_array("carrier_prod")
                .sel(carriers=CARRIER)
                .sum(dim="timesteps")) / scaling_factor


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
        return pd.Series([0], index=[GERMANY])
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
    charge = data.get_formatted_array("carrier_con").sel(techs=tech, carriers=CARRIER).sum(dim="techs")
    storage_cap = data.get_formatted_array("storage_cap").sel(techs=tech).sum(dim="techs")
    return (charge.sum(dim="timesteps") * (-1) / storage_cap)


def _full_load_hours_tech(data, tech, scaling_factors):
    generation = _excavate_generation_xarray(data, scaling_factors["power"]).sel(techs=tech).sum(dim="techs")
    capacity = _excavate_capacity_xarray(data, scaling_factors["power"]).sel(techs=tech).sum(dim="techs")
    return generation / (capacity * NUMBER_HOURS)


if __name__ == "__main__":
    excavate_all_results(
        paths_to_scenarios=snakemake.input.scenarios,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
