"""Collect all time series results in standardised format."""
from pathlib import Path
from dataclasses import dataclass

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
CARRIER = "carriers"
HOURS = list(range(1, 8761))
NUMBER_HOURS = 8760
SCENARIO_NAME_MAP = {
    "baseline": "Baseline battery costs|Full geographic coverage",
    "low-cost": "50percent battery costs|Full geographic coverage",
    "lowest-cost": "25percent battery costs|Full geographic coverage",
    "baseline-germany": "Baseline battery costs|Germany only",
    "low-cost-germany": "50percent battery costs|Germany only",
    "lowest-cost-germany": "25percent battery costs|Germany only"
}
MODEL_DIM = "Model"
SCENARIO_DIM = "Scenario"
REGION_DIM = "Region"
VARIABLE_DIM = "Variable"
TIME_DIM = "Hour"
TIME_INT_DIM = "Hour"
DATA_NAME = "Value"


@dataclass
class Variable:
    name: str
    value_function: callable = None
    aggregation_function: callable = np.sum


def excavate_all_ts_results(paths_to_scenarios, scaling_factors, path_to_output):
    """Collect all time series results in standardised format"""
    scenarios = {
        Path(path_to_scenario).parent.stem: xr.open_dataset(path_to_scenario)
        for path_to_scenario in paths_to_scenarios
    }
    variables = _set_up_variables()
    scaling_factors = _prepare_scaling_factors(scaling_factors)
    data = _excavate_data(scenarios, variables, scaling_factors)
    new_regions = [_rename_region(region) for region in data[REGION_DIM].to_index()]
    new_scenarios = [SCENARIO_NAME_MAP[scenario] for scenario in data[SCENARIO_DIM].to_index()]
    new_timestamps = range(1, len(data[TIME_DIM]) + 1)
    data = data.assign_coords(Region=new_regions, Scenario=new_scenarios, Hour=new_timestamps)
    data = data.expand_dims(Model=[MODEL_NAME], axis=0)
    data = data.to_dataframe(name=DATA_NAME).dropna().reset_index() # for CSV writing use pandas
    data[[MODEL_DIM, SCENARIO_DIM, REGION_DIM, VARIABLE_DIM, TIME_DIM, DATA_NAME]].to_csv(
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
    scenario_coords = list(scenarios.keys())
    region_coords = ALL_COUNTRIES + [TOTAL_COVERED_REGION]
    variable_coords = [variable.name for variable in variables]
    hour_coords = list(scenarios.values())[0].timesteps
    data = xr.DataArray(
        np.zeros(shape=(len(scenario_coords), len(region_coords), len(variable_coords), len(hour_coords))) * np.nan,
        coords=[scenario_coords, region_coords, variable_coords, hour_coords],
        dims=[SCENARIO_DIM, REGION_DIM, VARIABLE_DIM, TIME_DIM]
    )
    for variable in variables:
        if variable.value_function:
            data = _excavate_variable(data, scenarios, variable, scaling_factors)
    return data


def _excavate_variable(data, scenarios, variable, scaling_factors):
    values = {
        scenario_name: (variable.value_function(scenario_data, scaling_factors)
                                .rename(locs=REGION_DIM, timesteps=TIME_DIM)
                                .expand_dims(Scenario=[scenario_name], axis=0)
                                .expand_dims(Variable=[variable.name], axis=2))
        for scenario_name, scenario_data in scenarios.items()
    }
    for scenario, scenario_values in values.items():
        assert len(scenario_values.dims) == 4, scenario_values
        if "germany" not in scenario:
            data = data.combine_first(scenario_values)
        total_scenario_values = (scenario_values.reduce(variable.aggregation_function, dim="Region")
                                                .expand_dims(Region=[TOTAL_COVERED_REGION]))
        data = data.combine_first(total_scenario_values)
    return data


def _prepare_scaling_factors(scaling_factors):
    scaling_factors["energy"] = scaling_factors["power"] / 1e-6 # from MWh to TWh
    scaling_factors["stored_energy"] = scaling_factors["power"] / 1e-3 # from MWh to GWh
    scaling_factors["specific_cost"] = scaling_factors["monetary"] / scaling_factors["power"]
    scaling_factors["power"] = scaling_factors["power"] / 1e-3 # from MW to GW
    return scaling_factors


def _set_up_variables():
    return [
        Variable("Hourly usage|Electricity|Storage|Liion|Level", lambda data, sf: _level(data, "battery", sf), np.mean),
        Variable("Hourly usage|Electricity|Storage|Liion|Storage loading", lambda data, sf: _charging(data, "battery", sf)),
        Variable("Hourly usage|Electricity|Storage|Liion|Storage discharging", lambda data, sf: _discharging(data, "battery", sf)),
        Variable("Hourly usage|Electricity|Storage|Pumped hydro|Level", lambda data, sf: _level(data, "pumped_hydro", sf), np.mean),
        Variable("Hourly usage|Electricity|Storage|Pumped hydro|Storage loading", lambda data, sf: _charging(data, "pumped_hydro", sf)),
        Variable("Hourly usage|Electricity|Storage|Pumped hydro|Storage discharging", lambda data, sf: _discharging(data, "pumped_hydro", sf)),
        Variable("Hourly usage|Electricity|Storage|Other|Level"),
        Variable("Hourly usage|Electricity|Storage|Other|Storage loading"),
        Variable("Hourly usage|Electricity|Storage|Other|Storage discharging"),
    ]


def _level(data, tech, scaling_factors):
    stored_energy = _split_loc_techs(data["storage"]).sel(techs=tech)
    capacity = _split_loc_techs(data["storage_cap"]).sel(techs=tech)
    return stored_energy / capacity


def _charging(data, tech, scaling_factors):
    return _split_loc_techs(data["carrier_prod"]).squeeze(CARRIER).sel(techs=tech) / scaling_factors["power"]


def _discharging(data, tech, scaling_factors):
    return _split_loc_techs(data["carrier_prod"]).squeeze(CARRIER).sel(techs=tech) / scaling_factors["power"]


def _split_loc_techs(data_var, as_='DataArray'):
    # copy from Calliope source code, copyright Calliope authors
    # I am using the copy here because Calliope does not support xarray 0.12.
    # which I need for the data handling. Thus, I cannot import calliope here.
    """
    Get a DataArray with locations technologies, and possibly carriers
    split into separate coordinates.

    Parameters
    ----------
    data_var : xarray DataArray
        Variable from Calliope model_data, to split loc_techs dimension
    as_ : string
        'DataArray' to return xarray DataArray or 'Series' to return pandas
        Series with dimensions as a MultiIndex

    Returns
    -------
    updated_data_var : xarray DataArray of pandas Series
    """

    # Separately find the loc_techs(_carriers) dimension and all other dimensions
    loc_tech_dim = [i for i in data_var.dims if 'loc_tech' in i]
    if not loc_tech_dim:
        loc_tech_dim = [i for i in data_var.dims if 'loc_carrier' in i]
    non_loc_tech_dims = list(set(data_var.dims).difference(loc_tech_dim))

    if not loc_tech_dim:
        if as_ == 'Series':
            return data_var.to_series()
        elif as_ == 'DataArray':
            return data_var
        else:
            raise ValueError('`as_` must be `DataArray` or `Series`, '
                             'but `{}` given'.format(as_))

    elif len(loc_tech_dim) > 1:
        raise ValueError("Cannot split loc_techs or loc_techs_carrier dimension "
                         "for DataArray {}".format(data_var.name))

    loc_tech_dim = loc_tech_dim[0]
    # xr.Datarray -> pd.Series allows for string operations
    data_var_df = data_var.to_series().unstack(non_loc_tech_dims)
    index_list = data_var_df.index.str.split('::').tolist()

    # carrier_prod, carrier_con, and carrier_export will return an index_list
    # of size 3, all others will be an index list of size 2
    possible_names = ['loc', 'tech', 'carrier']
    names = [i + 's' for i in possible_names if i in loc_tech_dim]

    data_var_df.index = pd.MultiIndex.from_tuples(index_list, names=names)

    # If there were no other dimensions other than loc_techs(_carriers) then
    # nothing was unstacked on creating data_var_df, so nothing is stacked now
    if isinstance(data_var_df, pd.Series):
        data_var_series = data_var_df
    else:
        data_var_series = data_var_df.stack(non_loc_tech_dims)

    if as_ == "Series":
        return data_var_series

    elif as_ == "DataArray":
        updated_data_var = xr.DataArray.from_series(data_var_series)
        updated_data_var.attrs = data_var.attrs
        updated_data_var.name = data_var.name

        return updated_data_var

    else:
        raise ValueError('`as_` must be `DataArray` or `Series`, '
                         'but `{}` given'.format(as_))


if __name__ == "__main__":
    excavate_all_ts_results(
        paths_to_scenarios=snakemake.input.scenarios,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
