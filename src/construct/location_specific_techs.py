"""Create renewable techs that are location specific and allocate."""
import jinja2
import pandas as pd
from calliope import AttrDict

TEMPLATE = """
tech_groups:
    open_field_pv_group:
        essentials:
            name: Open field PV
            parent: pv
        constraints:
            resource_area_per_energy_cap: 0.0000125 # [km^2/kW] from (Gagnon:2016, Klauser:2016, Wirth:2017)
            resource: file=capacityfactors-open-field-pv.csv
            resource_unit: energy_per_cap
    roof_mounted_pv_group:
        essentials:
            name: Roof mounted PV
            parent: pv
        constraints:
            resource_area_per_energy_cap: 0.00000625 # [km^2/kW] from (Gagnon:2016, Klauser:2016)
            resource: file=capacityfactors-rooftop-pv.csv
            resource_unit: energy_per_cap
    wind_onshore_group:
        essentials:
            name: Onshore wind
            parent: wind
        constraints:
            resource_area_per_energy_cap: 0.000125 # [km^2/kW] from (European Environment Agency, 2009)
            resource: file=capacityfactors-wind-onshore.csv
            resource_unit: energy_per_cap
        costs:
            monetary:
                energy_cap: 1456 # from IRENA (Figure 5.2), valid for 2016, given in USD2016/kW
                om_annual: 55 # from IRENA (Table 5.1), given in USD2016/kW
    wind_offshore_group:
        essentials:
            name: Offshore wind
            parent: wind
        constraints:
            resource_area_per_energy_cap: 0.000066667 # [km^2/kW] from (European Environment Agency, 2009)
            resource: file=capacityfactors-wind-offshore.csv
            resource_unit: energy_per_cap
        costs:
            monetary:
                energy_cap: 4697 # from IRENA p.101, valid for 2016, given in USD2016/kW
                om_annual: 125 # from IRENA p.109, given in US/kW


techs:
    {% for location in locations %}
        open_field_pv_{{ location }}:
            essentials:
                name: Open field PV in {{ location }}
                parent: open_field_pv_group
        roof_mounted_pv_{{ location }}:
            essentials:
                name: Roof mounted PV in {{ location }}
                parent: roof_mounted_pv_group
        wind_onshore_{{ location }}:
            essentials:
                name: Onshore wind in {{ location }}
                parent: wind_onshore_group
        wind_offshore_{{ location }}:
            essentials:
                name: Offshore wind in {{ location }}
                parent: wind_offshore_group
    {% endfor %}

overrides:
    location_specific_techs:
        {% for location in locations %}
        locations.{{ location }}_roof_mounted_pv.techs.roof_mounted_pv.exists: False
        locations.{{ location }}_pv_or_wind_farm.techs.wind_onshore.exists: False
        locations.{{ location }}_pv_or_wind_farm.techs.open_field_pv.exists: False
        locations.{{ location }}_wind_onshore.techs.wind_onshore.exists: False
        locations.{{ location }}_wind_offshore.techs.wind_offshore.exists: False

        locations.{{ location }}_roof_mounted_pv.techs.roof_mounted_pv_{{ location }}.constraints.resource: file=capacityfactors-rooftop-pv.csv:{{ location }}
        locations.{{ location }}_pv_or_wind_farm.techs.wind_onshore_{{ location }}.constraints.resource: file=capacityfactors-wind-onshore.csv:{{ location }}
        locations.{{ location }}_pv_or_wind_farm.techs.open_field_pv_{{ location }}.constraints.resource: file=capacityfactors-open-field-pv.csv:{{ location }}
        locations.{{ location }}_wind_onshore.techs.wind_onshore_{{ location }}.constraints.resource: file=capacityfactors-wind-onshore.csv:{{ location }}
        locations.{{ location }}_wind_offshore.techs.wind_offshore_{{ location }}.constraints.resource: file=capacityfactors-wind-offshore.csv:{{ location }}
        {% endfor %}
"""


def generate_location_specific_techs(path_to_locations, path_to_output):
    """Create renewable techs that are location specific and allocate to locations."""
    locations = _read_locations(path_to_locations)

    capacities = jinja2.Template(TEMPLATE).render(
        locations=locations
    )
    with open(path_to_output, "w") as result_file:
        result_file.write(capacities)


def _read_locations(path_to_data):
    data = AttrDict.from_yaml(path_to_data)
    return pd.Series([location[:3] for location in data.locations.keys()]).unique()


if __name__ == "__main__":
    generate_location_specific_techs(
        path_to_locations=snakemake.input.locations,
        path_to_output=snakemake.output[0]
    )
