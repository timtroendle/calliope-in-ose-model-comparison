"""Create Calliope file defining bounds of generation capacities."""
import jinja2
import pandas as pd
import pycountry

BIOMASS_CAPACITY_FACTOR = 7500 / 8760 # ASSUME by DIW

TECH_MAP = {
    "Biofuels": "biofuels", # ASSUME inexistent
    "Gas": "ccgt", # ASSUME all gas to be CCGT
    "Hard coal": "coal",
    "Hydro-pump": "pumped_hydro",
    "Hydro-run": "hydro_run_of_river",
    "Hydro-turbine": "hydro_turbine", # ASSUME ignore
    "Lignite": "lignite",
    "Nuclear": "nuclear",
    "Oil": "oil", # ASSUME doesnt exist in electricity (was 1.9% in 2015)
    "Othernon-RES": "ccgt", # ASSUME ccgt
    "Other RES": "biomass", # ASSUME biomass
    "Solar-thermal": "solar_thermal", # ignored for electricity
    "Solar-PV": "roof_mounted_pv", # ASSUME all solar pv is roof mounted
    "Wind-on-shore": "wind_onshore_monopoly",
    "Wind-off-shore": "wind_offshore"
}

TEMPLATE = """
locations:
    {% for country, techs in capacities.iterrows() %}
    {{ country }}:
        techs:
            wind_onshore_monopoly:
                constraints:
                    energy_cap_min: {{ techs.wind_onshore_monopoly }} # [{{ unit_scaling_factor }} MW]
            wind_offshore:
                constraints:
                    energy_cap_min: {{ techs.wind_offshore }} # [{{ unit_scaling_factor }} MW]
            roof_mounted_pv:
                constraints:
                    energy_cap_min: {{ techs.roof_mounted_pv }} # [{{ unit_scaling_factor }} MW]
            hydro_run_of_river:
                constraints:
                    energy_cap_equals: {{ techs.hydro_run_of_river }} # [{{ unit_scaling_factor }} MW]
            biomass:
                constraints:
                    energy_cap_equals: {{ techs.biomass }} # [{{ unit_scaling_factor }} MW]
            pumped_hydro:
                constraints:
                    energy_cap_equals: {{ techs.pumped_hydro }} # [{{ unit_scaling_factor }} MW]
            coal:
                constraints:
                    energy_cap_max: {{ techs.coal }} # [{{ unit_scaling_factor }} MW]
            lignite:
                constraints:
                    energy_cap_max: {{ techs.lignite }} # [{{ unit_scaling_factor }} MW]
            ccgt:
                constraints:
                    energy_cap_max: {{ techs.ccgt }} # [{{ unit_scaling_factor }} MW]
            nuclear:
                constraints:
                    energy_cap_max: {{ techs.nuclear }} # [{{ unit_scaling_factor }} MW]
    {% endfor %}
"""


def generate_generation_capacities(path_to_data, scaling_factor, path_to_yaml, path_to_csv):
    """Create Calliope file defining bounds of generation capacities."""
    raw = _read_generation_capacities(path_to_data)
    raw.loc[:, "biomass"] = raw.loc[:, "biomass"] * BIOMASS_CAPACITY_FACTOR
    raw.to_csv(path_to_csv, index=True, header=True)

    capacities = jinja2.Template(TEMPLATE).render(
        capacities=raw.mul(scaling_factor),
        unit_scaling_factor=1 / scaling_factor
    )
    with open(path_to_yaml, "w") as result_file:
        result_file.write(capacities)


def _read_generation_capacities(path_to_data):
    data = pd.read_excel(path_to_data, sheet_name="2030 ST", skiprows=2, index_col=0, usecols="A:P")
    data.drop(
        index=[
            "DEkf", # TYNDP 2018 data, but we are using NEP
            "DE", # TYNDP 2018 data, but we are using NEP
            "DK_total", # data for sub nodes available, danger of counting data twice
            "IT_total", # data for sub nodes available, danger of counting data twice
            "LU_total", # data for sub nodes available, danger of counting data twice
            "NO_total", # data for sub nodes available, danger of counting data twice
            "SE_total", # data for sub nodes available, danger of counting data twice
            "TR", # Turkey out of scope
            "IL00", # Israel out of scope
            "TN00", # Tunesia out of scope
            " IS00", # Israel out of scope
            "MT", # Malta out of scope
        ],
        inplace=True
    )
    data.rename(
        index=lambda name: name if name != "NI" else "GB_NI",
        inplace=True
    )
    data.rename(
        index=lambda name: pycountry.countries.lookup(name.strip()[:2]).alpha_3,
        inplace=True
    )
    data.rename(
        columns=lambda name: TECH_MAP[name.strip().replace('\n', '').replace('\r', '')],
        inplace=True
    )
    return data.groupby(data.index).sum().groupby(data.columns, axis="columns").sum()


if __name__ == "__main__":
    generate_generation_capacities(
        path_to_data=snakemake.input.capacity,
        path_to_yaml=snakemake.output.yaml,
        path_to_csv=snakemake.output.csv,
        scaling_factor=snakemake.params.scaling_factor
    )
