"""Create Calliope file defining bounds of generation capacities."""
import jinja2
import pandas as pd
import pycountry

TECH_MAP = {
    "Biofuels": "biofuels",
    "Gas": "gas",
    "Hard coal": "hard_coal",
    "Hydro-pump": "hydro_pump",
    "Hydro-run": "hydro_run",
    "Hydro-turbine": "hydro_turbine",
    "Lignite": "lignite",
    "Nuclear": "nuclear",
    "Oil": "oil",
    "Othernon-RES": "other_non_res",
    "Other RES": "other_res",
    "Solar-thermal": "solar_thermal",
    "Solar-PV": "solar_pv",
    "Wind-on-shore": "wind_onshore",
    "Wind-off-shore": "wind_offshore"
}

TEMPLATE = """
locations:
    {% for country, techs in capacities.iterrows() %}
    {{ country }}_wind_onshore:
        techs:
            wind_onshore:
                constraints:
                    energy_cap_min: {{ techs.wind_onshore }}
    {% endfor %}
"""


def generate_generation_capacities(path_to_data, path_to_yaml, path_to_csv):
    """Create Calliope file defining bounds of generation capacities."""
    raw = _read_generation_capacities(path_to_data)
    raw.to_csv(path_to_csv, index=True, header=True)

    capacities = jinja2.Template(TEMPLATE).render(
        capacities=raw
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
    # TODO drop and properly name techs
    data.rename(
        columns=lambda name: TECH_MAP[name.strip().replace('\n', '').replace('\r', '')],
        inplace=True
    )
    return data.groupby(data.index).sum() * 1000 # from MW to kW


if __name__ == "__main__":
    generate_generation_capacities(
        path_to_data=snakemake.input.capacity,
        path_to_yaml=snakemake.output.yaml,
        path_to_csv=snakemake.output.csv
    )
