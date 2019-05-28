"""Create Calliope file defining charge rate of pumped hydro."""
import jinja2
import numpy as np
import pandas as pd

CHARGE_RATES = { # from [@Geth:2015], given in [1/h]
    "AUT": 3.246 / 132.41,
    "BEL": 1.196 / 5.71,
    "BGR": 0.93 / 11.13,
    "HRV": 0.246 / 2.34,
    "CYP": np.nan,
    "CZE": 1.145 / 5.72,
    "DNK": np.nan,
    "EST": np.nan,
    "FIN": np.nan,
    "FRA": 4.317 / 83.37,
    "DEU": 1 / 7, # ASSUME DIW assumption
    "GRC": 0.735 / 4.97,
    "HUN": np.nan,
    "IRL": 0.292 / 1.8,
    "ITA": 7.64 / 68.27,
    "LVA": np.nan,
    "LTU": 0.88 / 10.8,
    "LUX": 1.050 / 4.92,
    "MLT": np.nan,
    "NLD": np.nan,
    "POL": 1.647 / 7.96,
    "PRT": 1.279 / 40.77,
    "ROU": 0.2 / 10.2,
    "SVK": 0.79 / 3.63,
    "SVN": 0.180 / 0.5,
    "ESP": 5.859 / 70,
    "SWE": 0.091 / 72.12,
    "GBR": 2.65 / 26.7,
    "NOR": 0.892 / 399.39,
    "CHE": 1.512 / 311.48,
    "ALB": np.nan, # has no pumped hydro
    "BIH": 0.735 / 4.97, # ASSUME like Greece
    "MKD": np.nan, # has no pumped hydro
    "MNE": np.nan, # has no pumped hydro
    "SRB": 0.735 / 4.97, # ASSUME like Greece
}

TEMPLATE = """
locations:
    {% for country, charge_rate in charge_rates.items() %}
    {{ country }}:
        techs:
            pumped_hydro:
                constraints:
                    energy_cap_per_storage_cap_equals: {{ charge_rate }} # [1/h]
    {% endfor %}
"""


def generate_pumped_hydro(path_to_yaml):
    """Create Calliope file defining charge rate of pumped hydro."""
    rates = {country: rate for country, rate in CHARGE_RATES.items() if not pd.isnull(rate)}
    charge_rates = jinja2.Template(TEMPLATE).render(
        charge_rates=rates
    )
    with open(path_to_yaml, "w") as result_file:
        result_file.write(charge_rates)


if __name__ == "__main__":
    generate_pumped_hydro(
        path_to_yaml=snakemake.output.yaml
    )
