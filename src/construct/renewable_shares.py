"""Create Calliope file defining minimal renewable shares."""
import jinja2
import pandas as pd
import pycountry


TEMPLATE = """
group_constraints:
    {% for country, share in shares.iteritems() %}
    renewable_share_{{ country }}:
        locs: ["{{ country }}"]
        techs: ["wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore",
                "open_field_pv", "roof_mounted_pv", "hydro_run_of_river", "hydro_reservoir",
                "biomass"]
        demand_share_min:
            electricity: {{ share }}
    {% endfor %}
"""


def generate_renewable_shares(path_to_shares, path_to_yaml, path_to_csv):
    """Create Calliope file defining minimal renewable shares."""
    raw = _read_shares(path_to_shares)
    raw.to_csv(path_to_csv, index=True, header=True)

    with open(path_to_yaml, "w") as result_file:
        result_file.write(jinja2.Template(TEMPLATE).render(shares=raw["renewable_share"]))


def _read_shares(path_to_data):
    data = pd.read_excel(
        path_to_data,
        sheet_name="bounds on RES and CO2",
        index_col=0,
        usecols="A:B"
    ).iloc[:-3, :]
    data.drop(index=["MT"], inplace=True)
    data.rename(index={"UK": "GB"}, inplace=True)
    data.rename(
        index=lambda iso2: pycountry.countries.lookup(iso2).alpha_3,
        inplace=True
    )
    data.rename(
        columns={"2030 RES-E (%)-RES in Gross Final Electricity Consumption": "renewable_share"},
        inplace=True
    )
    return data / 100 # make unit less


if __name__ == "__main__":
    generate_renewable_shares(
        path_to_shares=snakemake.input.shares,
        path_to_yaml=snakemake.output.yaml,
        path_to_csv=snakemake.output.csv
    )
