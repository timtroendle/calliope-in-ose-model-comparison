"""Create Calliope file defining minimal renewable shares."""
import jinja2
import pandas as pd
import pycountry


TEMPLATE = """
model:
    group_share:
        {% for country, share in shares.iteritems() %}
        wind_onshore_{{ country }},wind_offshore_{{ country }},open_field_pv_{{ country }},roof_mounted_pv_{{ country }}:
            carrier_prod_min:
                electricity: {{ share }}
        {% endfor %}
"""


def generate_renewable_shares(path_to_shares, path_to_demand, from_date, to_date, path_to_yaml, path_to_csv):
    """Create Calliope file defining minimal renewable shares."""
    raw = _read_shares(path_to_shares)
    raw.to_csv(path_to_csv, index=True, header=True)

    share_of_global_demand = _read_share_of_global_demand(path_to_demand, from_date, to_date).reindex(raw.index)
    global_share = raw["renewable_share"] * share_of_global_demand

    with open(path_to_yaml, "w") as result_file:
        result_file.write(jinja2.Template(TEMPLATE).render(shares=global_share))


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


def _read_share_of_global_demand(path_to_demand, from_date, to_date):
    return pd.read_csv(path_to_demand, index_col=0, parse_dates=True).loc[from_date:to_date, :].sum(
        axis="index"
    ).transform(
        lambda x: x / x.sum()
    )


if __name__ == "__main__":
    generate_renewable_shares(
        path_to_shares=snakemake.input.shares,
        path_to_demand=snakemake.input.demand,
        from_date=snakemake.config["from_date"],
        to_date=snakemake.config["to_date"],
        path_to_yaml=snakemake.output.yaml,
        path_to_csv=snakemake.output.csv
    )
