"""Create Calliope file defining CO2 caps."""
import jinja2
import pandas as pd
import pycountry


TEMPLATE = """
group_constraints:
    {% for country, cap in caps.iteritems() %}
    co2_cap_{{ country }}:
        locs: ["{{ country }}"]
        cost_max:
            co2: {{ cap }}
    {% endfor %}
"""


def generate_co2_caps(path_to_caps, path_to_yaml, path_to_csv):
    """Create Calliope file defining CO2 caps."""
    raw = _read_caps(path_to_caps)
    raw.to_csv(path_to_csv, index=True, header=True)

    with open(path_to_yaml, "w") as result_file:
        result_file.write(jinja2.Template(TEMPLATE).render(caps=raw["co2_cap"]))


def _read_caps(path_to_data):
    data = pd.read_excel(
        path_to_data,
        sheet_name="bounds on RES and CO2",
        index_col=0,
        usecols="A,C"
    ).iloc[:-3, :]
    data.drop(index=["MT"], inplace=True)
    data.rename(index={"UK": "GB"}, inplace=True)
    data.rename(
        index=lambda iso2: pycountry.countries.lookup(iso2).alpha_3,
        inplace=True
    )
    data.rename(
        columns={"2030 Power generation/District heating (Mt of CO2 eq.)": "co2_cap"},
        inplace=True
    )
    return data * 1e3 # from Mt to kt


if __name__ == "__main__":
    generate_co2_caps(
        path_to_caps=snakemake.input.caps,
        path_to_yaml=snakemake.output.yaml,
        path_to_csv=snakemake.output.csv
    )
