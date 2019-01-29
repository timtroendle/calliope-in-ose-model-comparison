"""Create Calliope link file."""
import jinja2
import pandas as pd
import pycountry

TEMPLATE = """
links:
    {% for row_id, row in ntcs.iterrows() %}
    {{ row.fromLocation }},{{ row.toLocation }}.techs:
        ac_transmission:
            constraints:
                energy_cap_equals: {{ row.capacity_go }}
                one_way: true
    {{ row.toLocation }},{{ row.fromLocation }}.techs:
        ac_transmission:
            constraints:
                energy_cap_equals: {{ row.capacity_return }}
                one_way: true
    {% endfor %}
"""

# TODO allow for transmission capacity extension?


def generate_links(path_to_ntc, path_to_result):
    """Generate a file that represents links in Calliope."""
    ntcs = _read_ntcs(path_to_ntc)
    links = jinja2.Template(TEMPLATE).render(
        ntcs=ntcs
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _read_ntcs(path_to_ntc):
    data = pd.read_excel(path_to_ntc, sheet_name="NTC", skiprows=1)
    data["fromLocation"] = data["Border"].map(
        lambda x: pycountry.countries.lookup(x.split("-")[0][:2]).alpha_3
    )
    data["toLocation"] = data["Border"].map(
        lambda x: pycountry.countries.lookup(x.split("-")[1][:2]).alpha_3
    )
    data["capacity_go"] = data["=>"] * 1e3 # from MW to kW
    data["capacity_return"] = data["<="] * 1e3 # from MW to kW
    # FIXME Northern Ireland is broken
    data.drop(columns=["Border", "=>", "<="], inplace=True)
    data.drop(index=data[data["fromLocation"] == data["toLocation"]].index, inplace=True)
    return (data.groupby(["fromLocation", "toLocation"])
                .agg({"capacity_go": "sum", "capacity_return": "sum"})
                .reset_index())


if __name__ == "__main__":
    generate_links(
        snakemake.input.ntc,
        snakemake.output[0]
    )
