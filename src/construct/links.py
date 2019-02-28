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
                energy_cap_equals: {{ row.capacity_go }} # MW
                one_way: true
    {{ row.toLocation }},{{ row.fromLocation }}.techs:
        ac_transmission:
            constraints:
                energy_cap_equals: {{ row.capacity_return }} # MW
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
    data.replace(
        to_replace={
            "NIC": "GBR", # NI is not Northern Ireland, but Nicaragua
            "TUR": None, # out of scope
            "MLT": None, # out of scope
            "TUN": None # out of scope
        },
        inplace=True
    )
    data.dropna(axis="index", how="any", subset=["fromLocation", "toLocation"], inplace=True)
    data.drop(index=data[data["fromLocation"] == data["toLocation"]].index, inplace=True)
    data["capacity_go"] = data["=>"]
    data["capacity_return"] = data["<="]
    data.drop(columns=["Border", "=>", "<="], inplace=True)

    # invert power line from IRL to GBR otherwise it will appear twice in YAML
    idx = data[(data["fromLocation"] == "IRL") & (data["toLocation"] == "GBR")].index
    assert len(idx) == 1
    idx = idx[0]
    from_irl_to_gbr = data.loc[idx, "capacity_go"]
    from_gbr_to_irl = data.loc[idx, "capacity_return"]
    data.loc[idx, "fromLocation"] = "GBR"
    data.loc[idx, "toLocation"] = "IRL"
    data.loc[idx, "capacity_go"] = from_gbr_to_irl
    data.loc[idx, "capacity_return"] = from_irl_to_gbr

    return (data.groupby(["fromLocation", "toLocation"])
                .agg({"capacity_go": "sum", "capacity_return": "sum"})
                .reset_index())


if __name__ == "__main__":
    generate_links(
        snakemake.input.ntc,
        snakemake.output[0]
    )
