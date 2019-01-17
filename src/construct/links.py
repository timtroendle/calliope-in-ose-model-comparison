"""Create Calliope link file."""
from pathlib import Path

import pandas as pd
import pycountry
from calliope import AttrDict


def generate_links(path_to_ntc, path_to_result):
    """Generate a file that represents links in Calliope."""
    ntcs = _read_ntcs(path_to_ntc)
    locations = _generate_links(ntcs)
    locations.to_yaml(Path(path_to_result))


def _read_ntcs(path_to_ntc):
    data = pd.read_excel(path_to_ntc, sheet_name="NTC", skiprows=1)
    data["fromLocation"] = data["Border"].map(
        lambda x: pycountry.countries.lookup(x.split("-")[0][:2]).alpha_3
    )
    data["toLocation"] = data["Border"].map(
        lambda x: pycountry.countries.lookup(x.split("-")[1][:2]).alpha_3
    )
    data["capacity"] = data.min(axis=1) * 1e3 # from MW to # TODO make directional
    data.drop(columns=["Border", "=>", "<="], inplace=True)
    data.drop(index=data[data["fromLocation"] == data["toLocation"]].index, inplace=True)
    return data.groupby(["fromLocation", "toLocation"]).capacity.sum().reset_index()


def _generate_links(ntcs):
    config = AttrDict()
    for _, row in ntcs.iterrows():
        config.union(
            AttrDict.from_yaml_string(
                f"""
                links:
                    {row.fromLocation},{row.toLocation}.techs:
                        ac_transmission:
                            constraints:
                                energy_cap_equals: {row.capacity}
                """
                # TODO allow for transmission capacity extension?
            )
        )
    return config


if __name__ == "__main__":
    generate_links(
        snakemake.input.ntc,
        snakemake.output[0]
    )
