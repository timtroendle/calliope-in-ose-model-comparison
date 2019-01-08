"""Preprocess the renewables data to be compatible with Calliope."""
import pandas as pd
import pycountry


def preprocess_load_data(path_to_raw_load, assumed_year, path_to_result):
    """Preprocess the load data to be compatible with Calliope."""
    data = pd.read_excel(path_to_raw_load, sheet_name="load (MW)", index_col=0)
    data.index = pd.date_range(
        start=f"{assumed_year}",
        end=f"{assumed_year + 1}",
        freq="H",
        closed="left"
    )
    data["FR"] = data["FR"] + data["FR15"] # france is lacking a "total" column
    # for all other countries, choose "total" header and delete all subnodes with name
    # longer than two letters (e.g. DKw)
    data.rename(
        columns=lambda name: name if "_total" not in name else name.split("_total")[0],
        inplace=True
    )
    data.drop(
        columns=[col for col in data.columns if len(col) > 2],
        inplace=True
    )
    data.rename(
        columns=lambda iso2: pycountry.countries.lookup(iso2).alpha_3,
        inplace=True
    )
    data = data * (-1) * 1000 # from MW to kW
    data.to_csv(path_to_result, index=True, header=True)


if __name__ == "__main__":
    preprocess_load_data(
        path_to_raw_load=snakemake.input.raw,
        assumed_year=snakemake.params.year,
        path_to_result=snakemake.output[0]
    )
