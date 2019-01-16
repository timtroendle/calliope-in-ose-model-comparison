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
    data.drop(
        columns=[
            "DK_total", # data for sub nodes available, danger of counting data twice
            "IT_total", # data for sub nodes available, danger of counting data twice
            "LU_total", # data for sub nodes available, danger of counting data twice
            "NO_total", # data for sub nodes available, danger of counting data twice
            "SE_total", # data for sub nodes available, danger of counting data twice
            "TR", # Turkey out of scope
            "IL00", # Israel out of scope
            "TN00", # Tunesia out of scope
            "IS00", # Israel out of scope
            "MT", # Malta out of scope
        ],
        inplace=True
    )
    data.rename(
        columns=lambda name: name if name != "NI" else "GB_NI",
        inplace=True
    )
    data.rename(
        columns=lambda name: pycountry.countries.lookup(name.strip()[:2]).alpha_3,
        inplace=True
    )
    data = data.groupby(data.columns, axis='columns').sum()
    data = data * (-1) * 1000 # from MW to kW
    data.to_csv(path_to_result, index=True, header=True)


if __name__ == "__main__":
    preprocess_load_data(
        path_to_raw_load=snakemake.input.raw,
        assumed_year=snakemake.params.year,
        path_to_result=snakemake.output[0]
    )
