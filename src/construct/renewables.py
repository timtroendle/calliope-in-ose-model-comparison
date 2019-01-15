"""Preprocess the renewables data to be compatible with Calliope."""
import pandas as pd
import pycountry


def preprocess_renewables_data(path_to_raw_data, assumed_year, technology, path_to_result):
    """Preprocess the renewables data to be compatible with Calliope."""
    if "pv" in technology:
        data = pd.read_excel(path_to_raw_data, sheet_name="pv", index_col=0)
    elif technology == "wind-onshore":
        data = _read_wind(path_to_raw_data, sheet_name="wind onshore")
    elif technology == "wind-offshore":
        data = _read_wind(path_to_raw_data, sheet_name="wind offshore")
    else:
        raise ValueError(f"Unknown technology: {technology}.")
    data.index = pd.date_range(
        start=f"{assumed_year}",
        end=f"{assumed_year + 1}",
        freq="H",
        closed="left"
    )
    data.rename(
        columns=lambda iso2: pycountry.countries.lookup(iso2).alpha_3,
        inplace=True
    )
    data.to_csv(path_to_result, index=True, header=True)


def _read_wind(path_to_raw_data, sheet_name):
    combined = pd.read_excel(path_to_raw_data, sheet_name="wind on&off combined", index_col=0)
    combined["RS"] = 0 # FIXME this is only because Serbia data is missing
    combined["ME"] = 0 # FIXME this is only because Montenegro data is missing
    combined["BA"] = 0 # FIXME this is only because Bosnia and Herzegovina data is missing
    combined["MT"] = 0 # FIXME this is only because Malta data is missing
    combined["AL"] = 0 # FIXME this is only because Albania data is missing
    specific = pd.read_excel(path_to_raw_data, sheet_name=sheet_name, index_col=0)
    specific = specific.reindex(columns=combined.columns)
    return specific.where(~pd.isna(specific), other=combined)


if __name__ == "__main__":
    preprocess_renewables_data(
        path_to_raw_data=snakemake.input.raw,
        assumed_year=snakemake.params.year,
        technology=snakemake.wildcards.technology,
        path_to_result=snakemake.output[0]
    )
