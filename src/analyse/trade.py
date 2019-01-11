"""Excavate trade between countries from the model results."""
import calliope
import pandas as pd


def excavate_trade(path_to_results, path_to_trade_amounts):
    """Excavate trade between countries from the model results."""
    model = calliope.read_netcdf(path_to_results)

    imports = _read_imports(model) / 1e6 # from kWh to GWh
    exports = _read_exports(model) / 1e6 # from kWh to GWh

    pd.DataFrame(
        index=imports.index,
        data={
            "import": imports,
            "export": exports.reindex(imports.index),
            "net": imports + exports.reindex(imports.index)
        }
    ).to_csv(
        path_to_trade_amounts,
        index=True,
        header=True,
        float_format="%.1f"
    )


def _read_imports(model):
    data = model.get_formatted_array("carrier_prod").to_dataframe(name="carrier_prod").reset_index()
    data.drop(index=data[~data.techs.str.contains("ac_transmission")].index, inplace=True)
    data.locs = data.locs.str[:3]
    return data.groupby("locs").carrier_prod.sum()


def _read_exports(model):
    data = model.get_formatted_array("carrier_con").to_dataframe(name="carrier_con").reset_index()
    data.drop(index=data[~data.techs.str.contains("ac_transmission")].index, inplace=True)
    data.locs = data.locs.str[:3]
    return data.groupby("locs").carrier_con.sum()


if __name__ == "__main__":
    excavate_trade(
        path_to_results=snakemake.input.results[0],
        path_to_trade_amounts=snakemake.output[0]
    )
