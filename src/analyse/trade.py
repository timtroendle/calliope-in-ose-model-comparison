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
    data = model.results["carrier_prod"].to_dataframe().reset_index()

    data["location"] = data["loc_tech_carriers_prod"].map(lambda x: x.split("::")[0])
    data["tech"] = data["loc_tech_carriers_prod"].map(lambda x: x.split("::")[1])

    data.drop(index=data[~data.tech.str.contains("ac_transmission")].index, inplace=True)

    return data.groupby("location").carrier_prod.sum()


def _read_exports(model):
    data = model.results["carrier_con"].to_dataframe().reset_index()

    data["location"] = data["loc_tech_carriers_con"].map(lambda x: x.split("::")[0])
    data["tech"] = data["loc_tech_carriers_con"].map(lambda x: x.split("::")[1])

    data.drop(index=data[~data.tech.str.contains("ac_transmission")].index, inplace=True)

    return data.groupby("location").carrier_con.sum()



if __name__ == "__main__":
    excavate_trade(
        path_to_results=snakemake.input.results[0],
        path_to_trade_amounts=snakemake.output[0]
    )
