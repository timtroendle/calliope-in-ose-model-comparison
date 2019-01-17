from pathlib import Path

import pytest
import calliope
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_RENEWABLE_SHARES = PATH_TO_BUILD / "input" / "renewable-shares.csv"
PATH_TO_MODEL = PATH_TO_BUILD / "output" / "results.nc"
EPSILON = 0.01
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore", "wind_offshore"]


@pytest.fixture(
    scope="module",
    params=pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).index
)
def country(request):
    return request.param


@pytest.fixture(scope="module")
def requested_shares():
    return pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).iloc[:, 0]


@pytest.fixture(scope="module")
def generated_electricity():
    model = calliope.read_netcdf(PATH_TO_MODEL)
    prod = model.get_formatted_array("carrier_prod").to_dataframe(name="carrier_prod")
    prod = prod.groupby(["locs", "techs"]).carrier_prod.sum().reset_index()
    prod.drop(index=prod[prod.techs.str.contains("transmission")].index, inplace=True)
    prod.locs = prod.locs.str[:3]
    prod.techs = prod.techs.map(lambda tech: tech if tech[-4] != "_" else tech[:-4])
    return prod.groupby(["locs", "techs"]).carrier_prod.sum().reset_index().pivot(
        index="locs",
        columns="techs",
        values="carrier_prod"
    )


@pytest.fixture(scope="module")
def consumption():
    model = calliope.read_netcdf(PATH_TO_MODEL)
    con = model.get_formatted_array("carrier_con").to_dataframe(name="carrier_con")
    con = con.groupby(["locs", "techs"]).carrier_con.sum().reset_index()
    con.drop(index=con[con.techs.str.contains("transmission")].index, inplace=True)
    con.locs = con.locs.str[:3]
    return con.groupby(["locs", "techs"]).carrier_con.sum().reset_index().pivot(
        index="locs",
        columns="techs",
        values="carrier_con"
    )


def test_minimal_capacity_is_installed(requested_shares, generated_electricity, consumption, country):
    re_prod = generated_electricity.loc[country, RE_TECHS].sum()
    re_share = re_prod / (consumption.loc[country, "demand_elec"] * -1)

    assert re_share + EPSILON >= requested_shares.loc[country]
