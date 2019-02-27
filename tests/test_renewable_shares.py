from pathlib import Path

import pytest
import calliope
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_RENEWABLE_SHARES = PATH_TO_BUILD / "input" / "renewable-shares.csv"
PATH_TO_OUTPUT = PATH_TO_BUILD / "output"
FILENAME_RESULTS = Path("results.nc")
EPSILON = 0.01
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore"]
NON_RE_TECHS = ["coal", "lignite", "ccgt", "nuclear"]
SCENARIOS = ["baseline", "low-cost"]


@pytest.fixture(
    scope="module",
    params=pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).index
)
def country(request):
    return request.param


@pytest.fixture(scope="module")
def requested_shares():
    return pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).iloc[:, 0]


@pytest.fixture(
    scope="session",
    params=SCENARIOS
)
def model_output(request):
    return calliope.read_netcdf(PATH_TO_OUTPUT / request.param / FILENAME_RESULTS)


@pytest.fixture(scope="module")
def generated_electricity(model_output):
    prod = model_output.get_formatted_array("carrier_prod").to_dataframe(name="carrier_prod")
    prod = prod.groupby(["locs", "techs"]).carrier_prod.sum().reset_index()
    prod.drop(index=prod[prod.techs.str.contains("transmission")].index, inplace=True)
    prod.locs = prod.locs.str[:3]
    prod.techs = prod.techs.map(lambda tech: tech if tech[-4] != "_" else tech[:-4])
    return prod.groupby(["locs", "techs"]).carrier_prod.sum().reset_index().pivot(
        index="locs",
        columns="techs",
        values="carrier_prod"
    )


@pytest.mark.xfail(reason="Not yet implemented.")
def test_minimal_capacity_is_installed(requested_shares, generated_electricity, country):
    re_prod = generated_electricity.loc[country, RE_TECHS].sum()
    non_re_prod = generated_electricity.loc[country, NON_RE_TECHS].sum()
    re_share = re_prod / (non_re_prod + re_prod)

    assert re_share + EPSILON >= requested_shares.loc[country]
