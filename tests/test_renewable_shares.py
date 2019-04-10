from pathlib import Path

import pytest
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_RENEWABLE_SHARES = PATH_TO_BUILD / "input" / "renewable-shares.csv"
EPSILON = 0.01
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly",
            "wind_onshore_competing", "wind_offshore", "hydro_run_of_river",
            "hydro_reservoir", "biomass"]


@pytest.fixture(scope="module")
def requested_shares(request):
    return pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).iloc[:, 0]


@pytest.fixture(scope="module")
def renewable_generation(model, variables):
    return (model.get_formatted_array("carrier_prod")
                 .sel(carriers="electricity", techs=RE_TECHS)
                 .sum(dim=["timesteps", "techs"])) / variables["scaling-factors"]["power"]


@pytest.fixture(scope="module")
def demand(model, variables):
    return (
        model.get_formatted_array("carrier_con")
             .sel(carriers="electricity", techs="demand_elec")
             .sum(dim=["timesteps"])
    ) / variables["scaling-factors"]["power"]


@pytest.fixture()
def re_share(renewable_generation, demand, country):
    return renewable_generation.sel(locs=country) / (demand.sel(locs=country) * -1)


@pytest.fixture()
def requested_share(requested_shares, country):
    return requested_shares.loc[country]


class Base:

    def test_renewable_share(self, re_share, requested_share):
        assert re_share + EPSILON >= requested_share


class TestAllEurope(Base):

    @pytest.fixture(
        scope="module",
        params=pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).index
    )
    def country(self, request):
        return request.param

    @pytest.fixture(scope="module")
    def model(self, europe_model):
        return europe_model


class TestGermanyOnly(Base):

    @pytest.fixture
    def country(self):
        return "DEU"

    @pytest.fixture(scope="module")
    def model(self, germany_model):
        return germany_model
