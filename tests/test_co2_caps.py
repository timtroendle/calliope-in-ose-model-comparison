from pathlib import Path

import pytest
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_CO2_CAPS = PATH_TO_BUILD / "input" / "co2-caps.csv"
PATH_TO_OUTPUT = PATH_TO_BUILD / "output"
FILENAME_RESULTS = Path("results.nc")
EPSILON = 0.001 # 1 t CO2eq


@pytest.fixture(scope="module")
def requested_caps():
    return pd.read_csv(PATH_TO_REQUESTED_CO2_CAPS, index_col=0).iloc[:, 0]


@pytest.fixture(scope="module")
def co2_emissions(model, scaling_factors):
    return (model.get_formatted_array("cost")
                 .sel(costs="co2")
                 .sum(dim="techs")) / scaling_factors["co2"]


class Base:

    def test_co2_caps(self, requested_caps, co2_emissions, country):
        assert co2_emissions.loc[country] <= requested_caps.loc[country] + EPSILON


class TestAllEurope(Base):

    @pytest.fixture(
        scope="module",
        params=pd.read_csv(PATH_TO_REQUESTED_CO2_CAPS, index_col=0).index
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
