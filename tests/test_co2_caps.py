from pathlib import Path

import pytest
import calliope
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_CO2_CAPS = PATH_TO_BUILD / "input" / "co2-caps.csv"
PATH_TO_OUTPUT = PATH_TO_BUILD / "output"
FILENAME_RESULTS = Path("results.nc")
EPSILON = 0.001 # 1 t CO2eq
EUROPE_SCENARIOS = ["baseline", "low-cost"]
GERMANY_SCENARIOS = ["baseline-germany", "low-cost-germany"]


@pytest.fixture(scope="module")
def requested_caps():
    return pd.read_csv(PATH_TO_REQUESTED_CO2_CAPS, index_col=0).iloc[:, 0]


@pytest.fixture(scope="session")
def model_output(scenario):
    return calliope.read_netcdf(PATH_TO_OUTPUT / scenario / FILENAME_RESULTS)


@pytest.fixture(scope="module")
def co2_emissions(model_output, variables):
    return (model_output.get_formatted_array("cost")
                        .to_dataframe()
                        .reset_index()
                        .groupby(["locs", "costs"])
                        .sum()
                        .reset_index()
                        .pivot(columns="costs", index="locs", values="cost")
                        .loc[:, "co2"]
                        .mul( 1 / variables["scaling-factors"]["co2"]))


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

    @pytest.fixture(scope="session", params=EUROPE_SCENARIOS)
    def scenario(self, request):
        return request.param


class TestGermanyOnly(Base):

    @pytest.fixture
    def country(self):
        return "DEU"

    @pytest.fixture(scope="session", params=GERMANY_SCENARIOS)
    def scenario(self, request):
        return request.param
