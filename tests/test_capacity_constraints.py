from pathlib import Path

import pytest
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_CAPACITY_CONSTRAINTS = PATH_TO_BUILD / "input" / "capacity.csv"
PATH_TO_OUTPUT_DIRECTORY = PATH_TO_BUILD / "output"
FILENAME_CAPACITY = Path("capacity-raw.csv")
EPSILON = 0.001 # 1 kW
EUROPE_SCENARIOS = ["baseline", "low-cost", "lowest-cost"]
GERMANY_SCENARIOS = ["baseline-germany", "low-cost-germany", "lowest-cost-germany"]


@pytest.fixture(scope="module")
def capacity_constraints():
    return pd.read_csv(PATH_TO_CAPACITY_CONSTRAINTS, index_col=0)


@pytest.fixture(scope="module")
def installed_capacity(scenario):
    return pd.read_csv(
        PATH_TO_OUTPUT_DIRECTORY / scenario / FILENAME_CAPACITY,
        index_col=0
    ) * 1e3 # from GW to MW


class Base:

    @pytest.mark.parametrize("tech", [
        ("wind_onshore_monopoly"), ("wind_offshore"), ("roof_mounted_pv")
    ])
    def test_minimal_capacity_is_installed(self, capacity_constraints, installed_capacity, country, tech):
        assert installed_capacity.loc[country, tech] + EPSILON >= capacity_constraints.loc[country, tech]

    @pytest.mark.parametrize("tech", [
        ("hydro_run_of_river"), ("hydro_reservoir"), ("biomass"), ("pumped_hydro")
    ])
    def test_exact_capacity_is_installed(self, capacity_constraints, installed_capacity, country, tech):
        assert installed_capacity.loc[country, tech] == pytest.approx(
            capacity_constraints.loc[country, tech],
            abs=EPSILON
        )

    @pytest.mark.parametrize("tech", [
        ("coal"), ("lignite"), ("ccgt"), ("nuclear")
    ])
    def test_maximum_capacity_is_not_exceeded(self, capacity_constraints, installed_capacity, country, tech):
        assert installed_capacity.loc[country, tech] - EPSILON <= capacity_constraints.loc[country, tech]


class TestAllEurope(Base):

    @pytest.fixture(
        scope="module",
        params=pd.read_csv(PATH_TO_CAPACITY_CONSTRAINTS, index_col=0).index
    )
    def country(self, request):
        return request.param

    @pytest.fixture(scope="module", params=EUROPE_SCENARIOS)
    def scenario(self, request):
        return request.param


class TestGermanyOnly(Base):

    @pytest.fixture
    def country(self):
        return "DEU"

    @pytest.fixture(scope="module", params=GERMANY_SCENARIOS)
    def scenario(self, request):
        return request.param
