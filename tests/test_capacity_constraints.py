from pathlib import Path

import pytest
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_CAPACITY_CONSTRAINTS = PATH_TO_BUILD / "input" / "capacity.csv"
PATH_TO_OUTPUT_DIRECTORY = PATH_TO_BUILD / "output"
FILENAME_CAPACITY = Path("capacity-raw.csv")
EPSILON = 0.001 # 1 kW
SCENARIOS = ["baseline", "low-cost"]


@pytest.fixture(
    scope="module",
    params=pd.read_csv(PATH_TO_CAPACITY_CONSTRAINTS, index_col=0).index
)
def country(request):
    return request.param


@pytest.fixture(scope="module")
def capacity_constraints():
    return pd.read_csv(PATH_TO_CAPACITY_CONSTRAINTS, index_col=0)


@pytest.fixture(
    scope="module",
    params=SCENARIOS
)
def installed_capacity(request):
    return pd.read_csv(
        PATH_TO_OUTPUT_DIRECTORY / request.param / FILENAME_CAPACITY,
        index_col=0
    ) * 1e3 # from GW to MW


@pytest.mark.parametrize("tech", [
    ("wind_onshore_monopoly"), ("wind_offshore"), ("roof_mounted_pv")
])
def test_minimal_capacity_is_installed(capacity_constraints, installed_capacity, country, tech):
    assert installed_capacity.loc[country, tech] + EPSILON >= capacity_constraints.loc[country, tech]


@pytest.mark.parametrize("tech", [
    ("coal"), ("lignite"), ("ccgt"), ("nuclear")
])
def test_maximum_capacity_is_not_exceeded(capacity_constraints, installed_capacity, country, tech):
    assert installed_capacity.loc[country, tech] - EPSILON <= capacity_constraints.loc[country, tech]
