from pathlib import Path

import pytest
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_CAPACITY_CONSTRAINTS = PATH_TO_BUILD / "input" / "capacity.csv"
PATH_TO_INSTALLED_CAPACITY = PATH_TO_BUILD / "output" / "capacity-raw.csv"
EPSILON = 1 # kW


@pytest.fixture(params=pd.read_csv(PATH_TO_CAPACITY_CONSTRAINTS, index_col=0).index)
def country(request):
    return request.param


@pytest.fixture()
def capacity_constraints():
    return pd.read_csv(PATH_TO_CAPACITY_CONSTRAINTS, index_col=0)


@pytest.fixture()
def installed_capacity():
    return pd.read_csv(PATH_TO_INSTALLED_CAPACITY, index_col=0) * 1e6 # from GW to kW


@pytest.mark.parametrize("tech", [
    ("wind_onshore"), ("wind_offshore"), ("roof_mounted_pv")
])
def test_minimal_capacity_is_installed(capacity_constraints, installed_capacity, country, tech):
    assert installed_capacity.loc[country, tech] + EPSILON >= capacity_constraints.loc[country, tech]
