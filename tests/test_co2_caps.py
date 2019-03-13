from pathlib import Path

import pytest
import calliope
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_CO2_CAPS = PATH_TO_BUILD / "input" / "co2-caps.csv"
PATH_TO_OUTPUT = PATH_TO_BUILD / "output"
FILENAME_RESULTS = Path("results.nc")
SCENARIOS = ["baseline", "low-cost", "baseline-germany", "low-cost-germany"]


@pytest.fixture(
    scope="module",
    params=pd.read_csv(PATH_TO_REQUESTED_CO2_CAPS, index_col=0).index
)
def country(request):
    return request.param


@pytest.fixture(scope="module")
def requested_caps():
    return pd.read_csv(PATH_TO_REQUESTED_CO2_CAPS, index_col=0).iloc[:, 0]


@pytest.fixture(
    scope="session",
    params=SCENARIOS
)
def model_output(request):
    return calliope.read_netcdf(PATH_TO_OUTPUT / request.param / FILENAME_RESULTS)


@pytest.fixture(scope="module")
def co2_emissions(model_output):
    return (model_output.get_formatted_array("cost")
                        .to_dataframe()
                        .reset_index()
                        .groupby(["locs", "costs"])
                        .sum()
                        .reset_index()
                        .pivot(columns="costs", index="locs", values="cost")
                        .loc[:, "co2"])


def test_co2_caps(requested_caps, co2_emissions, country):
    assert co2_emissions.loc[country] <= requested_caps.loc[country]
