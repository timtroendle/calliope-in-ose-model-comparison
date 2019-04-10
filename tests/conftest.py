import json

import calliope
import pytest

# There is no other way to parameterise a fixture from a file than hard coding
# the file here.
VARIABLES_FILE = "variables.json"

with open(VARIABLES_FILE) as variables_file:
    scenarios = json.load(variables_file)["scenarios"]


@pytest.fixture(scope="session", params=[scenario for scenario in scenarios if "germany" not in scenario])
def europe_model(request):
    return calliope.read_netcdf(request.param)


@pytest.fixture(scope="session", params=[scenario for scenario in scenarios if "germany" in scenario])
def germany_model(request):
    return calliope.read_netcdf(request.param)
