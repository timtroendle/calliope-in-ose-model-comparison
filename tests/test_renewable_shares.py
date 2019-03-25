from pathlib import Path

import pytest
import calliope
import pandas as pd

PATH_TO_BUILD = Path(__file__).parent / ".." / "build"
PATH_TO_REQUESTED_RENEWABLE_SHARES = PATH_TO_BUILD / "input" / "renewable-shares.csv"
PATH_TO_OUTPUT = PATH_TO_BUILD / "output"
FILENAME_RESULTS = Path("results.nc")
EPSILON = 0.01
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly",
            "wind_onshore_competing", "wind_offshore", "hydro_run_of_river", "biomass"]
EUROPE_SCENARIOS = ["baseline", "low-cost"]
GERMANY_SCENARIOS = ["baseline-germany", "low-cost-germany"]


@pytest.fixture(scope="module")
def requested_shares(request):
    return pd.read_csv(PATH_TO_REQUESTED_RENEWABLE_SHARES, index_col=0).iloc[:, 0]


@pytest.fixture(scope="session")
def model_output(scenario):
    return calliope.read_netcdf(PATH_TO_OUTPUT / scenario / FILENAME_RESULTS)


@pytest.fixture(scope="module")
def generated_electricity(model_output, variables):
    prod = model_output.get_formatted_array("carrier_prod").to_dataframe(name="carrier_prod")
    prod = prod.groupby(["locs", "techs"]).carrier_prod.sum().reset_index()
    prod.drop(index=prod[prod.techs.str.contains("transmission")].index, inplace=True)
    prod.locs = prod.locs.str[:3]
    prod.techs = prod.techs.map(lambda tech: tech if tech[-4] != "_" else tech[:-4])
    return (prod.groupby(["locs", "techs"])
                .carrier_prod
                .sum()
                .reset_index()
                .pivot(index="locs", columns="techs", values="carrier_prod")
                .mul(1 / variables["scaling-factors"]["power"]))


@pytest.fixture(scope="module")
def consumption(model_output, variables):
    con = model_output.get_formatted_array("carrier_con").to_dataframe(name="carrier_con")
    con = con.groupby(["locs", "techs"]).carrier_con.sum().reset_index()
    con.drop(index=con[con.techs.str.contains("transmission")].index, inplace=True)
    con.locs = con.locs.str[:3]
    return (con.groupby(["locs", "techs"])
               .carrier_con
               .sum()
               .reset_index()
               .pivot(index="locs", columns="techs", values="carrier_con")
               .mul(1 / variables["scaling-factors"]["power"]))


@pytest.fixture()
def re_share(generated_electricity, consumption, country):
    re_prod = generated_electricity.loc[country, RE_TECHS].sum()
    return re_prod / (consumption.loc[country, "demand_elec"] * -1)


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
