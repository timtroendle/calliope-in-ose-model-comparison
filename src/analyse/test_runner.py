import sys

import pytest
import calliope


def run_tests(plugin, path_to_output):
    exit_code = pytest.main(
        [
            f"--html={path_to_output}",
            f"--self-contained-html",
        ],
        plugins=[plugin]
    )
    sys.exit(exit_code)


def _create_config_plugin(scenario_results, scaling_factors):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(
            scope="session",
            params=[scenario for scenario in scenario_results if "germany" not in scenario]
        )
        def europe_model(self, request):
            return calliope.read_netcdf(request.param)

        @pytest.fixture(
            scope="session",
            params=[scenario for scenario in scenario_results if "germany" in scenario]
        )
        def germany_model(self, request):
            return calliope.read_netcdf(request.param)

        @pytest.fixture(scope="session")
        def scaling_factors(self):
            return scaling_factors

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    plugin = _create_config_plugin(
        scenario_results=snakemake.input.results,
        scaling_factors=snakemake.params.scaling_factors
    )
    run_tests(
        plugin=plugin,
        path_to_output=snakemake.output[0]
    )
