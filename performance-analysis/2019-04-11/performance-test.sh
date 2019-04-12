#!/bin/sh

mkdir -p ../../build/output/minus6
mkdir -p ../../build/output/minus5
mkdir -p ../../build/output/minus4
mkdir -p ../../build/output/minus3
mkdir -p ../../build/output/minus2

bsub -oo ../../build/output/minus6/01.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus6/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus6/02.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus6/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus6/03.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus6/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-5}"

bsub -oo ../../build/output/minus5/01.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus5/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-5, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus5/02.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus5/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-5, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus5/03.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus5/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-5, run.solver_options.OptimalityTol: 1e-5}"

bsub -oo ../../build/output/minus4/01.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus4/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus4/02.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus4/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus4/03.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus4/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-5}"

bsub -oo ../../build/output/minus3/01.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus3/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-3, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus3/02.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus3/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-3, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus3/03.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus3/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-3, run.solver_options.OptimalityTol: 1e-5}"

bsub -oo ../../build/output/minus2/01.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus2/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus2/02.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus2/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-5}"
bsub -oo ../../build/output/minus2/03.log -W 1439 -n 4 -R "rusage[mem=24000]" calliope run ../../build/model/model.yaml --save_netcdf ../../build/output/minus2/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-5}"
