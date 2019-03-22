#!/bin/sh
bsub -oo build/output/minus6/01.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/02.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/03.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/04.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/04.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/05.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/05.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/06.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/06.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/07.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/07.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/08.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/08.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/09.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/09.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"
bsub -oo build/output/minus6/10.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus6/10.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-6, run.solver_options.OptimalityTol: 1e-6}"

bsub -oo build/output/minus4/01.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/02.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/03.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/04.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/04.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/05.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/05.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/06.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/06.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/07.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/07.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/08.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/08.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/09.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/09.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"
bsub -oo build/output/minus4/10.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus4/10.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-4, run.solver_options.OptimalityTol: 1e-4}"

bsub -oo build/output/minus2/01.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/01.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/02.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/02.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/03.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/03.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/04.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/04.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/05.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/05.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/06.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/06.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/07.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/07.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/08.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/08.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/09.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/09.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
bsub -oo build/output/minus2/10.log -W 1439 -n 4 -R "rusage[mem=12000]" -B -N calliope run build/model/model.yaml --save_netcdf build/output/minus2/10.nc --scenario baseline --override_dict "{run.solver_options.FeasibilityTol: 1e-2, run.solver_options.OptimalityTol: 1e-2}"
