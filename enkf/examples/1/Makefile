NP_CALC = 2
NP_POST = 2

enoi: clean enoi-prep enoi-calc enoi-post

enkf: clean enkf-prep enkf-calc enkf-post

enoi-prep:
	./enkf_prep enoi.prm 2>&1 | tee prep.out

enoi-calc: observations.nc
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm 2>&1 | tee calc.out

enoi-post: w.nc
	nice -n +19 mpirun -np $(NP_POST) ./enkf_post --separate-output --calculate-spread enoi.prm 2>&1 | tee post.out

enkf-prep: 
	./enkf_prep enkf.prm 2>&1 | tee prep.out

enkf-calc: observations.nc
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enkf.prm 2>&1 | tee calc.out

enkf-post: X5.nc
	nice -n +19 mpirun -np $(NP_POST) ./enkf_post --separate-output --calculate-spread enkf.prm 2>&1 | tee post.out
	rm -f ensemble_6565/mem*_*-???.nc

clean:
	rm -f X5.nc w.nc enkf_diag.nc observations*.nc *.out ensemble_6565/*.analysis ensemble_6565/*.increment background_6565/*.analysis background_6565/*.increment ensemble_6565/ens_*-???.nc background_6565/bg_*-???.nc spread.nc ensemble*/spread*.nc pointlog_*,*.nc
