NP_CALC = 2
NP_UPDATE = 2

enkf: clean enkf-prep enkf-calc enkf-update

enoi: clean enoi-prep enoi-calc enoi-update

enkf-prep: enkf_prep conf obs prm
	./enkf_prep enkf.prm | tee prep.out

enkf-calc: enkf_calc observations.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enkf.prm | tee calc.out

enkf-update: enkf_update X5-0.nc X5-1.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation enkf.prm | tee update.out

enoi-prep: enkf_prep conf obs prm
	./enkf_prep enoi.prm | tee prep.out

enoi-calc: enkf_calc observations.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm | tee calc.out

enoi-update: enkf_update w-0.nc w-1.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread enoi.prm | tee update.out

enkf_prep:
	cp ../../bin/enkf_prep .

enkf_calc:
	cp ../../bin/enkf_calc .

enkf_update:
	cp ../../bin/enkf_update .

prm: enkf.prm enoi.prm obstypes.prm obs.prm

conf obs ensemble_6565 ensemble_static background_6565 enkf.prm enoi.prm obstypes.prm obs.prm:
	ln -s ../1/$@ .

clean:
	@rm -f X5*.nc w*.nc enkf_diag*.nc badbatches.txt observations*.nc *.out ensemble_6565/*.analysis* ensemble_6565/*.increment* ensemble_6565/*-???.nc background_6565/*.analysis* background_6565/*.increment* ensemble_6565/ens_*-???.nc background_6565/bg_*-???.nc *.nc enkf_* conf enkf.prm enoi.prm ensemble_6565 obs obs.prm obstypes.prm ensemble_6565 ensemble_static background_6565
