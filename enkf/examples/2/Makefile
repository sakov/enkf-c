NP_CALC = 2
NP_UPDATE = 2

enkf: clean enkf-prep enkf-calc enkf-update

enoi: clean enoi-prep enoi-calc enoi-update

enkf-prep: enkf_prep conf obs prm
	./enkf_prep enkf.prm 2>&1 | tee prep.out

enkf-calc: enkf_calc observations.nc conf conf ensemble_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enkf.prm 2>&1 | tee calc.out

enkf-update: enkf_update X5.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation enkf.prm 2>&1 | tee update.out

enoi-prep: enkf_prep conf obs enoi.prm
	./enkf_prep enoi.prm 2>&1 | tee prep.out

enoi-calc: enkf_calc observations.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm 2>&1 | tee calc.out

enoi-update: enkf_update w.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread enoi.prm 2>&1 | tee update.out

enkf_prep:
	cp ../../bin/enkf_prep .

enkf_calc:
	cp ../../bin/enkf_calc .

enkf_update:
	cp ../../bin/enkf_update .

conf:
	ln -s ../1/conf .
obs:
	ln -s ../1/obs .

ensemble_6565:
	ln -s ../1/ensemble_6565 .

ensemble_static:
	ln -s ../1/ensemble_static .

background_6565:
	ln -s ../1/background_6565 .

prm: enkf.prm enoi.prm model.prm obstypes.prm obs.prm

enkf.prm:
	ln -s ../1/enkf.prm .

enoi.prm:
	ln -s ../1/enoi.prm .

model.prm:
	ln -s ../1/model.prm .

obstypes.prm:
	ln -s ../1/obstypes.prm .

obs.prm:
	ln -s ../1/obs.prm .
clean:
	@rm -f X5*.nc w*.nc enkf_diag*.nc badbatches.txt observations*.nc *.out ensemble_6565/*.analysis* ensemble_6565/*.increment* ensemble_6565/*-???.nc background_6565/*.analysis* background_6565/*.increment* ensemble_6565/ens_*-???.nc background_6565/bg_*-???.nc *.nc enkf_* conf enkf.prm enoi.prm ensemble_6565 model.prm obs obs.prm obstypes.prm
