NP_CALC = 4
NP_UPDATE = 4

enkf: clean enkf-prep enkf-calc enkf-update

enoi: clean enoi-prep enoi-calc enoi-update

enkf-prep: enkf_prep conf obs prm
	./enkf_prep enkf.prm 2>&1 | tee prep.out

enkf-calc: enkf_calc observations.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enkf.prm 2>&1 | tee calc.out

enkf-update: enkf_update transforms.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation enkf.prm 2>&1 | tee update.out

enoi-prep: enkf_prep conf obs prm
	./enkf_prep enoi.prm 2>&1 | tee prep.out

enoi-calc: enkf_calc observations.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm 2>&1 | tee calc.out

enoi-update: enkf_update transforms.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread enoi.prm 2>&1 | tee update.out

enkf_prep:
	cp ../../bin/enkf_prep .

enkf_calc:
	cp ../../bin/enkf_calc .

enkf_update:
	cp ../../bin/enkf_update .

prm: enkf.prm enoi.prm model.prm grid.prm

conf obs ensemble_6565 ensemble_static background_6565 enkf.prm enoi.prm model.prm grid.prm:
	ln -s ../1/$@ .

clean:
	@rm -f transforms*.nc enkf_diag*.nc obsbatches.txt badobsbatches.txt observations*.nc *.out ensemble_6565/*.analysis* ensemble_6565/*.increment* ensemble_6565/*-???.nc background_6565/*.analysis* background_6565/*.increment* ensemble_6565/ens_*-???.nc background_6565/bg_*-???.nc *.nc enkf_* conf enkf.prm enoi.prm model.prm grid.prm ensemble_6565 obs ensemble_6565 ensemble_static background_6565
	@rm -rf .enkftmp
