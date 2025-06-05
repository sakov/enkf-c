NP_CALC = 4
NP_UPDATE = 4

enkf: clean enkf-prep enkf-calc enkf-update

enoi: clean enoi-prep enoi-calc enoi-update

stats: clean stats-prep stats-calc

singleob: clean singleob-calc singleob-update

pointlogs: pointlogs-calc pointlogs-update

global: clean enkf-prep enkf-calc-global enkf-update-global

hybrid: clean hybrid-prep hybrid-calc hybrid-update

enkf-prep: enkf_prep conf obs enkf.prm
	./enkf_prep enkf.prm 2>&1 | tee prep.out

enkf-calc: enkf_calc observations.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enkf.prm 2>&1 | tee calc.out

enkf-update: enkf_update transforms.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation --output-increment enkf.prm 2>&1 | tee update.out

enoi-prep: enkf_prep conf obs enoi.prm
	./enkf_prep enoi.prm 2>&1 | tee prep.out

enoi-calc: enkf_calc observations.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm 2>&1 | tee calc.out

enoi-update: enkf_update transforms.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --output-increment enoi.prm 2>&1 | tee update.out

stats-prep: enkf_prep conf obs stats.prm
	./enkf_prep stats.prm 2>&1 | tee prep.out

stats-calc: enkf_calc observations.nc conf background_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --forecast-stats-only stats.prm 2>&1 | tee calc.out

singleob-calc: enkf_calc ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --single-observation 156.0 -32.0 0 SST 0.5 0.5 singleob.prm 2>&1 | tee calc.out

singleob-update: enkf_update ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --output-increment singleob.prm 2>&1 | tee update.out

pointlogs-calc: enkf_calc observations.nc conf ensemble_6565 transforms.nc
	cp enkf.prm pointlogs.prm
	echo "POINTLOG 160 -29" >> pointlogs.prm
	echo "POINTLOG 152 -36" >> pointlogs.prm
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --point-logs-only pointlogs.prm 2>&1 | tee calc-pointlogs.out

pointlogs-update: enkf_update transforms.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --no-fields-write pointlogs.prm 2>&1 | tee update-pointlogs.out

enkf-calc-global: enkf_calc observations.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enkf-global.prm 2>&1 | tee calc.out

enkf-update-global: enkf_update transforms.nc conf ensemble_6565
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation --output-increment enkf-global.prm 2>&1 | tee update.out
	rm -f ensemble_6565/mem*_*-???.nc

hybrid-prep: enkf_prep conf obs hybrid.prm
	./enkf_prep enkf.prm 2>&1 | tee prep.out

hybrid-calc: enkf_calc observations.nc conf ensemble_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc hybrid.prm 2>&1 | tee calc.out

hybrid-update: enkf_update transforms.nc conf ensemble_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation --output-increment hybrid.prm 2>&1 | tee update.out

enkf_prep:
	cp ../../bin/enkf_prep .

enkf_calc:
	cp ../../bin/enkf_calc .

enkf_update:
	cp ../../bin/enkf_update .

clean:
	@rm -f transforms*.nc w*.nc enkf_diag*.nc badobsbatches.txt obsbatches.txt observations*.nc *.out ensemble_6565/*.analysis* ensemble_6565/*.increment* ensemble_6565/*-???.nc background_6565/*.analysis* background_6565/*.increment* ensemble_6565/ens_*-???.nc background_6565/bg_*-???.nc *.nc enkf_* pointlogs.prm
	@rm -rf .enkftmp
