NP_CALC = 2
NP_UPDATE = 2

enoi: clean enoi-prep enoi-calc enoi-update

stats: clean stats-prep stats-calc

enoi-prep: enkf_prep conf obs enoi.prm
	./enkf_prep enoi.prm 2>&1 | tee prep.out

enoi-calc: enkf_calc observations.nc conf background ensemble
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm 2>&1 | tee calc.out

enoi-update: enkf_update conf background ensemble transforms-0.nc transforms-1.nc 
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --output-increment enoi.prm 2>&1 | tee update.out

stats-prep: enkf_prep conf obs stats.prm
	./enkf_prep stats.prm 2>&1 | tee prep.out

stats-calc: enkf_calc observations.nc conf background
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --forecast-stats-only stats.prm 2>&1 | tee calc.out

enkf_prep:
	cp ../../bin/enkf_prep .

enkf_calc:
	cp ../../bin/enkf_calc .

enkf_update:
	cp ../../bin/enkf_update .

clean:
	@rm -f *.nc *.out *.txt background/*.analysis background/*.increment enkf_*
	@rm -rf .enkftmp
