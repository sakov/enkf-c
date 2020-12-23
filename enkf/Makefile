# Select flags common for PREP and CALC from the following list
# (default = none):
#
# -DINTERNAL_QSORT_R
PREPCALC_FLAGS =

# Select flags for CALC from the following list (default = none):
#
# -DSHUFFLE_ROWS
# -DHE_VIASHMEM
# -DOBS_SHUFFLE
# -DMINIMISE_ALLOC
# -DUSE_SHMEM
# -DTW_VIAFILE
# -DDEFLATE_ALL
CALC_FLAGS = -DUSE_SHMEM -DMINIMISE_ALLOC

# Select flags for UPDATE
#
# -DNCW_SKIPSINGLE
# -DDEFLATE_ALL
UPDATE_FLAGS =

# Set MPI status of EnKF-C programs either to "" or "MPI".
# Note: prep is single-threaded only.
#
MPISTATUS_PREP =
MPISTATUS_CALC = MPI
MPISTATUS_UPDATE = MPI

### do not edit below ###

include ./make.inc

INCS_COMMON = -I. -I./common $(INC_NC)
LIBS_COMMON = $(LD_PATH) $(LIB_NC) $(LIB_M)

INCS_PREP =
LIBS_PREP =
INCS_CALC =
LIBS_CALC = $(LIB_LAPACK)
INCS_UPDATE =
LIBS_UPDATE = $(LIB_LAPACK)

SRC_COMMON =\
common/ncw.c\
common/stringtable.c\
common/version.c\
common/global.c\
common/utils.c\
common/ncutils.c\
common/enkfprm.c\
common/grid.c\
common/gridprm.c\
common/model.c\

HDR_COMMON =\
common/ncw.h\
common/stringtable.h\
common/version.h\
common/definitions.h\
common/utils.h\
common/ncutils.h\
common/enkfprm.h\
common/grid.h\
common/gridprm.h\
common/model.h

SRC_PREP =\
common/hash.c\
common/gxy_curv.c\
prep/prep_utils.c\
prep/obsprm.c\
common/kdtree.c\
common/obstypes.c\
common/observations.c\
prep/allreaders.c\
prep/reader_navo.c\
prep/reader_windsat.c\
prep/reader_cars.c\
prep/reader_mmt.c\
prep/reader_amsr2.c\
prep/reader_amsre.c\
prep/reader_h8.c\
prep/reader_scattered.c\
prep/reader_xy_gridded.c\
prep/reader_xyz_gridded.c\
prep/reader_xyh_gridded.c\
prep/reader_cmems.c\
prep/reader_z.c\
prep/enkf_prep.c

HDR_PREP =\
common/hash.h\
common/kdtree.h\
common/gxy_curv.h\
prep/prep_utils.h\
prep/obsprm.h\
common/obstypes.h\
common/observations.h\
prep/allreaders.h

SRC_CALC =\
common/hash.c\
common/distribute.c\
common/gxy_curv.c\
common/kdtree.c\
common/pointlog.c\
common/obstypes.c\
common/observations.c\
common/dasystem.c\
calc/obsstats.c\
calc/allhs.c\
calc/model2obs.c\
calc/calcs.c\
calc/ensobs.c\
calc/transforms.c\
calc/enkf_calc.c

HDR_CALC =\
common/hash.h\
common/distribute.h\
common/gxy_curv.h\
common/pointlog.h\
common/dasystem.h\
common/kdtree.h\
common/obstypes.h\
common/observations.h\
calc/allhs.h\
calc/model2obs.h\
calc/calcs.h

SRC_UPDATE =\
common/hash.c\
common/distribute.c\
common/pointlog.c\
common/dasystem.c\
update/update.c\
update/diags.c\
update/enkf_update.c

HDR_UPDATE =\
common/hash.h\
common/distribute.h\
common/pointlog.h\
common/dasystem.h\
update/diags.h

PROGRAMS =\
bin/enkf_prep\
bin/enkf_calc\
bin/enkf_update

ALL = $(PROGRAMS)

default: $(PROGRAMS)

all: $(ALL)

bin/enkf_prep: $(SRC_PREP) $(SRC_COMMON) $(HDR_PREP) $(HDR_COMMON)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_PREP}} -o bin/enkf_prep -DENKF_PREP $(PREPCALC_FLAGS) $(CFLAGS$(MPISTATUS_PREP)) $(INCS_PREP) $(INCS_COMMON) $(SRC_PREP) $(SRC_COMMON) $(LIBS_PREP) $(LIBS_COMMON)

bin/enkf_calc: $(SRC_COMMON) $(SRC_CALC) $(HDR_COMMON) $(HDR_CALC)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_CALC}} -o bin/enkf_calc -DENKF_CALC $(PREPCALC_FLAGS) $(CALC_FLAGS) $(CFLAGS$(MPISTATUS_CALC)) $(INCS_CALC) $(INCS_COMMON) $(SRC_CALC) $(SRC_COMMON) $(LIBS_CALC) $(LIBS_COMMON)

bin/enkf_update: $(SRC_UPDATE) $(SRC_COMMON) $(HDR_UPDATE) $(HDR_COMMON)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_UPDATE}} -o bin/enkf_update -DENKF_UPDATE $(UPDATE_FLAGS) $(CFLAGS$(MPISTATUS_UPDATE)) $(INCS_UPDATE) $(INCS_COMMON) $(SRC_UPDATE) $(SRC_COMMON) $(LIBS_UPDATE) $(LIBS_COMMON)

clean:
	rm -f bin/*

indent:
	indent -T FILE -T field -T dasystem -T model -T variable -T kdtree -T kdnode -T kdset -T resnode -T modeldata -T grid -T gridprm -T gxy_simple -T gxy_curv -T gz_simple -T gz_sigma -T gz_hybrid -T gxy_curv -T stringtable -T tm -T obsmeta -T metastd -T metapar -T obsreader_entry -T H_entry -T model_entry -T enkfprm -T gnll -T rangeentry -T observation -T obstype -T observations -T pointlog -T region -T obsdomain -T zint -T calcstats -T hashtable -T ht_bucket -T badbatchspec -T nctype2str -T nc_type -T size_t -T otdesc -T uint32_t -T uint16_t */*.[ch]
	rm -f `find . -name "*.*~"`

# following are 3 convenience targets for compiling particular binaries

.PHONY: prep
prep:
	@rm -f bin/enkf_prep
	@$(MAKE) --no-print-directory bin/enkf_prep

.PHONY: calc
calc:
	@rm -f bin/enkf_calc
	@$(MAKE) --no-print-directory bin/enkf_calc

.PHONY: update
update:
	@rm -f bin/enkf_update
	@$(MAKE) --no-print-directory bin/enkf_update
