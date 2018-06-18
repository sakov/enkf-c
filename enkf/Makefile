# Select flags common for PREP and CALC from the following list
# (default = none):
#
# -DGRIDNODES_WRITE
# -DINTERNAL_QSORT_R
# -DNO_GRIDUTILS
PREPCALC_FLAGS =

# Select flags for CALC from the following list (default = none):
#
# -DCHECK_X5
# -DCHECK_G
# -DSHUFFLE_ROWS
# -DHE_VIAFILE | -DHE_VIASHMEM
# -DOBS_SHUFFLE
CALC_FLAGS = -DHE_VIASHMEM

# There are currently no flags for UPDATE
#
UPDATE_FLAGS =

# Set MPI status of EnKF-C programs either to "" or "MPI".
# Note: prep is single-threaded only.
#
MPISTATUS_PREP =
MPISTATUS_CALC = MPI
MPISTATUS_UPDATE = MPI

### do not edit below ###

include ./make.inc

INCS = -I. -I./common $(INC_NC) $(INC_GU)
LIBS = $(LD_PATH) $(LIB_NC) $(LIB_GU) $(LIB_M)

SRC_COMMON =\
common/ncw.c\
common/stringtable.c\
common/version.c\
common/global.c\
common/utils.c\
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
common/enkfprm.h\
common/grid.h\
common/gridprm.h\
common/model.h

SRC_PREP =\
common/hash.c\
prep/prep_utils.c\
prep/obsprm.c\
common/obstypes.c\
common/observations.c\
prep/allreaders.c\
prep/reader_navo.c\
prep/reader_windsat.c\
prep/reader_rads.c\
prep/reader_cars.c\
prep/reader_mmt.c\
prep/reader_amsr2.c\
prep/reader_amsre.c\
prep/reader_pathfinder.c\
prep/reader_h8.c\
prep/reader_viirs.c\
prep/reader_xy_scattered.c\
prep/reader_xy_gridded.c\
prep/reader_xyz_scattered.c\
prep/reader_xyz_gridded.c\
prep/reader_xyh_gridded.c\
prep/enkf_prep.c

HDR_PREP =\
prep/prep_utils.h\
common/hash.h\
prep/obsprm.h\
common/obstypes.h\
common/observations.h\
prep/allreaders.h

SRC_CALC =\
common/hash.c\
common/distribute.c\
common/pointlog.c\
common/obstypes.c\
common/observations.c\
common/dasystem.c\
common/kdtree.c\
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
update/enkf_update.c

HDR_UPDATE =\
common/hash.h\
common/distribute.h\
common/pointlog.h\
common/dasystem.h

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
	${CC${MPISTATUS_PREP}} -o bin/enkf_prep -DENKF_PREP $(PREPCALC_FLAGS) $(CFLAGS$(MPISTATUS_PREP)) $(INCS) $(SRC_PREP) $(SRC_COMMON) $(LIBS)

bin/enkf_calc: LIBS := $(LIB_LAPACK) $(LIBS)
bin/enkf_calc: $(SRC_COMMON) $(SRC_CALC) $(HDR_COMMON) $(HDR_CALC)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_CALC}} -o bin/enkf_calc -DENKF_CALC $(PREPCALC_FLAGS) $(CALC_FLAGS) $(CFLAGS$(MPISTATUS_CALC)) $(INCS) $(SRC_CALC) $(SRC_COMMON) $(LIBS)

bin/enkf_update: LIBS := $(LIB_LAPACK) $(LIBS)
bin/enkf_update: $(SRC_UPDATE) $(SRC_COMMON) $(HDR_UPDATE) $(HDR_COMMON)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_UPDATE}} -o bin/enkf_update -DENKF_UPDATE $(UPDATE_FLAGS) $(CFLAGS$(MPISTATUS_UPDATE)) $(INCS) $(SRC_UPDATE) $(SRC_COMMON) $(LIBS)

clean:
	rm -f bin/*

indent:
	indent -T FILE -T field -T dasystem -T grid -T kdtree -T kdnode -T kdset -T model -T modeldata -T gridprm -T gxy_simple -T gxy_curv -T gz_simple -T gz_sigma -T gz_hybrid -T gridnodes -T gridmap -T stringtable -T resnode -T tm -T obsmeta -T metastd -T metapar -T obsreader_entry -T H_entry -T model_entry -T enkfprm -T gnll -T rangeentry -T observation -T obstype -T observations -T pointlog -T region -T obsdomain -T zint -T calcstats -T hashtable -T ht_bucket -T badbatchspec -T nc_type -T size_t -T otdesc -T uint32_t -T uint16_t */*.[ch]
	rm -f `find . -name "*.*~"`
