# The compile flags are described in section 2.16.1 of EnKF-C user guide
# available from enkf-c/enkf/doc/enkf-userguide.pdf.

# Select flags common for PREP and CALC from the following list
# (default = none):
#
# -DINTERNAL_QSORT_R
PREPCALC_FLAGS =

# Select flags for CALC from the following list (default = none):
#
# -DSHUFFLE_ROWS
# -DOBS_SHUFFLE
# -DMINIMISE_ALLOC
# -DUSE_SHMEM
# -DTW_VIAFILE
# -DDEFLATE_ALL
# -DUSE_MPIQUEUE
#
# NOTE: When using just a few CPUs remove "-DUSE_MPIQUEUE". (This will also
#       remove automatic definition of -DTW_VIAFILE in definitions.h.)
#
CALC_FLAGS = -DUSE_SHMEM -DMINIMISE_ALLOC -DUSE_MPIQUEUE

# Select flags for UPDATE
#
# -DNCW_SKIPSINGLE
# -DDEFLATE_ALL
# -DUSE_MPIQUEUE
#
# NOTE: when using just a few CPUs remove "-DUSE_MPIQUEUE". (This will also
#       remove automatic definition of -DTW_VIAFILE in definitions.h.)
#
UPDATE_FLAGS = -DUSE_MPIQUEUE

DIAG_FLAGS =

# Set MPI status of EnKF-C programs either to "" or "MPI".
# Note: prep is single-threaded only.
#
MPISTATUS_PREP =
MPISTATUS_CALC = MPI
MPISTATUS_UPDATE = MPI
MPISTATUS_DIAG = MPI

### do not edit below ###

include ./make.inc

INCS_COMMON = -I./common -I./grid
LIBS_COMMON = $(LD_PATH) $(LIB_NC) $(LIB_M)

LIBS_PREP =
LIBS_CALC = $(LIB_LAPACK)
LIBS_UPDATE = $(LIB_LAPACK)
LIBS_DIAG =

SRC_COMMON =\
common/ncw.c\
common/stringtable.c\
common/version.c\
common/global.c\
common/utils.c\
common/ncutils.c\
common/enkfprm.c\
common/model.c\
grid/grid.c\
grid/gridprm.c\

HDR_COMMON =\
common/ncw.h\
common/stringtable.h\
common/version.h\
common/definitions.h\
common/utils.h\
common/ncutils.h\
common/enkfprm.h\
common/model.h\
grid/grid.h\
grid/gridprm.h

SRC_PREP =\
common/hash.c\
common/kdtree.c\
common/obstypes.c\
common/observations.c\
grid/triangulation.c\
grid/gxy_rect.c\
grid/gxy_curv.c\
grid/gxy_curv2.c\
grid/gxy_unstr.c\
grid/vgrid.c\
grid/hgrid.c\
prep/prep_utils.c\
prep/obsprm.c\
prep/allreaders.c\
prep/reader_navo.c\
prep/reader_windsat.c\
prep/reader_mmt.c\
prep/reader_amsr2.c\
prep/reader_amsre.c\
prep/reader_scattered.c\
prep/reader_gridded_xy.c\
prep/reader_gridded_xyz.c\
prep/reader_gridded_xyh.c\
prep/reader_cmems.c\
prep/reader_z.c\
prep/reader_en4.c\
prep/enkf_prep.c

HDR_PREP =\
common/hash.h\
common/kdtree.h\
common/obstypes.h\
common/observations.h\
grid/triangulation.h\
grid/gxy_rect.h\
grid/gxy_curv.h\
grid/gxy_curv2.h\
grid/gxy_unstr.h\
grid/vgrid.h\
grid/hgrid.h\
prep/prep_utils.h\
prep/obsprm.h\
prep/allreaders.h

SRC_CALC =\
common/hash.c\
common/distribute.c\
common/mpiqueue.c\
common/kdtree.c\
common/pointlog.c\
common/obstypes.c\
common/observations.c\
common/dasystem.c\
grid/triangulation.c\
grid/gxy_unstr.c\
grid/gxy_rect.c\
grid/gxy_curv.c\
grid/gxy_curv2.c\
grid/vgrid.c\
grid/hgrid.c\
calc/obsstats.c\
calc/allhs.c\
calc/model2obs.c\
calc/calcs.c\
calc/ensobs.c\
calc/transforms.c\
calc/enkf_calc.c

HDR_CALC =\
common/hash.h\
common/kdtree.h\
common/distribute.h\
common/mpiqueue.h\
common/pointlog.h\
common/dasystem.h\
common/obstypes.h\
common/observations.h\
grid/triangulation.h\
grid/gxy_rect.h\
grid/gxy_curv.h\
grid/gxy_curv2.h\
grid/gxy_unstr.h\
grid/vgrid.h\
grid/hgrid.h\
calc/allhs.h\
calc/model2obs.h\
calc/calcs.h

SRC_UPDATE =\
common/hash.c\
common/distribute.c\
common/mpiqueue.c\
common/pointlog.c\
common/dasystem.c\
grid/gxy_rect.c\
grid/vgrid.c\
grid/hgrid.c\
update/update.c\
diag/diags.c\
update/enkf_update.c

HDR_UPDATE =\
common/hash.h\
common/distribute.h\
common/mpiqueue.h\
grid/vgrid.h\
grid/hgrid.h\
common/pointlog.h\
common/dasystem.h\
diag/diags.h

SRC_DIAG =\
common/hash.c\
common/distribute.c\
common/dasystem.c\
grid/gxy_rect.c\
grid/vgrid.c\
grid/hgrid.c\
diag/diags.c\
diag/ens_diag.c

HDR_DIAG =\
common/hash.h\
common/distribute.h\
grid/vgrid.h\
grid/hgrid.h\
common/dasystem.h\
diag/diags.h

PROGRAMS =\
bin/enkf_prep\
bin/enkf_calc\
bin/enkf_update\
bin/ens_diag

ALL = $(PROGRAMS)

VERSION := $(shell grep ENKF_VERSION common/version.c | awk -F '"' '{print $$2}')

default: $(PROGRAMS)

all: $(ALL)

bin/enkf_prep: $(SRC_PREP) $(SRC_COMMON) $(HDR_PREP) $(HDR_COMMON)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_PREP}} -o bin/enkf_prep -DENKF_PREP $(PREPCALC_FLAGS) $(CFLAGS$(MPISTATUS_PREP)) $(INCS_COMMON) $(SRC_PREP) $(SRC_COMMON) $(LIBS_PREP) $(LIBS_COMMON)

bin/enkf_calc: $(SRC_COMMON) $(SRC_CALC) $(HDR_COMMON) $(HDR_CALC)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_CALC}} -o bin/enkf_calc -DENKF_CALC $(PREPCALC_FLAGS) $(CALC_FLAGS) $(CFLAGS$(MPISTATUS_CALC)) $(INCS_COMMON) $(SRC_CALC) $(SRC_COMMON) $(LIBS_CALC) $(LIBS_COMMON)

bin/enkf_update: $(SRC_UPDATE) $(SRC_COMMON) $(HDR_UPDATE) $(HDR_COMMON)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_UPDATE}} -o bin/enkf_update -DENKF_UPDATE $(UPDATE_FLAGS) $(CFLAGS$(MPISTATUS_UPDATE)) $(INCS_COMMON) $(SRC_UPDATE) $(SRC_COMMON) $(LIBS_UPDATE) $(LIBS_COMMON)

bin/ens_diag: $(SRC_DIAG) $(SRC_COMMON) $(HDR_DIAG) $(HDR_COMMON)
	@echo
	@mkdir -p bin
	${CC${MPISTATUS_DIAG}} -o bin/ens_diag -DENS_DIAG $(DIAG_FLAGS) $(CFLAGS$(MPISTATUS_UPDATE)) $(INCS_COMMON) $(SRC_DIAG) $(SRC_COMMON) $(LIBS_DIAG) $(LIBS_COMMON)

clean:
	rm -f bin/*

indent:
	indent -T FILE -T uint32_t -T uint16_t -T int8_t -T int16_t -T int32_t -T int64_t -T badbatchspec -T region -T obsregion -T zint -T field -T dasystem -T model -T variable -T modeldata -T gridprm -T grid -T vgrid -T hgrid -T gxy_1d -T gxy_2d -T gxy_rect -T gxy_curv -T gxy_curv2 -T gxy_unstr -T point -T triangle -T triangle_neighbours -T triangulation -T gz_z -T gz_sigma -T gz_hybrid -T gz_numeric -T obsmeta -T metastd -T metapar -T obsreader_entry -T H_entry -T model_entry -T enkfprm -T observation -T obstype -T observations -T otdesc -T pointlog -T calcstats -T stringtable -T kdtree -T kdnode -T kdresult -T hashtable -T ht_bucket -T nctype2str -T nc_type -T size_t -T tm -T mpiqueue */*.[ch]
	rm -f `find . -name "*.*~"`

tar:
	cd ..; tar --exclude=examples --exclude=bin -czvf enkf-c-v$(VERSION).tar.gz enkf; echo "  ->../enkf-c-v$(VERSION).tar.gz"

# following are 4 convenience targets for compiling particular binaries

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

.PHONY: diag
diag:
	@rm -f bin/ens_diag
	@$(MAKE) --no-print-directory bin/ens_diag
