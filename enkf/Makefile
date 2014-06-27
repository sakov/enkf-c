# set MPI status either to "" or "MPI"

ENKF_FLAGS = # -DHE_VIAFILE -DCHECK_X5 -DCHECK_G -DSHUFFLE_ROWS

PROGRAMS =\
bin/enkf_prep\
bin/enkf_calc\
bin/enkf_post

ALL = $(PROGRAMS)

bin/enkf_prep: MPISTATUS =
bin/enkf_calc: MPISTATUS = MPI
bin/enkf_post: MPISTATUS = MPI

# do not edit below

include ./make.inc

PREP_INCS = -I./prep
CALC_INCS = -I./calc
INCS = -I. -I./common -I./models $(INC_NC) $(INC_GU)
LIBS = $(LD_PATH) $(LIB_NC) $(LIB_GU) $(LIB_M)

SRC_COMMON =\
common/ncw.c\
common/stringtable.c\
common/kdtree.c\
common/version.c\
common/global.c\
common/utils.c\
common/enkfprm.c\
common/observations.c\
common/grid.c\
common/model.c\
common/allmodels.c

HDR_COMMON =\
common/nan.h\
common/ncw.h\
common/stringtable.h\
common/kdtree.h\
common/version.h\
common/definitions.h\
common/utils.h\
common/enkfprm.h\
common/observations.h\
common/grid.h\
common/model.h

SRC_PREP =\
prep/prep_utils.c\
prep/obsmeta.c\
prep/allreaders.c\
prep/reader_navo.c\
prep/reader_windsat.c\
prep/reader_rads.c\
prep/reader_cars.c\
prep/reader_mmt.c\
prep/enkf_prep.c

HDR_PREP =\
prep/prep.h\
prep/obsmeta.h\
prep/allreaders.h

SRC_CALC =\
common/distribute.c\
common/pointlog.c\
common/dasystem.c\
calc/obsstats.c\
calc/allhs.c\
calc/model2obs.c\
calc/calcs.c\
calc/ensobs.c\
calc/transforms.c\
calc/enkf_calc.c

HDR_CALC =\
common/distribute.h\
common/allmodels.h\
common/pointlog.h\
common/dasystem.h\
calc/allhs.h\
calc/model2obs.h\
calc/calcs.h

SRC_POST =\
common/distribute.c\
common/pointlog.c\
common/dasystem.c\
post/update.c\
post/enkf_post.c

HDR_POST =\
common/distribute.h\
common/pointlog.h\
common/dasystem.h

SRC_MODEL =\
models/mom4/mom4.c

HDR_MODEL =\
models/mom4.h

HDR_MKUNITOBS =\
common/distribute.h\
common/dasystem.h\
calc/model2obs.h

default: $(PROGRAMS)

all: $(ALL)

bin/enkf_prep: INCS += $(PREP_INCS)
bin/enkf_prep: $(SRC_PREP) $(SRC_COMMON) $(SRC_MODEL) $(HDR_PREP) $(HDR_COMMON) $(HDR_MODEL)
	@echo
	${CC${MPISTATUS}} -o bin/enkf_prep $(ENKF_FLAGS) $(CFLAGS$(MPISTATUS)) $(INCS) $(SRC_PREP) $(SRC_COMMON) $(SRC_MODEL) $(LIBS)

bin/enkf_calc: INCS += $(CALC_INCS)
bin/enkf_calc: LIBS += $(LIB_LAPACK)
bin/enkf_calc: $(SRC_COMMON) $(SRC_MODEL) $(SRC_CALC) $(HDR_COMMON) $(HDR_MODEL) $(HDR_CALC)
	@echo
	${CC${MPISTATUS}} -o bin/enkf_calc $(ENKF_FLAGS) $(CFLAGS$(MPISTATUS)) $(INCS) $(SRC_CALC) $(SRC_COMMON) $(SRC_MODEL) $(LIBS)

bin/enkf_post: INCS += $(POST_INCS)
bin/enkf_post: LIBS += $(LIB_LAPACK)
bin/enkf_post: $(SRC_POST) $(SRC_COMMON) $(SRC_MODEL) $(HDR_POST) $(HDR_COMMON) $(HDR_MODEL)
	@echo
	${CC${MPISTATUS}} -o bin/enkf_post $(ENKF_FLAGS) $(CFLAGS$(MPISTATUS)) $(INCS) $(SRC_POST) $(SRC_COMMON) $(SRC_MODEL) $(LIBS)

clean:
	rm -f bin/*

indent:
	indent -T FILE -T field -T dasystem -T grid -T kdtree -T kdnode -T kdset -T model -T modeldata -T gridprm -T mom4prm -T stringtable -T resnode -T tm -T obsmeta -T obsreader_entry -T H_entry -T model_entry -T enkfprm -T gnll -T rangeentry -T measurement -T obstype -T observations -T gridprm -T pointlog -T region -T calcstats -T obstypedesc -T nc_type -T size_t */*.[ch]
	rm -f `find . -name "*.*~"`
