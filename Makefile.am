# This file is part of submoduled p4est fork plattenschieber/particleChase
# Makefile.am in particleChase
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += particleChase/particleChase
particleChase_particleChase_SOURCES = example/simple/ptest2.c

LINT_CSOURCES += $(particleChase_particleChase_SOURCES)
endif

# only needed if 3D version is available
#if P4EST_ENABLE_BUILD_3D
#bin_PROGRAMS += example/simple/p8est_simple
#example_simple_p8est_simple_SOURCES = example/simple/simple3.c
#
#LINT_CSOURCES += $(example_simple_p8est_simple_SOURCES)
#endif
