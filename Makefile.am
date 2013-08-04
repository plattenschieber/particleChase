<<<<<<< HEAD

# This file is part of particleChase.
# Makefile.am in particleChase
# included non-recursively from toplevel directory

bin_PROGRAMS += particleChase/particleChase
particleChase_particleChase_SOURCES = particleChase.c world.c list.c 
LINT_CSOURCES += $(particleChase_particleChase_SOURCES) 
=======
# This file is part of submoduled p4est fork plattenschieber/particleChase
# Makefile.am in particleChase
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += particleChase/particleChase
particleChase_particleChase_SOURCES = particleChase/particleChase.c particleChase/list.c particleChase/world.c 

LINT_CSOURCES += $(particleChase_particleChase_SOURCES)
endif

# only needed if 3D version is available
#if P4EST_ENABLE_BUILD_3D
#bin_PROGRAMS += example/simple/p8est_simple
#example_simple_p8est_simple_SOURCES = example/simple/simple3.c
#
#LINT_CSOURCES += $(example_simple_p8est_simple_SOURCES)
#endif
>>>>>>> origin/p4estIntegration
