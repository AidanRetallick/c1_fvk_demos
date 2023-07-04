# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS = clamped_square_inflation \
		  rotated_square \
		  circular_disc \
		  circular_sector
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Local sources that each code depends on:
clamped_square_inflation_SOURCES = \
 clamped_square_inflation.cc
rotated_square_SOURCES = \
 rotated_square.cc
circular_disc_SOURCES = \
 circular_disc.cc
circular_sector_SOURCES = \
 circular_sector.cc
#---------------------------------------------------------------------------

clamped_square_inflation_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

rotated_square_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_disc_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_sector_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
#jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
# jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.

AM_CPPFLAGS += -I@includedir@ -Wall -ggdb3 -Og
#EXTRA_DIST = Circle1.1.ele Circle1.1.node Circle1.1.poly
