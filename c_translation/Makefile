# Makefile for the C port of fmm2d.
#
# This is a standalone makefile, separate from the top-level fmm2d
# makefile. It builds the C translations of selected Fortran sources
# and the per-file differential tests that compare each C routine to
# its Fortran original.
#
# Usage:
#   make             - build all C objects and all diff tests
#   make tests       - build and run all diff tests
#   make clean       - remove build artifacts
#
# Most diff tests pin BOTH the C translation under test and the
# corresponding Fortran original at -O0. This is necessary because
# gfortran -O3 and gcc -O3 do not produce bit-identical floating-point
# results for complex arithmetic — the optimizers reorder operations
# (vectorization, FMA, etc.) in subtly different ways. At -O0 the two
# compilers do strictly the operations the source code requests, in
# the order it requests them, and we get exact bit-for-bit agreement.
# This makes the diff tests a clean correctness signal: any nonzero
# difference indicates a real bug in the translation, not a numerical
# artifact of compiler optimization.
#
# The parent fmm2d static library (built at -O3 by the top-level
# makefile) is also linked, so any cross-file dependencies that the
# routine under test happens to call are resolved against the
# trusted Fortran versions of those dependencies.

CC      = gcc
FC      = gfortran

# Production flags (used to build the C translation as it would ship).
CFLAGS_PROD = -std=c99 -O3 -march=native -funroll-loops -fPIC \
              -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable

# Diff-test flags (used for the bit-identical correctness check).
CFLAGS_TEST = -std=c99 -O0 -g -fPIC \
              -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable

CFLAGS  = $(CFLAGS_TEST)
FFLAGS  = -O3 -march=native -funroll-loops -std=legacy -fPIC -w
FFLAGS_REF = -O0 -g -std=legacy -fPIC -w

ROOT    = ..
INCDIR  = include
SRCDIR  = src
TESTDIR = test
BUILD   = build

CPPFLAGS = -I$(INCDIR)

# Static lib built by the top-level fmm2d makefile.
FMMLIB  = $(ROOT)/lib-static/libfmm2d.a

# Test helpers (hkrand, dlaran) live as standalone .o files outside the
# library. The diff test drivers use hkrand for randomized inputs.
TESTOBJS = $(ROOT)/src/common/hkrand.o $(ROOT)/src/common/dlaran.o

# Map from translated-file basename to the original Fortran source path.
# Each new translation must be added here so the Makefile can build the
# matching -O0 reference object for the diff test.
FORT_SRC_hndiv2d         = $(ROOT)/src/helmholtz/hndiv2d.f
FORT_SRC_cumsum          = $(ROOT)/src/common/cumsum.f
FORT_SRC_fmmcommon2d     = $(ROOT)/src/common/fmmcommon2d.f
FORT_SRC_l2dterms        = $(ROOT)/src/laplace/l2dterms.f
FORT_SRC_cauchykernels2d = $(ROOT)/src/laplace/cauchykernels2d.f
FORT_SRC_laprouts2d      = $(ROOT)/src/laplace/laprouts2d.f
FORT_SRC_tree_routs2d    = $(ROOT)/src/common/tree_routs2d.f
FORT_SRC_pts_tree2d      = $(ROOT)/src/common/pts_tree2d.f
FORT_SRC_cfmm2d          = $(ROOT)/src/laplace/cfmm2d.f
FORT_SRC_cfmm2d_ndiv     = $(ROOT)/src/laplace/cfmm2d_ndiv.f
FORT_SRC_rfmm2d_ndiv     = $(ROOT)/src/laplace/rfmm2d_ndiv.f

# Each translated file is listed here as it gets added.
C_SRCS  = $(wildcard $(SRCDIR)/*.c)
C_OBJS  = $(patsubst $(SRCDIR)/%.c,$(BUILD)/%.o,$(C_SRCS))

# Each diff test driver is one Fortran main per translated file,
# named test_<file>.f. Build them all as standalone executables.
TEST_SRCS = $(wildcard $(TESTDIR)/test_*.f)
TEST_BINS = $(patsubst $(TESTDIR)/test_%.f,$(BUILD)/test_%,$(TEST_SRCS))

.PHONY: all tests clean dropin e2e

all: $(C_OBJS) $(TEST_BINS)

# Drop-in build: compile each C source at -O3 with -DFMM2D_DROP_IN so
# the resulting objects export the bare Fortran symbol names (e.g.
# `cfmm2d_` instead of `cfmm2d_c_`). These objects can then be linked
# in place of the corresponding objects from libfmm2d.a as a true
# Fortran-ABI drop-in replacement.
DROPIN_DIR = $(BUILD)/dropin
DROPIN_OBJS = $(patsubst $(SRCDIR)/%.c,$(DROPIN_DIR)/%.o,$(C_SRCS))
DROPIN_CFLAGS = -std=c99 -O3 -march=native -funroll-loops -fPIC \
                -DFMM2D_DROP_IN \
                -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable

$(DROPIN_DIR):
	mkdir -p $(DROPIN_DIR)

$(DROPIN_DIR)/%.o: $(SRCDIR)/%.c | $(DROPIN_DIR)
	$(CC) $(DROPIN_CFLAGS) $(CPPFLAGS) -c $< -o $@

dropin: $(DROPIN_OBJS)

# End-to-end smoke test: build the existing Fortran test_rfmm2d.f
# (which calls rfmm2d -> cfmm2d -> ...) against the C drop-in
# library, plus a custom test driver that exercises the matlab path
# rfmm2d_ndiv directly. Both should run cleanly and produce errors at
# the FMM precision target.
$(BUILD)/end2end_rfmm2d: $(DROPIN_OBJS) $(FMMLIB) | $(BUILD)
	$(FC) -O3 -march=native -funroll-loops -std=legacy -w \
	    $(ROOT)/test/laplace/test_rfmm2d.f \
	    $(DROPIN_OBJS) $(TESTOBJS) $(FMMLIB) -lm \
	    -o $@

$(BUILD)/end2end_rfmm2d_ndiv: $(DROPIN_OBJS) $(FMMLIB) e2e/test_rfmm2d_e2e.f | $(BUILD)
	$(FC) -O3 -march=native -funroll-loops -std=legacy -w \
	    e2e/test_rfmm2d_e2e.f \
	    $(DROPIN_OBJS) $(TESTOBJS) $(FMMLIB) -lm \
	    -o $@

e2e: $(BUILD)/end2end_rfmm2d $(BUILD)/end2end_rfmm2d_ndiv
	@echo "==> $(BUILD)/end2end_rfmm2d (existing test_rfmm2d.f via C drop-in)"
	@$(BUILD)/end2end_rfmm2d > $(BUILD)/end2end_rfmm2d.log 2>&1 \
	    && echo "    PASS (output in $(BUILD)/end2end_rfmm2d.log)" \
	    || (echo "    FAIL"; cat $(BUILD)/end2end_rfmm2d.log; exit 1)
	@echo "==> $(BUILD)/end2end_rfmm2d_ndiv (matlab path smoke test)"
	@$(BUILD)/end2end_rfmm2d_ndiv

# Compile a C source at the diff-test optimization level.
$(BUILD)/%.o: $(SRCDIR)/%.c | $(BUILD)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

# Compile the Fortran original of a translated source at -O0 so it
# matches the C translation bit-for-bit. One explicit rule per file —
# the source paths vary across the source tree so a single pattern
# rule isn't possible without secondary expansion gymnastics.
$(BUILD)/hndiv2d_ref.o: $(FORT_SRC_hndiv2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/cumsum_ref.o: $(FORT_SRC_cumsum) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/fmmcommon2d_ref.o: $(FORT_SRC_fmmcommon2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/l2dterms_ref.o: $(FORT_SRC_l2dterms) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/cauchykernels2d_ref.o: $(FORT_SRC_cauchykernels2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/laprouts2d_ref.o: $(FORT_SRC_laprouts2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/tree_routs2d_ref.o: $(FORT_SRC_tree_routs2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/pts_tree2d_ref.o: $(FORT_SRC_pts_tree2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/cfmm2d_ref.o: $(FORT_SRC_cfmm2d) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/cfmm2d_ndiv_ref.o: $(FORT_SRC_cfmm2d_ndiv) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@
$(BUILD)/rfmm2d_ndiv_ref.o: $(FORT_SRC_rfmm2d_ndiv) | $(BUILD)
	$(FC) $(FFLAGS_REF) -c $< -o $@

.SECONDARY:

# Build a diff test. The link order matters: the -O0 Fortran reference
# (`%_ref.o`) is named on the command line BEFORE libfmm2d.a, so the
# linker uses our deterministic copy of the routine under test instead
# of the -O3 copy that lives in the .a. The C translation is in
# `%.o` and exports its symbols with a `_c_` suffix so there is no
# collision.
#
# Static pattern rule (over the explicit list of test binaries) to
# avoid GNU Make's "implicit rule recursion" rejection that happens
# when chaining ordinary pattern rules through %.o intermediates.
$(TEST_BINS): $(BUILD)/test_%: $(TESTDIR)/test_%.f $(BUILD)/%.o $(BUILD)/%_ref.o $(FMMLIB) | $(BUILD)
	$(FC) $(FFLAGS_REF) $< $(BUILD)/$*.o $(BUILD)/$*_ref.o \
	    $(TESTOBJS) $(FMMLIB) -lm -o $@

$(BUILD):
	mkdir -p $(BUILD)

tests: $(TEST_BINS)
	@set -e; for t in $(TEST_BINS); do \
	  echo "==> $$t"; \
	  $$t; \
	done

clean:
	rm -rf $(BUILD)
