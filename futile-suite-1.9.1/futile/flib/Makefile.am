# To be installed:
noinst_LIBRARIES = libflib.a libfmalloc.a
module_DATA      = $(modules)

#malloc override tools
dynlibdir        = $(libdir)
dynlib_DATA      = $(fmalloc_dynlib)
noinst_DATA      = $(fmalloc_linklib)

libfmalloc_a_SOURCES = override_malloc.c

if CC_AUTOBUILD_DYNAMIC_LIBS
fmalloc_dynlib = $(fmalloc_dynamic_library)
fmalloc_linklib = $(fmalloc_library_soname)
else
fmalloc_dynlib = 
fmalloc_linklib = 
endif

fmalloc_dynamic_library = libfmalloc-1.so.@FUTILE_MINOR_VERSION@.0.@FUTILE_MICRO_VERSION@
fmalloc_library_soname = libfmalloc-1.so.@FUTILE_MINOR_VERSION@

libfmalloc.a: $(libfmalloc_a_OBJECTS)
	$(AR) $(ARFLAGS) $@ $(libfmalloc_a_OBJECTS)
	$(RANLIB) $@

override_malloc.o: override_malloc.c
	$(CC) -I. -I$(top_builddir) $(CFLAGS) @CC_PIC_FLAG@ -c $< -o $@

$(fmalloc_dynamic_library): $(libfmalloc_a_OBJECTS)
	$(CC) $(CFLAGS) @CC_BUILD_SHARED@ $(LDFLAGS) @CC_RPATH@$(DESTDIR)$(dynlibdir) -Wl,-soname=$(fmalloc_library_soname) -o $@ $(libfmalloc_a_OBJECTS) -ldl
	@chmod a+x $@

$(fmalloc_library_soname): $(fmalloc_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libfmalloc-1.so

install-data-hook: 
	if test -n "$(fmalloc_dynlib)" ; then \
	cd $(DESTDIR)$(dynlibdir); \
	chmod a+x $(fmalloc_dynamic_library); \
	ln -fs $(fmalloc_dynamic_library) $(fmalloc_library_soname); \
	ln -fs $(fmalloc_library_soname) libfmalloc-1.so ;\
	fi

if USE_CUDA_GPU
CUDA_PROF = nvtx_fake.f90 #nvtx.f90
CUDA_MODS = 
NVTX_DEP = nvtx_fake.o #nvtx.o
CUDA_SOURCES = cudall.cu
else
CUDA_PROF = nvtx_fake.f90
NVTX_DEP = nvtx_fake.o
CUDA_SOURCES = cudall_fake.f90
endif

if USE_SIMGRID_SHARED_ALLOCS
SIMGRID_SHARED_SOURCE = simgrid_shared.f90 #nvtx.f90
SIMGRID_SHARED_DEP = simgrid_shared.o
else
SIMGRID_SHARED_SOURCE = simgrid_shared_fake.f90
SIMGRID_SHARED_DEP = simgrid_shared_fake.o
endif

# Module files.
if CAPITALIZE
CUDA_MODS = NVTX.@MODULE_EXT@
else
CUDA_MODS = nvtx.@MODULE_EXT@
endif

# Source code.
if HAVE_FC_FLUSH_SUB
flush_files = flush_sub.f90
flush_object= flush_sub.o
else
if HAVE_FC_FLUSH
flush_files = flush.f90
flush_object= flush.o
else
flush_files = flush_fake.f90
flush_object= flush_fake.o
endif
endif

if HAVE_FC_GET_COMMAND_ARGUMENT
getarg_files =
else
getarg_files = get_command_argument_fake.f90
endif
libflib_a_SOURCES = utils.c utils.h utilsadd.c f_plugin.c malloc_repad.c \
	f_malloc.f90 f_enums.f90 \
	f_utils.f90 f_refcnts.f90 f_trees.f90 performances.f90 \
	yaml_output.f90 mem_profiling.f90 dynamic_memory.f90 \
	razero.f90 getadd.f90 time.f90 f_bibliography.f90 f_iostream.f90 \
	yaml_parser.c yaml_parse.f90 f_jmp.f90 f_regtests.f90 \
	f_input_file.f90 randomData.f90 \
	f_ternary.f90 f_objects.f90 fpython.f90 f_environment.f90 f_arrays.f90\
	$(flush_files) $(getarg_files) $(CUDA_PROF) $(SIMGRID_SHARED_SOURCE) $(CUDA_SOURCES)

DYN_INCLUDES =  deallocate-inc.f90 allocate-inc.f90 \
	deallocate-c-inc.f90 allocate-c-inc.f90 \
	deallocate-profile-inc.f90 allocate-profile-inc.f90 \
	malloc_templates-inc.f90 \
	deallocate-multiple-inc.f90 deallocate-multiple-ptr-inc.f90 \
	f_memcpy-inc.f90 f_memcpy-base-inc.f90 f_maxdiff-base-inc.f90 \
	allocate-ptr-inc.f90 deallocate-ptr-inc.f90 f_map-inc.f90 \
	allocate-simgrid-inc.f90 allocate-aligned-inc.f90\
	allocate-base-inc.f90 allocate-end-inc.f90 allocation-validation-inc.f90

F_MALLOC_INCLUDES = f_malloc-inc.f90 f_malloc-extra-inc.f90 \
	f_malloc-base-inc.f90 f_malloc-null-inc.f90 \
	f_malloc-bound-inc.f90 \
	f_malloc-bounds-inc.f90 f_malloci-simple-inc.f90 f_mallocli-simple-inc.f90 \
	f_malloc-ptr-inc.f90 f_malloc-check-inc.f90 \
	f_malloc-null-ptr-inc.f90 f_malloc-total-inc.f90 \
	f_malloc-buf-inc.f90 f_malloc-buf-base-inc.f90

YAML_INCLUDES = yaml_map-inc.f90 yaml_map-arr-inc.f90 \
	yaml_map-mat-inc.f90 yaml_open-inc.f90
GETADD_INCLUDES = getadd-c-inc.f90
NOOMP_INCLUDES = remove_omp-inc.f90
REGTESTS_INCLUDES = compare-scalar-inc.f90 compare-array-inc.f90

EXTRA_DIST = $(DYN_INCLUDES) $(YAML_INCLUDES) $(F_MALLOC_INCLUDES) $(NOOMP_INCLUDES) $(GETADD_INCLUDES) $(REGTESTS_INCLUDES)

# Module files.
if CAPITALIZE
modules = YAML_OUTPUT.@MODULE_EXT@ \
	F_ENUMS.@MODULE_EXT@ \
	F_ENVIRONMENT.@MODULE_EXT@ \
	YAML_PARSE.@MODULE_EXT@ \
	METADATA_INTERFACES.@MODULE_EXT@ \
	MODULE_F_MALLOC.@MODULE_EXT@ \
	DYNAMIC_MEMORY.@MODULE_EXT@ \
	DYNAMIC_MEMORY_BASE.@MODULE_EXT@ \
	F_RANDOM.@MODULE_EXT@ \
	MEMORY_PROFILING.@MODULE_EXT@ \
	F_UTILS.@MODULE_EXT@ \
	F_IOSTREAM.@MODULE_EXT@ \
	F_PERFS.@MODULE_EXT@ \
	F_BIBLIOGRAPHY.@MODULE_EXT@ \
	F_REFCNTS.@MODULE_EXT@ \
	F_TREES.@MODULE_EXT@ \
	F_INPUT_FILE.@MODULE_EXT@ \
	F_JMP.@MODULE_EXT@ \
	TIME_PROFILING.@MODULE_EXT@ \
	F_TERNARY.@MODULE_EXT@ \
	MODULE_F_OBJECTS.@MODULE_EXT@ \
	F_PYTHON.@MODULE_EXT@ \
	SMPI_SHARED.@MODULE_EXT@ \
	F_ARRAYS.@MODULE_EXT@ \
	F_REGTESTS.@MODULE_EXT@ \
	MODULE_RAZERO.@MODULE_EXT@ $(CUDA_MODS)
else
modules = f_enums.@MODULE_EXT@ \
	f_environment.@MODULE_EXT@ \
	yaml_output.@MODULE_EXT@ \
	yaml_parse.@MODULE_EXT@ \
	metadata_interfaces.@MODULE_EXT@ \
	module_f_malloc.@MODULE_EXT@ \
	dynamic_memory.@MODULE_EXT@ \
	dynamic_memory_base.@MODULE_EXT@ \
	f_random.@MODULE_EXT@ \
	f_bibliography.@MODULE_EXT@ \
	memory_profiling.@MODULE_EXT@ \
	f_utils.@MODULE_EXT@ \
	f_iostream.@MODULE_EXT@ \
	f_perfs.@MODULE_EXT@ \
	f_input_file.@MODULE_EXT@ \
	f_refcnts.@MODULE_EXT@ \
	f_trees.@MODULE_EXT@ \
	f_ternary.@MODULE_EXT@ \
	module_f_objects.@MODULE_EXT@ \
	f_python.@MODULE_EXT@ \
	time_profiling.@MODULE_EXT@ \
	f_jmp.@MODULE_EXT@ \
	smpi_shared.@MODULE_EXT@ \
	f_arrays.@MODULE_EXT@ \
	f_regtests.@MODULE_EXT@ \
	module_razero.@MODULE_EXT@ $(CUDA_MODS)
endif

CLEANFILES = libfmalloc-1.so* *.@MODULE_EXT@ memstatus.yaml

DISTCLEANFILES = f_utils.inc

SUFFIXES = .cu

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/dicts @MPI_INCLUDE@
AM_CFLAGS = -I$(top_srcdir)/dicts @LIB_YAML_CFLAGS@
AM_CPPFLAGS = -I. -I$(top_builddir) -I@CUDA_PATH@/include/ 

NVCC = @NVCC@
NVCC_FLAGS = @NVCC_FLAGS@
.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@


# Dependencies
$(modules): $(libflib_a_OBJECTS)

yaml_output.o : f_utils.o f_enums.o f_bibliography.o $(YAML_INCLUDES)
yaml_parse.o : yaml_output.o
f_utils.o: razero.o f_malloc.o
f_regtests.o: yaml_output.o

f_objects.o: yaml_output.o
$(NVTX_DEP): $(NOOMP_INCLUDES)
$(SIMGRID_DEP): $(NOOMP_INCLUDES)
f_malloc.o: $(SIMGRID_DEP) $(F_MALLOC_INCLUDES)
f_input_file.o: yaml_output.o dynamic_memory.o f_utils.o
dynamic_memory.o: getadd.o yaml_output.o yaml_parse.o mem_profiling.o f_malloc.o f_utils.o $(DYN_INCLUDES) f_environment.o performances.o\
	$(SIMGRID_SHARED_DEP)
mem_profiling.o: yaml_output.o yaml_parse.o f_environment.o
f_lib_highlev.o: yaml_output.o dynamic_memory.o yaml_parse.o $(NVTX_DEP) $(SIMGRID_DEP) time.o f_utils.o f_input_file.o
f_refcnts.o: yaml_output.o
f_trees.o: f_refcnts.o yaml_output.o yaml_parse.o
getadd.o: $(GETADD_INCLUDES) f_malloc.o
razero.o: randomData.o
f_arrays.o: dynamic_memory.o f_refcnts.o
$(flush_object): yaml_output.o
fpython.o: yaml_output.o f_objects.o

time.o: $(NOOMP_INCLUDES) dynamic_memory.o f_utils.o $(NVTX_DEP) f_jmp.o 

f_jmp.o: dynamic_memory.o
f_iostream.o: f_utils.o
