# To be installed:
noinst_LIBRARIES = libflib.a
module_DATA      = $(modules)

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
libflib_a_SOURCES = utils.c utilsadd.c f_plugin.c \
	yaml_strings.f90 dictionaries_base.f90 dictionaries.f90 \
	callbacks.f90 f_malloc.f90 \
	f_enums.f90 f_precisions.f90 \
	f_utils.f90 f_refcnts.f90 f_trees.f90 \
	yaml_output.f90 mem_profiling.f90 dynamic_memory.f90 \
	razero.f90 getadd.f90 time.f90 \
	yaml_parser.c yaml_parse.f90 f_jmp.f90 \
	f_input_file.f90 randomData.f90 \
	f_ternary.f90 f_bind.f90 f_objects.f90 fpython.f90 \
	$(flush_files) $(getarg_files) $(CUDA_PROF) $(CUDA_SOURCES)

DYN_INCLUDES =  deallocate-inc.f90 allocate-inc.f90 \
	deallocate-c-inc.f90 allocate-c-inc.f90 \
	deallocate-profile-inc.f90 allocate-profile-inc.f90 \
	malloc_templates-inc.f90 \
	deallocate-multiple-inc.f90 deallocate-multiple-inc-ptr.f90 \
	f_memcpy-inc.f90 f_memcpy-base-inc.f90 f_maxdiff-base-inc.f90 \
	allocate-ptr-inc.f90 f_map-inc.f90

F_MALLOC_INCLUDES = f_malloc-inc.f90 f_malloc-extra-inc.f90 \
	f_malloc-base-inc.f90 f_malloc-null-inc.f90 \
	f_malloc-bound-inc.f90 \
	f_malloc-bounds-inc.f90 f_malloci-simple-inc.f90 f_mallocli-simple-inc.f90 \
	f_malloc-ptr-inc.f90 f_malloc-check-inc.f90 \
	f_malloc-null-ptr-inc.f90 f_malloc-total-inc.f90 \
	f_malloc-buf-inc.f90 f_malloc-buf-base-inc.f90

YAML_INCLUDES = yaml_map-inc.f90 yaml_map-arr-inc.f90 \
	yaml_map-mat-inc.f90 yaml_open-inc.f90
STRINGS_INCLUDES = yaml_toa-inc.f90 yaml_toa-arr-inc.f90 yaml_plus-inc.f90
ERR_INCLUDES = get_err-inc.f90
DICT_INCLUDES = error_handling.f90 $(ERR_INCLUDES) dict_add-inc.f90 dict_getvec-inc.f90 \
	set_arr-inc.f90 dict_cont-inc.f90 dict_cont_arr-inc.f90
GETADD_INCLUDES = getadd-c-inc.f90
NOOMP_INCLUDES = remove_omp-inc.f90 halt_omp-inc.f90

EXTRA_DIST = $(DYN_INCLUDES) $(YAML_INCLUDES) $(STRINGS_INCLUDES) \
	$(DICT_INCLUDES) $(GETADD_INCLUDES) $(F_MALLOC_INCLUDES) $(NOOMP_INCLUDES)

# Module files.
if CAPITALIZE
modules = YAML_STRINGS.@MODULE_EXT@ \
	DICTIONARIES.@MODULE_EXT@ \
	DICTIONARIES_BASE.@MODULE_EXT@ \
	YAML_OUTPUT.@MODULE_EXT@ \
	F_ENUMS.@MODULE_EXT@ \
	F_PRECISIONS.@MODULE_EXT@ \
	YAML_PARSE.@MODULE_EXT@ \
	METADATA_INTERFACES.@MODULE_EXT@ \
	MODULE_F_MALLOC.@MODULE_EXT@ \
	DYNAMIC_MEMORY.@MODULE_EXT@ \
	DYNAMIC_MEMORY_BASE.@MODULE_EXT@ \
	RANDOMDATA.@MODULE_EXT@ \
	RANDOM.@MODULE_EXT@ \
	MEMORY_PROFILING.@MODULE_EXT@ \
	F_UTILS.@MODULE_EXT@ \
	F_REFCNTS.@MODULE_EXT@ \
	F_TREES.@MODULE_EXT@ \
	F_INPUT_FILE.@MODULE_EXT@ \
	EXCEPTION_CALLBACKS.@MODULE_EXT@ \
	F_JMP.@MODULE_EXT \
	TIME_PROFILING.@MODULE_EXT@ \
	F_TERNARY.@MODULE_EXT@ \
	MODULE_F_BIND.@MODULE_EXT@ \
	MODULE_F_OBJECTS.@MODULE_EXT@ \
	F_PYTHON.@MODULE_EXT@ \
	MODULE_RAZERO.@MODULE_EXT@ $(CUDA_MODS)
else
modules = yaml_strings.@MODULE_EXT@ \
	dictionaries.@MODULE_EXT@ \
	dictionaries_base.@MODULE_EXT@ \
	f_enums.@MODULE_EXT@ \
	f_precisions.@MODULE_EXT@ \
	yaml_output.@MODULE_EXT@ \
	yaml_parse.@MODULE_EXT@ \
	metadata_interfaces.@MODULE_EXT@ \
	module_f_malloc.@MODULE_EXT@ \
	dynamic_memory.@MODULE_EXT@ \
	dynamic_memory_base.@MODULE_EXT@ \
	randomdata.@MODULE_EXT@ \
	random.@MODULE_EXT@ \
	memory_profiling.@MODULE_EXT@ \
	f_utils.@MODULE_EXT@ \
	f_input_file.@MODULE_EXT@ \
	f_refcnts.@MODULE_EXT@ \
	f_trees.@MODULE_EXT@ \
	f_ternary.@MODULE_EXT@ \
	module_f_bind.@MODULE_EXT@ \
	module_f_objects.@MODULE_EXT@ \
	f_python.@MODULE_EXT@ \
	exception_callbacks.@MODULE_EXT@ \
	time_profiling.@MODULE_EXT@ \
	f_jmp.@MODULE_EXT@ \
	module_razero.@MODULE_EXT@ $(CUDA_MODS)
endif

CLEANFILES = *.@MODULE_EXT@ memstatus.yaml

SUFFIXES = .cu

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@
AM_CFLAGS = @LIB_YAML_CFLAGS@
AM_CPPFLAGS = -I. -I$(top_builddir) -I@CUDA_PATH@/include/ 

NVCC = @NVCC@
NVCC_FLAGS = @NVCC_FLAGS@
.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@


# Dependencies
$(modules): $(libflib_a_OBJECTS)
yaml_strings.o: $(STRINGS_INCLUDES) f_precisions.o
f_enums.o: yaml_strings.o
f_bind.o: f_precisions.o
f_objects.o: f_precisions.o dictionaries.o yaml_output.o
f_ternary.o: f_precisions.o yaml_strings.o
dictionaries_base.o: yaml_strings.o
dictionaries.o: yaml_strings.o dictionaries_base.o $(DICT_INCLUDES) callbacks.o
dictionaries_base.o: yaml_strings.o
callbacks.o: $(NOOMP_INCLUDES)
$(NVTX_DEP): $(NOOMP_INCLUDES)
#mem_profiling.o : error_handling.o
yaml_output.o : dictionaries.o f_utils.o f_enums.o $(YAML_INCLUDES)
yaml_parse.o : dictionaries.o yaml_output.o
#error_handling.o: dictionaries.o callbacks.o yaml_output.o $(ERR_INCLUDES)
f_malloc.o: dictionaries.o $(F_MALLOC_INCLUDES)
f_input_file.o: dictionaries.o yaml_output.o dynamic_memory.o f_utils.o
dynamic_memory.o: getadd.o yaml_strings.o yaml_output.o yaml_parse.o mem_profiling.o f_malloc.o f_utils.o $(DYN_INCLUDES)
mem_profiling.o: yaml_output.o yaml_parse.o yaml_strings.o
f_lib_highlev.o: dictionaries.o yaml_output.o dynamic_memory.o yaml_parse.o $(NVTX_DEP) time.o f_utils.o f_input_file.o
f_utils.o: dictionaries.o yaml_strings.o f_precisions.o razero.o
f_refcnts.o: dictionaries.o yaml_output.o f_precisions.o
f_trees.o: f_refcnts.o dictionaries.o yaml_output.o
getadd.o: yaml_strings.o $(GETADD_INCLUDES) f_malloc.o
razero.o: randomData.o f_precisions.o
$(flush_object): yaml_output.o dictionaries.o
fpython.o: dictionaries.o yaml_output.o f_objects.o

time.o: $(NOOMP_INCLUDES) dynamic_memory.o f_utils.o $(NVTX_DEP) f_jmp.o 

f_jmp.o: f_precisions.o dynamic_memory.o
