#!/bin/bash

rm -rf futile-tmp
mkdir futile-tmp
mkdir futile-tmp/futile
mkdir futile-tmp/yaml-0.1.6
mkdir futile-tmp/futile/dicts
mkdir futile-tmp/futile/c-bindings
mkdir futile-tmp/futile/src
mkdir futile-tmp/futile/flib
mkdir futile-tmp/futile/wrappers
mkdir futile-tmp/futile/wrappers/fft
mkdir futile-tmp/futile/wrappers/dgemmsy
mkdir futile-tmp/futile/wrappers/mpi

dicts=(err.h err.c dictionaries_base.f90 callbacks.f90 f_precisions.f90 dictionaries.f90 \
    addresses.c dictf.f90 dict.c yaml_strings.f90 err-fapi.h dict.h errf.f90 dict-fapi.h \
    addresses.h)
dicts_h=(yaml_toa-inc.f90 get_err-inc.f90 yaml_toa-arr-inc.f90 dict_cont-inc.f90 dict_cont_arr-inc.f90 \
    set_arr-inc.f90 dict_add-inc.f90 dict_getvec-inc.f90 yaml_plus-inc.f90 halt_omp-inc.f90 error_handling.f90)

cbindings=(misc.c treef.f90 python.c yaml_output_c.c futile.h tree.h yaml_output_bind.f90 f_progress_bar.h \
    futileC_objects.h enum.h futileC_objects.c tree.c yaml_output.h f_string.h misc.h)

src=(f_lib_highlev.f90)

flib_h=(allocate-simgrid-inc.f90 f_malloc-bounds-inc.f90 deallocate-multiple-ptr-inc.f90 \
    allocate-base-inc.f90 deallocate-c-inc.f90 f_malloc-check-inc.f90 f_malloc-bound-inc.f90 \
    f_map-inc.f90 f_malloc-null-ptr-inc.f90 remove_omp-inc.f90 getadd-c-inc.f90 \
    deallocate-multiple-inc.f90 f_malloc-buf-base-inc.f90 allocate-end-inc.f90 yaml_open-inc.f90 \
    f_maxdiff-base-inc.f90 deallocate-inc.f90 allocate-inc.f90 allocate-aligned-inc.f90 \
    f_malloc-ptr-inc.f90 f_malloc-null-inc.f90 yaml_map-inc.f90 allocation-validation-inc.f90 \
    f_mallocli-simple-inc.f90 f_malloc-extra-inc.f90 malloc_templates-c-inc.f90 \
    allocate-profile-inc.f90 f_memcpy-base-inc.f90 allocate-c-inc.f90 deallocate-profile-inc.f90 \
    yaml_map-arr-inc.f90 f_memcpy-inc.f90 malloc_templates-inc.f90 f_malloc-buf-inc.f90 \
    f_malloci-simple-inc.f90 f_malloc-total-inc.f90 f_malloc-base-inc.f90 compare-scalar-inc.f90 \
    deallocate-ptr-inc.f90 f_malloc-inc.f90 yaml_map-mat-inc.f90 compare-array-inc.f90 \
    allocate-ptr-inc.f90)

flib=(simgrid_shared_fake.f90 f_input_file.f90 yaml_parse.f90 f_refcnts.f90 \
    performances.f90 f_malloc.f90 f_objects.f90 f_plugin.c time.f90 razero.f90 f_iostream.f90 \
    yaml_parser.c f_ternary.f90 getadd.f90 yaml_output.f90 fpython.f90 dynamic_memory.f90 \
    flush.f90 f_enums.f90 malloc_repad.c utilsadd.c cudall_fake.f90 f_arrays.f90 \
    randomData.f90 nvtx_fake.f90 f_regtests.f90 utils.h f_jmp.f90 utils.c mem_profiling.f90 \
    f_environment.f90 f_trees.f90 override_malloc.c get_command_argument_fake.f90 f_utils.f90 \
    f_bibliography.f90)

wrappers=(interface_ocl_fake.f90 op2p_module.f90 cublas_fake.f90 linalg.f90 \
    interface_dgemmsy_fake.f90 f_blas.f90)

wrappers_fft=(fft_factors.f90 fft3d.f90 fft2d.f90)

wrappers_dgemmsy=(dgemmsy.h visitors.h gemm_block_c.c patterns.h visitors.c \
    utils.h patterns.c gemm_block_c.h)

wrappers_mpi=(mpi.F90 MPI3fake.c scalars.f90 gather.f90 allreduce.f90 allgather.f90 \
    test_mpi_wrappers.f90 bcast.f90 onesided.f90 sendrecv.f90 alltoall.f90)

wrappers_mpi_h=(mpi_get_alltoallv-inc.f90 maxdiff-inc.f90 bcast-inc.f90 gather-inc.f90 \
    allreduce-decl-inc.f90 scatterv-decl-inc.f90 maxdiff-decl-inc.f90 get-inc.f90 \
    maxdiff-arr-inc.f90 maxdiff-routines-inc.f90 bcast-decl-inc.f90 scatterv-inc.f90 \
    allgather-inc.f90 allreduce-core-inc.f90 win-create-inc.f90 get-decl-inc.f90 \
    alltoallv-inc.f90 get-end-inc.f90 allreduce-inc.f90 bcast-decl-arr-inc.f90 \
    scatter-inc.f90 control-routines-inc.f90 gather-inner-inc.f90 allreduce-arr-inc.f90 \
    allreduce-multi-inc.f90 maxdiff-end-inc.f90)

#------ Directory dicts --------------------------------------------------------
for file in ${dicts[@]}; do
  cp -v $FUTILE_DIR/futile/dicts/$file futile-tmp/futile/dicts
done
for file in ${dicts_h[@]}; do
  cp -v $FUTILE_DIR/futile/dicts/$file futile-tmp/futile/dicts/"$file".h
done
cp -v $FUTILE_DIR/futile/dicts/futile_cst.h.in futile-tmp/futile/dicts/futile_cst.h
sed -i 's/@GLIB_TRUE@/\/\*/g'     futile-tmp/futile/dicts/futile_cst.h
sed -i 's/@GLIB_END_TRUE@/\*\//g' futile-tmp/futile/dicts/futile_cst.h
sed -i 's/@GLIB_FALSE@//g'        futile-tmp/futile/dicts/futile_cst.h
sed -i 's/@GLIB_END_FALSE@//g'    futile-tmp/futile/dicts/futile_cst.h
#------ Directory c-bindings ---------------------------------------------------
for file in ${cbindings[@]}; do
  cp -v $FUTILE_DIR/futile/c-bindings/$file futile-tmp/futile/c-bindings/$file
done
#------ Directory src ----------------------------------------------------------
for file in ${src[@]}; do
  cp -v $FUTILE_DIR/futile/src/$file futile-tmp/futile/src/$file
done
cp -v $FUTILE_DIR/futile/src/configure.inc.in futile-tmp/futile/src/configure.inc
sed -i 's/@PACKAGE_NAME@/Futile/g'                     futile-tmp/futile/src/configure.inc
sed -i 's/@PACKAGE_TARNAME@/futile/g'                  futile-tmp/futile/src/configure.inc
sed -i 's/@PACKAGE_VERSION@/1.8/g'                     futile-tmp/futile/src/configure.inc
sed -i 's/@PACKAGE_STRING@/Futile 1.8/g'               futile-tmp/futile/src/configure.inc
sed -i 's/@PACKAGE_BUGREPORT@/Damien.Caliste@cea.fr/g' futile-tmp/futile/src/configure.inc
#------ Directory flib ---------------------------------------------------------
for file in ${flib[@]}; do
  cp -v $FUTILE_DIR/futile/flib/$file futile-tmp/futile/flib/$file
done
for file in ${flib_h[@]}; do
  cp -v $FUTILE_DIR/futile/flib/$file futile-tmp/futile/flib/"$file".h
done
cp -v $FUTILE_DIR/futile/flib/f_utils.inc.in futile-tmp/futile/flib/f_utils.inc
sed -i 's/@RECL_INT_KIND@/16/g'                     futile-tmp/futile/flib/f_utils.inc
#------ Directory wrappers ---------------------------------------------------------
for file in ${wrappers[@]}; do
  cp -v $FUTILE_DIR/futile/wrappers/$file futile-tmp/futile/wrappers/$file
done
#------ Directory wrappers_fft ---------------------------------------------------------
for file in ${wrappers_fft[@]}; do
  cp -v $FUTILE_DIR/futile/wrappers/fft/$file futile-tmp/futile/wrappers/fft/$file
done
#------ Directory wrappers_dgemmsy ---------------------------------------------------------
for file in ${wrappers_dgemmsy[@]}; do
  cp -v $FUTILE_DIR/futile/wrappers/dgemmsy/$file futile-tmp/futile/wrappers/dgemmsy/$file
done
#------ Directory wrappers_mpi ---------------------------------------------------------
for file in ${wrappers_mpi[@]}; do
  cp -v $FUTILE_DIR/futile/wrappers/mpi/$file futile-tmp/futile/wrappers/mpi/$file
done
for file in ${wrappers_mpi_h[@]}; do
  cp -v $FUTILE_DIR/futile/wrappers/mpi/$file futile-tmp/futile/wrappers/mpi/"$file".h
done
#------ Directory ?????????????---------------------------------------------------------

find futile-tmp -type f -iname '*.f90' -exec sed -i "s/\.f90'/\.f90\.h'/g" {} \;
find futile-tmp -type f -iname '*.f90.h' -exec sed -i "s/\.f90'/\.f90\.h'/g" {} \;

cp files_for_fpm/config.h futile-tmp/

cd futile-tmp
ln -s futile/dicts/futile_cst.h futile_cst.h
ln -s futile/dicts/dict-fapi.h dict-fapi.h
ln -s futile/c-bindings/misc.h misc.h
ln -s futile/dicts/dict.h dict.h
ln -s futile/dicts/addresses.h addresses.h
cd -

sed -i '11i use f_bcast, only: fmpi_maxdiff' futile-tmp/futile/wrappers/mpi/allgather.f90

