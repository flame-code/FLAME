
# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_ACC_CUDA],
[dnl Use the GPU computing (only NVIDIA cards with CUDA libraries)
  ax_have_acc_cuda="no"
  AC_ARG_ENABLE(cuda-gpu, AS_HELP_STRING([--enable-cuda-gpu], [Use CUDA implementation  (disabled by default).]),
                           ax_have_acc_cuda=$enableval, ax_have_acc_cuda="no")

  if test x"$ax_have_acc_cuda" = "xyes" ; then
    dnl Optional Cuda installation dir.
    AC_ARG_WITH(cuda-path, AS_HELP_STRING([--with-cuda-path], [give the path to the NVidia Cuda tools (default is /usr/local/cuda).]),
                           [ax_acc_cuda_path=$withval], [ax_acc_cuda_path="/usr/local/cuda"])
    AC_ARG_WITH(cuda-libs, AS_HELP_STRING([--with-cuda-libs], [give the link flags for nvidia libraries]),
                           [ax_acc_cuda_libs=$withval], [ax_acc_cuda_libs="-lcudart -lcublas -lcufft"])
    if test x"$ax_acc_cuda_path" = x"no" ; then
      ax_acc_cuda_path="/usr"
    fi
    CUDA_PATH=$ax_acc_cuda_path
    AC_SUBST(CUDA_PATH)
    AC_MSG_CHECKING([looking for CUDA in])
    AC_MSG_RESULT([$ax_acc_cuda_path])
    AC_MSG_CHECKING([Linking flags for CUDA])
    AC_MSG_RESULT([$ax_acc_cuda_libs])

    dnl look for the cublas header
    AC_LANG_PUSH(C)
    CPPFLAGS_SVG=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS -I$ax_acc_cuda_path/include"
    AC_CHECK_HEADER([cublas.h],
                    [ax_have_cublas="yes"],
                    [ax_have_cublas="no"])
    CPPFLAGS=$CPPFLAGS_SVG
    if test "$ax_have_cublas" = "no"; then
      AC_MSG_WARN(["No 'CUDA' header found, CUDA parts disabled."])
      ax_have_acc_cuda="no"
    else
      CUDA_INCLUDE_PATH=$ax_acc_cuda_path"/include"
      AC_SUBST(CUDA_INCLUDE_PATH)
      CUDA_INCLUDE="-I"$ax_acc_cuda_path"/include"
      AC_SUBST(CUDA_INCLUDE)
    fi

    dnl Test the link on libcudart and libcublas.
    LIBS_OLD=$LIBS

    dnl look for the cublas library
    ac_cuda_lib_path="$ax_acc_cuda_path/lib"
    LIBS="-L$ac_cuda_lib_path $LIBS_OLD"
    AC_CHECK_LIB(cublas, cublasSsymv, withlibcublas=yes, withlibcublas=no, [-lcudart])
    if test "$withlibcublas" = "no"; then
      AC_CHECK_LIB(cublas, cublasSsymv_v2, withlibcublas=yes, withlibcublas=no, [-lcudart])
    fi
    if test "$withlibcublas" = "no"; then
      ac_cuda_lib_path="$ax_acc_cuda_path/lib64"
      LIBS="-L$ac_cuda_lib_path $LIBS_OLD"
      AC_CHECK_LIB(cublas, cublasScopy, withlibcublas=yes, withlibcublas=no, [-lcudart])
      if test "$withlibcublas" = "no"; then
        AC_CHECK_LIB(cublas, cublasScopy_v2, withlibcublas=yes, withlibcublas=no, [-lcudart])
      fi
    fi
    if test "$withlibcublas" = "no"; then
      AC_MSG_ERROR(["No 'cublas' library found, link will fail."])
    fi
    LIBS=$LIBS_OLD
    AC_LANG_POP(C)
  
    dnl Test the existence of things.
    AC_PATH_PROG([NVCC], [nvcc], [$ac_with_nvcc], [$PATH:${ax_acc_cuda_path}/bin])
    AC_SUBST(NVCC)
  
    dnl Add the flags.
    AC_MSG_CHECKING([for NVCC flags])
    if test -z "$NVCC_FLAGS" ; then
      AC_REQUIRE([AX_FLAG_PIC])
      NVCC_FLAGS="-arch sm_20 -O3 --compiler-options '-fno-strict-aliasing $ax_flag_pic'"
    fi
    AC_SUBST(NVCC_FLAGS)
    AC_MSG_RESULT([$NVCC_FLAGS])
    LIBCUDA_LIBS="$ax_acc_cuda_libs -lstdc++"
    LDFLAGS="$LDFLAGS -L$ac_cuda_lib_path"
  fi

  AC_MSG_CHECKING([for GPU computing with CUDA libraries])
  AC_MSG_RESULT([$ax_have_acc_cuda])

  AM_CONDITIONAL(USE_CUDA_GPU, test "$ax_have_acc_cuda" = "yes")
])

AC_DEFUN([AX_ACC_OCL],
[dnl Use the OpenCL support
  ax_have_acc_ocl="no"
  AC_ARG_ENABLE(opencl, AS_HELP_STRING([--enable-opencl], [Use OpenCL implementation for GPU convolutions (disabled by default).]),
                           ax_have_acc_ocl=$enableval, ax_have_acc_ocl="no")
  AC_MSG_CHECKING([for OpenCL support])
  AC_MSG_RESULT([$ax_have_acc_ocl])
  if test x"$ax_have_acc_ocl" = "xyes" ; then
    AC_LANG_PUSH(C)
    dnl Optional OpenCL installation dir.
    AC_ARG_WITH(ocl-path, AS_HELP_STRING([--with-ocl-path], [give the path to the OpenCL library (default is /usr).]),
                           [ac_ocl_path=$withval], [ac_ocl_path=""])
    dnl Test the link on libOpenCL.
    dnl look for the header.
    CPPFLAGS_SVG=$CPPFLAGS
    if test -n "$ac_ocl_path" ; then
       CPPFLAGS="$CPPFLAGS -I$ac_ocl_path/include"
    fi
    AC_CHECK_HEADER([CL/cl.h],
                    [ac_ocl_header="yes"],
                    [ac_ocl_header="no"])
    CPPFLAGS=$CPPFLAGS_SVG
    if test "$ac_ocl_header" = "no"; then
      AC_MSG_WARN(["No 'OpenCL' header found, OpenCL parts disabled."])
      ax_have_acc_ocl="no"
    else
      if test -n "$ac_ocl_path" ; then
         OCL_INCLUDE_PATH=$ac_ocl_path"/include"
         OCL_INCLUDE="-I"$ac_ocl_path"/include"
      fi
      AC_SUBST(OCL_INCLUDE_PATH)
      AC_SUBST(OCL_INCLUDE)
      dnl Now test for the library.
      LIBS_OLD=$LIBS
      if test -n "$ac_ocl_path" ; then
         LIBS="-L$ac_ocl_path/lib $LIBS"
      fi
      AC_CHECK_LIB(OpenCL, clCreateContext, withlibocl=yes, withlibocl=no)
      if test "$withlibocl" = "no"; then
        if test -n "$ac_ocl_path" ; then
          LIBS="-L$ac_ocl_path/lib64 $LIBS"
        fi
        AC_CHECK_LIB(OpenCL, clCreateContext, withlibocl=yes, withlibocl=no)
      fi
      if test "$withlibocl" = "no"; then
        AC_MSG_WARN(["No 'OpenCL' library found, OpenCL parts disabled."])
        ax_have_acc_ocl="no"
      else
        if test -n "$ac_ocl_path" ; then
           LDFLAGS="$LDFLAGS -L$ac_ocl_path/lib"
        fi
        LIBOCL_LIBS="-lOpenCL -lm -lstdc++"
      fi
      LIBS=$LIBS_OLD
    fi
    AC_LANG_POP(C)
  fi
  AM_CONDITIONAL(USE_OCL, test "$ax_have_acc_ocl" = "yes")
])

AC_DEFUN([AX_ACC_MIC],
[dnl Use the OpenCL support
  ax_have_acc_mic="no"
  AC_ARG_ENABLE(intel-mic, AS_HELP_STRING([--enable-intel-mic], [Use MIC implementation for GPU convolutions (disabled by default).]),
                           ax_have_acc_mic=$enableval, ax_have_acc_mic="no")
  AC_MSG_CHECKING([for Intel MIC support])
  AC_MSG_RESULT([$ax_have_acc_mic])
  AC_ARG_WITH(intel-mic-libs, AS_HELP_STRING([--with-intel-mic-libs], [Link MIC executable with the following addional libraries (default = empty).]), ac_mic_libs=$withval, ac_mic_libs=no)
  if test x"$ax_have_acc_mic" = "xyes" ; then
    if test x"$ac_mic_libs" != "xno" ; then
      LDFLAGS="-offload-option,mic,ld,\"$ac_mic_libs\" $LDFLAGS"
    fi
  fi
  AM_CONDITIONAL(USE_MIC, test "$ax_have_acc_mic" = "yes")
])
