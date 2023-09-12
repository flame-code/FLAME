#!/bin/bash
  
write_fpm_rsp_common() {
    echo "@$1" >>$FPM_RSP
    echo -n "option $1 --compiler $FC --c-compiler $CC --cxx-compiler $CXX" >>$FPM_RSP
    #--flag
    echo -n " --flag \"$FFLAGS -DMPI" >>$FPM_RSP
    if [ "$decision_BigDFT" == true ]; then
        echo -n " -DHAVE_BPS -Ibigdft_install/include" >>$FPM_RSP
    fi
    if [ "$decision_SPGLIB" == true ]; then
        #echo -n " -DSPGLIB -Ispglib_install" >>$FPM_RSP
        echo -n " -DSPGLIB" >>$FPM_RSP
    fi
    if [ "$decision_LAMMPS" == true ]; then
        #echo -n " -DHAVE_LAMMPS -Ilammps_install" >>$FPM_RSP
        echo -n " -DHAVE_LAMMPS" >>$FPM_RSP
    fi
    echo -n " -Ifutile-tmp/include\"" >>$FPM_RSP
    #--c-flag
    echo -n " --c-flag \"-Isrc -DHAVE_MKL $CFLAGS\"" >>$FPM_RSP
    #--cxx-flag
    echo -n " --cxx-flag \"-Isrc $CXXFLAGS -DQSC_STANDALONE" >>$FPM_RSP
    if [ "$decision_LAMMPS" == true ]; then
        echo -n " -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX=1 -DLAMMPS_LIB_MPI -DHAVE_LAMMPS -Ilammps_install" >>$FPM_RSP
    fi
    echo -n " \"" >>$FPM_RSP
    #--link-flag
    echo -n " --link-flag \"" >>$FPM_RSP
    if [ "$decision_SPGLIB" == true ]; then
        echo -n " -Lspglib_install" >>$FPM_RSP
    fi
    if [ "$decision_LAMMPS" == true ]; then
        echo -n " -Llammps_install" >>$FPM_RSP
    fi
    if [ "$decision_BigDFT" == true ]; then
        echo -n " -Lbigdft_install/lib" >>$FPM_RSP
    fi
    echo -n " -L$MKLPATH -Lfutile-tmp/lib $LDFLAGS\"" >>$FPM_RSP
}
write_fpm_rsp_build() {
    write_fpm_rsp_common build
    #--verbose
    echo " --verbose" >>$FPM_RSP
}
write_fpm_rsp_test() {
    write_fpm_rsp_common test
    echo >>$FPM_RSP
}
write_fpm_rsp() {
    FPM_RSP=fpm.rsp
    rm -f $FPM_RSP
    touch $FPM_RSP
    echo "# note quotes are generally required on options, must use optional" >>$FPM_RSP
    echo "# argument keys like --target" >>$FPM_RSP
    write_fpm_rsp_build
    write_fpm_rsp_test
}
#*****************************************************************************************
write_fpm_rsp_futile_install() {
    write_fpm_rsp_common install
    #--prefix
    echo " --prefix ." >>$FPM_RSP
}
write_fpm_rsp_futile() {
    FPM_RSP=fpm.rsp
    rm -f $FPM_RSP
    touch $FPM_RSP
    echo "# note quotes are generally required on options, must use optional" >>$FPM_RSP
    echo "# argument keys like --target" >>$FPM_RSP
    write_fpm_rsp_build
    write_fpm_rsp_futile_install
}
#*****************************************************************************************
