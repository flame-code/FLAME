#!/bin/bash


check_if_it_is_clean() {
    is_it_clean=true
    if [ -d futile-tmp ]; then
        is_it_clean=false
    fi
    if [ -d spglib_install ]; then
        is_it_clean=false
    fi
    if [ -d bigdft_install ]; then
        is_it_clean=false
    fi
    if [ -d lammps_install ]; then
        is_it_clean=false
    fi
    if [ -f fpm.rsp ]; then
        is_it_clean=false
    fi
    if [ -f fpm.toml ]; then
        is_it_clean=false
    fi
}

clean_for_fpm() {
    echo -e "${RED}Cleaning files and directories previosly created by setup_fpm.sh!${COLOR_OFF}"
    if [ -d futile-tmp ]; then
        echo -e -n "Deleting directory ${BLUE}futile-tmp${COLOR_OFF} "
        rm -rf futile-tmp
        echo "done."
    fi
    if [ -d spglib_install ]; then
        echo -e -n "Deleting directory ${BLUE}spglib_install${COLOR_OFF} "
        rm -rf spglib_install
        echo "done."
    fi
    if [ -f src/spglib_f08.f90 ]; then
        echo -e -n "Deleting file ${BLUE}src/spglib_f08.f90${COLOR_OFF} "
        rm -f src/spglib_f08.f90
        echo "done."
    fi
    if [ -d bigdft_install ]; then
        echo -e -n "Deleting directory ${BLUE}bigdft_install${COLOR_OFF} "
        rm -rf bigdft_install
        echo "done."
    fi
    if [ -d lammps_install ]; then
        echo -e -n "Deleting directory ${BLUE}lammps_install${COLOR_OFF} "
        rm -rf lammps_install
        echo "done."
    fi
    if [ -f fpm.rsp ]; then
        echo -e -n "Deleting file ${BLUE}fpm.rsp${COLOR_OFF} "
        rm -f fpm.rsp
        echo "done."
    fi
    if [ -f fpm.toml ]; then
        echo -e -n "Deleting file ${BLUE}fpm.toml${COLOR_OFF} "
        rm -f fpm.toml
        echo "done."
    fi
}
