#!/bin/bash

source ./utils/bash/util.sh
source ./utils/bash/inquire_SPGLIB.sh
source ./utils/bash/inquire_BIGDFTROOT.sh
source ./utils/bash/inquire_LAMMPS.sh
source ./utils/bash/inquire_FLAGS.sh
source ./utils/bash/write_fpm_toml.sh
source ./utils/bash/write_fpm_rsp.sh

FLAME_DIR=.
FUTILE_DIR=$FLAME_DIR/futile-suite-1.9.1

BLUE='\033[0;34m'
WHITE='\033[0;37m'
RED='\033[0;31m'
GREEN='\033[0;32m'
COLOR_OFF='\033[0m'

rm -rf bigdft_install
mkdir bigdft_install
mkdir bigdft_install/include
mkdir bigdft_install/lib

rm -rf spglib_install
mkdir spglib_install

rm -rf lammps_install
mkdir lammps_install

main() {
    decision_SPGLIB=$(inquire_yesno "Do you like to link FLAME with SPGLIB?")
    if [ "$decision_SPGLIB" == true ]; then
        inquire_SPGLIB SPGLIB
    fi
    decision_BigDFT=$(inquire_yesno "Do you like to link FLAME with BigDFT PSolver?")
    if [ "$decision_BigDFT" == true ]; then
        inquire_BIGDFTROOT BIGDFTROOT
    fi
    decision_LAMMPS=$(inquire_yesno "Do you like to link FLAME with LAMMPS?")
    if [ "$decision_LAMMPS" == true ]; then
        inquire_LAMMPS LAMMPS
    fi
    inquire_FLAGS
    write_fpm_toml
    write_fpm_rsp
    source ./utils/bash/prep_futile.sh
    cd futile-tmp
    write_fpm_toml_futile
    write_fpm_rsp_futile
    cd ..
    source ./utils/bash/get_input_variables_definition_c.sh
}

main

