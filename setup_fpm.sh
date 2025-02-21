#!/bin/bash

BLUE='\033[0;34m'
WHITE='\033[0;37m'
RED='\033[0;31m'
GREEN='\033[0;32m'
COLOR_OFF='\033[0m'

source ./utils/bash/util.sh
source ./utils/bash/inquire_SPGLIB.sh
source ./utils/bash/inquire_BIGDFTROOT.sh
source ./utils/bash/inquire_LAMMPS.sh
source ./utils/bash/inquire_FLAGS.sh
source ./utils/bash/write_fpm_toml.sh
source ./utils/bash/write_fpm_rsp.sh
source ./utils/bash/fpm_check_clean.sh

FLAME_DIR=.
FUTILE_DIR=$FLAME_DIR/futile-suite-1.9.1

helpFunction() {
    #echo ""
    echo -e "${RED}Usage:${COLOR_OFF} ${GREEN}$0 -f filename -h${COLOR_OFF}"
    echo -e "    ${GREEN}-f filename${COLOR_OFF} must contain required parameters, examples in arch!"
    echo -e "    ${GREEN}-h${COLOR_OFF} shows these descriptions!"
    echo -e "    Only one of arguments at a time: either -f or -h"
    echo -e "    If executed with no argument, required parameters from screen!"
    exit 1 #exit script after printing help
}

inquire_ALL() {
    decision_clean=false
    check_if_it_is_clean
    if [ "$is_it_clean" == false ]; then
        decision_clean=$(inquire_yesno "setup_fpm.sh is previously executed, do you like to clean?")
    fi
    if [ "$decision_clean" == true ]; then
        clean_for_fpm
    fi
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
}
main() {
    if [ "$arg_f" == false ]; then
        inquire_ALL
    fi
    if [ "$CLEAN_YESNO" == "yes" ]; then
        decision_clean=true
        clean_for_fpm
    else
        decision_clean=false
    fi
    if [ "$SPGLIB_YESNO" == "yes" ]; then
        decision_SPGLIB=true
        make_spglib_install
    else
        decision_SPGLIB=false
    fi
    if [ "$BigDFT_YESNO" == "yes" ]; then
        decision_BigDFT=true
        make_bigdft_install
    else
        decision_BigDFT=false
    fi
    if [ "$LAMMPS_YESNO" == "yes" ]; then
        decision_LAMMPS=true
        make_lammps_install
    else
        decision_LAMMPS=false
    fi
    write_fpm_toml
    write_fpm_rsp
    source ./utils/bash/prep_futile.sh
    cd futile-tmp
    write_fpm_toml_futile
    write_fpm_rsp_futile
    cd ..
    source ./utils/bash/get_input_variables_definition_c.sh
}

arg_h=false
arg_f=false
while getopts "f:h" arg; do
    #echo "${arg}"
    case "${arg}" in
        f)
            config_filename=${OPTARG}
            arg_f=true
            ;;
        h)
            arg_h=true
            ;;
        ?)
            echo -e "${RED}ERROR:${COLOR_OFF} something wrong with argument list!"
            echo -e "       please run :${GREEN}./setup_fpm.sh -h${COLOR_OFF} for help!"
            exit 1
            ;;
    esac
done
if [ "$arg_h" == true ]; then
   helpFunction
fi
#print helpFunction in case parameters are empty
if [ -z "$config_filename" ] && [ "$arg_f" == true ]; then
   echo -e "${RED}ERROR: input filename is missing!${COLOR_OFF}";
   helpFunction
fi
#echo "$config_filename"
if [ "$config_filename" == "-h" ]; then
    echo -e "${RED}ERROR:${COLOR_OFF} ${GREEN}-f${COLOR_OFF} requires an argument, namely the filename!"
    echo -e "       You entered ${GREEN}-h${COLOR_OFF}, do you mean help?"
    echo -e "       please run :${GREEN}./setup_fpm.sh -h${COLOR_OFF} for help!"
    exit 1
fi
if [ "$arg_f" == true ]; then
    if [ -f $config_filename ]; then
        source "$config_filename"
    else
        echo -e "${RED}ERROR:${COLOR_OFF} $config_filename does not exist!"
        exit 1
    fi
fi
main

