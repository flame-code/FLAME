#!/bin/bash

inquire_SPGLIB() {
    echo
    echo -e "FLAME needs ${GREEN}spglib_f08.f90 ${COLOR_OFF}and ${GREEN}libsymspg.a${COLOR_OFF} from SPGLIB."
    echo -e "Depending on your installation of SPGLIB, these two files"
    echo -e "may be on the same place or two locations."
    echo -e "For this reason, paths to these two files need to be inquired independently,"
    echo "provide the full path of the file spglib_f08.f90:"
    local SPGLIB_F08=""
    read -p "SPGLIB_F08=" SPGLIB_F08
    echo "provide the full path of the file libsymspg.a:"
    local LIBSYMSPG=""
    read -p "LIBSYMSPG=" LIBSYMSPG

    error_SPGLIB_F08=false
    if ! [ -f $SPGLIB_F08 ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find spglib_f08.f90 needed to link with SPGLIB!"
        error_SPGLIB_F08=true
    else
        echo -e "file found: ${GREEN}$SPGLIB_F08${COLOR_OFF}"
    fi
    error_LIBSYMSPG=false
    if ! [ -f $LIBSYMSPG ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find libsymspg.a needed to link with SPGLIB!"
        error_LIBSYMSPG=true
    else
        echo -e "file found: ${GREEN}$LIBSYMSPG${COLOR_OFF}"
    fi
    if [[ $error_SPGLIB_F08 == true  ||  $error_LIBSYMSPG == true ]] ; then
        read -p "Do you want to continue without linking to SPGLIB?" answer
        if [ "$answer" == "n" ]; then
            echo -e "${RED}setup quits!${COLOR_OFF}"
            echo -e "${BLUE}Please run the setup again when both aforementioned files are available!${COLOR_OFF}"
            exit
        else
            decision_SPGLIB=false
        fi
    else
        ln -s $SPGLIB_F08 src/spglib_f08.f90
        ln -s $LIBSYMSPG spglib_install/libsymspg.a
    fi
}
