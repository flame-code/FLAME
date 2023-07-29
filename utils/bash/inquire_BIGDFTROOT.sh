#!/bin/bash

inquire_BIGDFTROOT() {
    echo "Provide path to the BigDFT installation directory!"
    local BIGDFTROOT=""
    read -p "BIGDFTROOT=" BIGDFTROOT
    echo
    #echo $BIGDFTROOT
    #find $BIGDFTROOT/install/lib -type f -name 'libatlab-1.a'
    error_BIGDFT_lib=false
    if ! [ -f $BIGDFTROOT/install/lib/libatlab-1.a ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find libatlab-1.a needed to link with BigDFT PSolver!"
        error_BIGDFT_lib=true
    else
        echo -e "file found: ${GREEN}$BIGDFTROOT/install/lib/libatlab-1.a${COLOR_OFF}"
    fi
    if ! [ -f $BIGDFTROOT/install/lib/libPSolver-1.a ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find libPSolver-1.a needed to link with BigDFT PSolver!"
        error_BIGDFT_lib=true
    else
        echo -e "file found: ${GREEN}$BIGDFTROOT/install/lib/libatlab-1.a${COLOR_OFF}"
    fi
    if $error_BIGDFT_lib ; then
        echo -e "BigDFT library files are expected to be at $BIGDFTROOT/install/lib"
    fi
    echo
    error_BIGDFT_mod=false
    if ! [ -f $BIGDFTROOT/install/include/at_domain.mod ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find at_domain.mod needed to link with BigDFT PSolver!"
        error_BIGDFT_mod=true
    else
        echo -e "file found: ${GREEN}$BIGDFTROOT/install/include/at_domain.mod${COLOR_OFF}"
    fi
    if ! [ -f $BIGDFTROOT/install/include/poisson_solver.mod ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find poisson_solver.mod needed to link with BigDFT PSolver!"
        error_BIGDFT_mod=true
    else
        echo -e "file found: ${GREEN}$BIGDFTROOT/install/include/poisson_solver.mod${COLOR_OFF}"
    fi
    if $error_BIGDFT_mod ; then
        echo -e "BigDFT module files are expected to be at $BIGDFTROOT/install/include"
    fi
    echo
    if [[ $error_BIGDFT_lib == true  ||  $error_BIGDFT_mod == true ]] ; then
        read -p "Do you want to continue without linking to BigDFT?" answer
        if [ "$answer" == "n" ]; then
            echo -e "${RED}setup quits!${COLOR_OFF}"
            echo -e "${BLUE}Please run the setup again after BigDFT installation is complete!${COLOR_OFF}"
            exit
        else
            decision_BigDFT=false
        fi
    else
        ln -s $BIGDFTROOT/install/include/at_domain.mod bigdft_install/include/at_domain.mod
        ln -s $BIGDFTROOT/install/include/poisson_solver.mod bigdft_install/include/poisson_solver.mod
        ln -s $BIGDFTROOT/install/lib/libatlab-1.a bigdft_install/lib/libatlab-1.a
        ln -s $BIGDFTROOT/install/lib/libPSolver-1.a bigdft_install/lib/libPSolver-1.a
    fi
    #echo $BIGDFTROOT
}
