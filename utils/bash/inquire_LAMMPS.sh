#!/bin/bash

inquire_LAMMPS() {
    echo
    echo -e "FLAME needs ${GREEN}LAMMPS header files ${COLOR_OFF}and ${GREEN}liblammps_*.a${COLOR_OFF} from LAMMPS."
    echo -e "Depending on your installation of LAMMPS, the library filename"
    echo -e "may be different."
    echo -e "For this reason, paths to directory of the header files and the library files are inquired independently,"
    echo "provide the full path of the directory of the header files:"
    local LAMMPS_HEADER=""
    read -p "LAMMPS_HEADER=" LAMMPS_HEADER
    echo "provide the full path of the library file, i.e. liblammps_*.a:"
    local LIBLAMMPS=""
    read -p "LIBLAMMPS=" LIBLAMMPS

    error_LAMMPS_HEADER=false
    for file in library.h lammps.h atom.h input.h fix.h fix_external.h compute.h modify.h error.h pointers.h lmptype.h
    do
        if ! [ -f $LAMMPS_HEADER/$file ]; then
            echo -e "${RED}ERROR:${COLOR_OFF} cannot find $file needed to link with LAMMPS!"
            error_LAMMPS_HEADER=true
        else
            echo -e "file found: ${GREEN}$LAMMPS_HEADER/$file ${COLOR_OFF}"
        fi
    done

    error_LIBLAMMPS=false
    if ! [ -f $LIBLAMMPS ]; then
        echo -e "${RED}ERROR:${COLOR_OFF} cannot find LAMMPS library needed to link with LAMMPS!"
        error_LIBLAMMPS=true
    else
        echo -e "file found: ${GREEN}$LIBLAMMPS${COLOR_OFF}"
    fi
    if [[ $error_LAMMPS_HEADER == true  ||  $error_LIBLAMMPS == true ]] ; then
        read -p "Do you want to continue without linking to LAMMPS?" answer
        if [ "$answer" == "n" ]; then
            echo -e "${RED}setup quits!${COLOR_OFF}"
            echo -e "${BLUE}Please run the setup again when both aforementioned files are available!${COLOR_OFF}"
            exit
        else
            decision_LAMMPS=false
        fi
    else
        for file in library.h lammps.h atom.h input.h fix.h fix_external.h compute.h modify.h error.h pointers.h lmptype.h
        do
            ln -s $LAMMPS_HEADER/$file lammps_install/$file
        done
        ln -s $LIBLAMMPS lammps_install/liblammps_mpi.a
    fi
}
