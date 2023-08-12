#!/bin/bash

inquire_FLAGS() {
    echo "Provide the following info: compilers and flags for compilation and linking,"
    read -p "FC=" FC
    read -p "CC=" CC
    read -p "CXX=" CXX
    read -p "FFLAGS=" FFLAGS
    read -p "CFLAGS=" CFLAGS
    read -p "CXXFLAGS=" CXXFLAGS
    read -p "LDFLAGS=" LDFLAGS
    read -p "MKLPATH=" MKLPATH
}
