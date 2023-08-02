#!/bin/bash


write_fpm_toml() {
    FPM_TOML=fpm.toml
    rm -f $FPM_TOML
    touch $FPM_TOML
    echo "name = \"FLAME\"" >>$FPM_TOML
    echo "license = \"GPLv3\"" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[build]" >>$FPM_TOML
    echo "auto-executables = true" >>$FPM_TOML
    echo "auto-tests = true" >>$FPM_TOML
    echo "auto-examples = true" >>$FPM_TOML
    echo -n "link = [" >>$FPM_TOML
    if [ "$decision_BigDFT" == true ]; then
        echo -n "\"PSolver-1\",\"atlab-1\"," >>$FPM_TOML
    fi
    if [ "$decision_LAMMPS" == true ]; then
        echo -n "\"lammps_mpi\"," >>$FPM_TOML
    fi
    if [ "$decision_SPGLIB" == true ]; then
        echo -n "\"symspg\"," >>$FPM_TOML
    fi
    echo -n "\"FUTILE\",\"yaml\"" >>$FPM_TOML
    #echo -n ",\"iomp5\",\"m\",\"dl\",\"rt\"" >>$FPM_TOML
    echo "]" >>$FPM_TOML
    echo "external-modules = [\"yaml_parse\",\"dictionaries\",\"dynamic_memory\",\"yaml_output\"," >>$FPM_TOML
    echo "    \"time_profiling\",\"wrapper_linalg\",\"poisson_solver\",\"futile\",\"wrapper_mpi\"," >>$FPM_TOML
    echo "    \"at_domain\",\"module_fft_sg\",\"mpi\",\"mod_processorsm\",\"m_siesta_init\",\"m_siesta_move\"," >>$FPM_TOML
    echo "    \"m_siesta_end\",\"m_siesta_forces\",\"siesta_geom\",\"atmfuncs\",\"parallel\",\"ifport\",\"m_energies\"," >>$FPM_TOML
    echo "    \"m_forces\",\"siesta_options\",\"basis_types\",\"atm_types\",\"f_utils\",\"yaml_strings\"," >>$FPM_TOML
    echo -n "    \"interface_tinker\",\"Poisson_Solver\"" >>$FPM_TOML
    if [ "$decision_SPGLIB" == false ]; then
        echo -n ",\"spglib_f08\"" >>$FPM_TOML
    fi
    echo "]" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[install]" >>$FPM_TOML
    echo "library = false" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[[executable]]" >>$FPM_TOML
    echo "name = \"flame\"" >>$FPM_TOML
    echo "source-dir = \"src\"" >>$FPM_TOML
    echo "main = \"flame.F90\"" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[[test]]" >>$FPM_TOML
    echo "name=\"fpm\"" >>$FPM_TOML
    echo "source-dir=\"tests-fpm\"" >>$FPM_TOML
    echo "main = \"main.F90\"" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[fortran]" >>$FPM_TOML
    echo "implicit-typing = true" >>$FPM_TOML
    echo "implicit-external = true" >>$FPM_TOML
    echo "source-form = \"default\"" >>$FPM_TOML
}

write_fpm_toml_futile() {
    FPM_TOML=fpm.toml
    rm -f $FPM_TOML
    touch $FPM_TOML
    echo "name = \"FUTILE\"" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[install]" >>$FPM_TOML
    echo "library = true" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[library]" >>$FPM_TOML
    echo "source-dir = \".\"" >>$FPM_TOML
    echo "include-dir = \".\"" >>$FPM_TOML
    echo >>$FPM_TOML
    echo "[preprocess]" >>$FPM_TOML
    echo "[preprocess.cpp]" >>$FPM_TOML
    echo "macros = [\"HAVE_CONFIG_H\"]" >>$FPM_TOML
}
