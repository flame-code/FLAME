!*****************************************************************************************
module mod_flm_futile
    !-----------------------------------------------------------------
    use dictionaries, only: dictionary, dict_new, dict_free
    use dictionaries, only: operator(.is.)
    !-----------------------------------------------------------------
    use dynamic_memory, only: f_malloc_set_status
    use dynamic_memory, only: f_routine, f_release_routine
    use dynamic_memory, only: f_malloc, f_malloc0, f_free
    use dynamic_memory, only: operator(.to.)
    use dynamic_memory, only: assignment(=)
    !-----------------------------------------------------------------
    use yaml_output, only: yaml_new_document, yaml_release_document
    use yaml_output, only: yaml_flush_document
    use yaml_output, only: yaml_set_stream, yaml_set_default_stream
    use yaml_output, only: yaml_map, yaml_comment, yaml_warning
    use yaml_output, only: yaml_mapping_open, yaml_mapping_close
    use yaml_output, only: yaml_scalar
    !-----------------------------------------------------------------
    use f_utils, only: f_get_free_unit
    !-----------------------------------------------------------------
    use yaml_strings, only: yaml_toa, yaml_date_and_time_toa
    !-----------------------------------------------------------------
    !use wrapper_mpi, only: mpi_environment, MPI_COMM_WORLD
    !-----------------------------------------------------------------
    implicit none
end module mod_flm_futile
!*****************************************************************************************
