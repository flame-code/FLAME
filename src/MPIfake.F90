        subroutine MPI_INIT(ierr)
        return
        end

        subroutine MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        return
        end

        subroutine MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        return
        end

        subroutine MPI_FINALIZE(ierr)
        return
        end

        subroutine MPI_ISEND(abuf,i,MPI_DOUBLE_PRECISION,jproc,mtag,MPI_COMM_WORLD,irequest,ierr)
        return
        end

        subroutine MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,flag,status_mpi,ierr)
        return
        end

        subroutine MPI_RECV(re_send,i,MPI_DOUBLE_PRECISION,i_source,i_tag,MPI_COMM_WORLD,status_mpi,ierr)
        return
        end

	subroutine MPI_barrier(MPI_COMM_WORLD,ierr)
        return
        end

        subroutine MPI_BCAST(nlmin,i,MPI_INTEGER,i_source,MPI_COMM_WORLD,ierr)
        return
        end
