#include <config.h>

#define _GNU_SOURCE

void FC_FUNC_(mpi_ialltoallv,MPI_IALLTOALLV)(void *sendbuf,void *sendcounts,void *senddspls,void *sendtype,
					     void *recvbuf,void *recvcounts,void *recvdspls,void *recvtype,
					     void *comm,void *request,void *ierr)
{
  FC_FUNC_(mpi_alltoallv,MPI_ALLTOALLV)(sendbuf,sendcounts,senddspls,sendtype,
					recvbuf,recvcounts,recvdspls,recvtype,
					comm,ierr);
}


void FC_FUNC_(mpi_iallreduce,MPI_IALLREDUCE)(void *sendbuf,void *recvbuf,void *count,void *type,void *op,
					     void *comm,void *request,void *ierr)
{
  FC_FUNC_(mpi_allreduce,MPI_ALLREDUCE)(sendbuf,recvbuf,count,type,op,comm,ierr);
}
