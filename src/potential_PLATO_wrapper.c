#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#ifdef PARALLEL
//#include <mpi.h>
//#endif
#define  _MAIN_
#include "tb2.h"
//#include "tb2_scalapack.h"

//int Cblacs_pinfo (int *, int *);

//int Banner ();
extern char *StandardFileName (char *, char *);
int SetUp ();
int DoSCF (int);
//int DoTest ();
//int DoRelax ();
//int DoMD ();
//int DoUserDefined ();
int WriteConfig ();
//int WriteBondLengths ();
//int WriteBondAngles ();
double WallTime ();
double CPUTime ();
int LogMemory (int, int, int);
int LogTime (int, int);
//void WriteBader ();
//#ifdef HAVE_LIBSCALAPACK
//int BLACS_Finalise();
//#endif
//#ifdef STEERING
//int InitialiseREGSteerer();
//int FinaliseREGSteerer();
//int ReadParameters ();
//int WriteParameters();
//#endif

//******************************************************************************
void set_plato_parameters_(double *ResidueTol_t) {
    ResidueTol=*ResidueTol_t;
}
//******************************************************************************
int
plato_initialize_(int *nat_p,double *cv_p,double *rat_p) //(int argc, char **argv)
{
//  double tWall, tCPU;
  char FileName[100];
  int i, j;
  
//#ifdef PARALLEL
///* initialise the BLACS communication routines */
//  MPI_Init (&argc, &argv);
//#ifdef HAVE_LIBSCALAPACK
//  Cblacs_pinfo (&rank, &size);
//#endif
//  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
//  MPI_Comm_size (MPI_COMM_WORLD, &size);
//#else
  rank = 0;
  size = 1;
//#endif
  
  /* Set clock ticking */  
  WallTime ();
  CPUTime ();

  /* Write banner */
//  Banner ();

  /* Get job name */
  JobName = (char *) malloc (1024 * sizeof (char));
  GenericFileName = (char *) malloc (1024 * sizeof (char));

//  if (argc == 2)
//    {
//      strcpy (JobName, argv[1]);
//      KPointParFlag = 0;
//    }
//  else if (argc == 3)
//    {
//      if (strcmp (argv[1], "-kptpar") == 0)
//	KPointParFlag = 1;
//      else 
//	{
//	  printf ("Usage tb2 [-kptpar] <job name>\n");
//#ifdef PARALLEL
//	  MPI_Finalize ();
//#endif
//	  exit(0);
//	}
//      strcpy (JobName, argv[2]);
//    }
//  else
//    {
//      KPointParFlag = 0;
//      printf ("\nPlease enter a job name >\n");
//      scanf ("%s", JobName);
//    }
    strcpy (JobName,"plato");
  strcpy (GenericFileName, JobName);

  /* Open the output file */

  if (rank == 0)
    {
      if ((outfp = fopen (StandardFileName (FileName, "out"), "a")) == 0)
	{
	  printf ("Unable to open output file %s.\n", FileName);
	  exit (1);
	}
    }

/* Setup the program */
  SetUp ();


  *nat_p=NAtom;
  for(i=0;i<3*(*nat_p);i++) rat_p[i]=Pos[i]*0.52917721;
  for(j=0;j<3;j++) for(i=0;i<3;i++) cv_p[i+3*j]=CellVec[j][i]*0.52917721;
  //printf("REZA  %3d  %3d  %3d\n",CellRepeat[0],CellRepeat[1],CellRepeat[2]);
  //printf("CELLVEC1  %24.15f  %24.15f  %24.15f\n",CellVec[0][0],CellVec[0][1],CellVec[0][2]);
  //printf("CELLVEC2  %24.15f  %24.15f  %24.15f\n",CellVec[1][0],CellVec[1][1],CellVec[1][2]);
  //printf("CELLVEC3  %24.15f  %24.15f  %24.15f\n",CellVec[2][0],CellVec[2][1],CellVec[2][2]);
//#ifdef STEERING
//  ReadParameters();
//#endif


//#ifdef STEERING
//  /* Set up steerer */
//  InitialiseREGSteerer();
//#endif 

  return 1;
}
//******************************************************************************
int plato_energy_forces_(int *nat,double *rat,double *epot,double *fat) {
    int i;
//  /* Do job */
//  switch (job)
    {
//    case 0:
//      DoTest ();
//      break;
//    case 1:
    for(i=0;i<3*(*nat);i++) {
        Pos[i]=rat[i]/0.52917721;
    }
      DoSCF (1);
    *epot=Etotal*13.6058;
    for(i=0;i<3*(*nat);i++) {
        fat[i]=Ftotal[i]*(13.6058/0.52917721);
    }
      WriteConfig ();
//      if (rank == 0)
//        {
//          WriteBondLengths ();
//          WriteBondAngles ();
//        }
//      break;
//    case 2:
//      DoRelax ();
//      break;
//    case 3:
//      DoMD ();
//      break;
//    case 99:
//      DoUserDefined ();
//      break;
//    default:
//      printf ("Unknown job.\n");
//#ifdef PARALLEL
//      MPI_Finalize ();
//#endif
//      exit (1);
    }
  return 1;
}
//******************************************************************************
int plato_finalize_() {
//#ifdef STEERING
//  /* Write out the steered parameters for restarting */
//  if (rank == 0)
//    WriteParameters ();
//#endif

//  if (WriteBaderFlag)
//    {
//      WriteBader ();
//    }

  /* Write out the total wall time taken to complete the simulation,
     and then close the output file */
  LogMemory (0, 0, 1);
  LogTime (0, 2);
//  tWall = WallTime ();
//  tCPU = CPUTime ();
  if (rank == 0)
    {
//      fprintf (outfp, "\n\nWall time = %9.1fs  %9.1fm  %9.1fh\n", tWall, tWall / 60.0, tWall / 3600.0);
//      fprintf (outfp, "CPU time  = %9.1fs  %9.1fm  %9.1fh\n", tCPU, tCPU / 60.0, tCPU / 3600.0);
      fclose (outfp);
    }

  /* Finalise the BLACS stuff for the ScaLAPACK diagonaliser */
//#ifdef HAVE_LIBSCALAPACK
//      if (Diagmethod > 5)
//	BLACS_Finalise();
//#endif

//#ifdef STEERING
//  /* Free the steerer variables */
//  FinaliseREGSteerer();
//#endif 

//#ifdef PARALLEL
//  MPI_Finalize ();
//#endif

  return 0;
}
//******************************************************************************
