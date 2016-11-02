/*
!> @file
!!  define the padding for the malloc function
!! @author
!!    Copyright (C) 2016-2016 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>   
#include <string.h>

static size_t repad_size_f_malloc = 0;
static char envstr[256]= "";

void FC_FUNC(togglepadding,TOGGLEPADDING)(long long int * totbytes)
{
repad_size_f_malloc = (size_t) *totbytes;
memset(envstr,0,sizeof(envstr));
snprintf(envstr, sizeof(envstr), "FUTILE_REPAD_SIZE_INTERNAL=%zu", repad_size_f_malloc);
/*fprintf(stderr, "%s \n",envstr );*/
int ierr=putenv(envstr);
if (ierr != 0) fprintf(stderr, "error on putenv %d \n",ierr );
/*test of the environment 
char *restest=getenv("FUTILE_REPAD_SIZE_INTERNAL");
printf("Got environment:%s \n",restest); */
}

size_t get_padding(void)
{
  return repad_size_f_malloc;
}
