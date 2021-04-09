/*
!> @file
!!  Routines to access easily to the filesystem and to other C goodies.
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/

#include <config.h>

#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <setjmp.h>
#ifdef HAVE_TIME_H
#include <time.h>
#endif

#include "utils.h"

/*#include "smpi.h"*/

#ifndef HAVE_STRNDUP
char* strndup(const char *src, size_t len);
#endif

static jmp_buf jb_main,jb_child;

#ifndef HAVE_CLOCK_GETTIME
#define CLOCK_REALTIME 0
static int clock_gettime(int clk_id, struct timespec *tp)
{
  struct timeval now;
  int rv;

  tp->tv_sec = 0;
  tp->tv_nsec = 0;

  rv = gettimeofday(&now, NULL);
  if (rv != 0)
    return rv;
        
  tp->tv_sec = now.tv_sec;
  tp->tv_nsec = now.tv_usec * 1000;

  return 0;
}
#endif

#ifndef HAVE_ALIGNED_ALLOC
void *aligned_alloc(size_t alignment, size_t size)
{
  void* ptr = NULL;
  //printf(" size, alignment %zu, %zu \n",size,alignment);
  posix_memalign(&ptr,alignment,size);
  //printf("ptr %p \n",ptr);
  return ptr;
}
#endif

void FC_FUNC(nanosec,NANOSEC)(unsigned long long int * t){
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  *t = time.tv_sec;
  *t *= 1000000000;
  *t += time.tv_nsec;
}


void FC_FUNC(csleep,CSLEEP)(int * secs){
    unsigned int seconds;
    seconds = (unsigned int) *secs;
    sleep(seconds);
  }

void FC_FUNC(getaddress, GETADDRESS)(void *ptr,char *address, int *lgaddress,
			     int* status)
{
  char buff[50]; //test buffer to check the length
  int lgt;

  memset(address,' ', sizeof(char) * (*lgaddress));

  lgt=sprintf(buff,"%p",(void*)ptr);
  if (lgt < *lgaddress)
    {
      *status = 0;
      memcpy(address, buff, sizeof(char) * lgt);
    }
  else
    *status = 1;
  return;

  //printf("\n test address = %p %d; \n", (void*)ptr,lgt);
  //return;
}

void FC_FUNC(getdir, GETDIR)(const char *dir, int *lgDir,
                             char *out, int *lgOut,
                             int *status)
{
  char *path;
  struct stat sb;
  int lgCpy;
#if defined _WIN32 || defined __CYGWIN__
  const char kPathSeparator = '\\';
  int lnsep=2;
#else
  const char kPathSeparator = '/';
  int lnsep=1;
#endif

  *status = 1;
  memset(out, ' ', sizeof(char) * (*lgOut));

  path = strndup(dir, (size_t)*lgDir);

  if (stat(path, &sb) == 0)
    {
      free(path);
      /* Dir exists already. */
      if (S_ISDIR(sb.st_mode))
        {
          *status = 0;
          lgCpy = ((*lgDir > *lgOut - 1)?*lgOut - 1:*lgDir);
          memcpy(out, dir, sizeof(char) * lgCpy);
          /* Add a separator if not already present */
          if (out[lgCpy-lnsep] != kPathSeparator) { out[lgCpy] = kPathSeparator; };
        }
      else
        *status = 1;
      return;
    }

  /* Try to create it. */
#if defined _WIN32 || defined __CYGWIN__
  if (mkdir(path) != 0)
#else
  if (mkdir(path, 0755) != 0)
#endif
    {
      if (errno != EEXIST)
	{
           free(path);
           *status = 2;
           return;
        }
     }

  free(path);
  lgCpy = ((*lgDir > *lgOut - 1)?*lgOut - 1:*lgDir);
  memcpy(out, dir, sizeof(char) * lgCpy);
  /* Add a '/' if not already present */
  if (out[lgCpy-lnsep] != kPathSeparator) { out[lgCpy] = kPathSeparator; };
  *status = 0;
  return;
}

void FC_FUNC(bindfree, BINDFREE)(long long int * fromadd)
{
  void *from = (void*)*fromadd;
  /*printf("\n test long address = %p %lli\n", (void*)from,*fromadd);*/
  free(from);
}

void mybindfree(long long int * fromadd)
{
  void *from = (void*)*fromadd;
  /*printf("\n test long address = %p %lli\n", (void*)from,*fromadd);*/
#ifdef HAVE_SIMGRID_SHARED_ALLOCS
  smpi_shared_free(from);
#else
  free(from);
#endif
}


void FC_FUNC(callsystem, CALLSYSTEM)(const char *cmd, int *lgcmd, int *status)
{
  char *command;

  command = strndup(cmd, (size_t)*lgcmd);
  *status = system(command);
  free(command);
}


void FC_FUNC(delete, DELETE)(const char *f, int *lgF, int *status)
{
  char *path;

  path = strndup(f, (size_t)*lgF);
  *status = unlink(path);
  free(path);
}

void FC_FUNC(deldir, DELDIR)(const char *f, int *lgF, int *status)
{
  char *path;

  path = strndup(f, (size_t)*lgF);
  *status = rmdir(path);
  free(path);
}


void FC_FUNC(movefile, MOVEFILE)(const char *oldfile, int *lgoldfile, const char *newfile, int *lgnewfile, int *status)
{
  char *oldpath;
  char *newpath;

  oldpath = strndup(oldfile, (size_t)*lgoldfile);
  newpath = strndup(newfile, (size_t)*lgnewfile);
  *status = rename(oldpath,newpath);
  free(oldpath);
  free(newpath);
}

void FC_FUNC(getfilecontent, GETFILECONTENT)(void **pt, long *pt_len, const char *fname, int *ln)
{
  FILE *f;
  char *buf, *fname_;
  size_t r;
  long s;

  fname_ = strndup(fname, (size_t)*ln);
  f = fopen(fname_, "rb");
  free(fname_);

  fseek(f, 0L, SEEK_END);
  s = ftell(f);
  rewind(f);

  buf = malloc(sizeof(char) * (s + 1));
  r = fread(buf, s, 1, f);
  buf[r * s] = '\0';

  fclose(f);

  *pt = (void*)buf;
  *pt_len = s;
}

// Take a Long integer as th enumber of bythes
void FC_FUNC(memsetzero,MEMSETZERO)(void *buf, long long int *ln)
{
  size_t nbytes = *ln;

  memset(buf, 0, nbytes);
}


void FC_FUNC(copycbuffer, COPYCBUFFER)(char *to, void **cbuf, long *ln)
{
  char *from = (char*)*cbuf;

  memcpy(to, from, sizeof(char) * *ln);
}

void FC_FUNC(freecbuffer, FREECBUFFER)(void **buf)
{
  free(*buf);
}

void FC_FUNC(getprocid, GETPROCID)(int *procid)
{
  *procid = (int) getpid();
}

void FC_FUNC(stdoutistty, STDOUTISTTY)(int *itis)
{
  int fd = STDOUT_FILENO;
  *itis = isatty(fd);
}

void FC_FUNC(getjmpbufsize, GETJMPBUFSIZE)(int *bufsize)
{
  *bufsize = (int) sizeof(jb_main);
}

void FC_FUNC(setandcpyjmpbuf, SETANDCPYJMPBUF)(int *signal, FFunc_void* routine,long long int* arg, int* last_signal)
{
  *signal = setjmp(jb_main);
    //fill the jmpbuf with the correct values
    //if (!*signal){ //!*signal){
    //  printf("setjmp first time, signal, %d\n",*signal);
    //  //memcpy(*jb, jbs, sizeof(jbs));
    //  //here call lngjmp to see what happens
    //  //FC_FUNC(longjmpwrap, LONGJMPWRAP)(jb,signal);
    //}
    //else 
    //  printf("restarting, signal %d, %d \n",*signal,*last_signal);   
    if (*signal != *last_signal && routine)
      (*routine)();
}

void FC_FUNC(longjmpwrap, LONGJMPWRAP)(jmp_buf* jb,int* signal)
{
  //memcpy(jbs,*jb, sizeof(jbs));

  //printf("before longjump, signal %d\n",*signal);   
  longjmp(jb_main,*signal);
}
