// Blocks visitors
//
// Copyright (c) 2010, Commissariat a l'Energie Atomique
// Eric Bainville, Mar 2010

#include <stdio.h>

#ifndef visitors_h
#define visitors_h

// Defines a callback function for a block visitor.  It will be called
// for all visited blocks, with the index of IA,IB (in 0..P-1), and the
// provided ARG.  IA multiple of 2, and IB multiple of 4.
typedef void (* BlockVisitor_2x4_Proc)(void * arg,int ia,int ib);

// For visitors with ZZZData structures, Call initZZZData() before running
// the pattern, then cleanupZZZData() after running the pattern to display stats and cleanup.

// Visitor doing nothing

void Nop_visitor(void * arg,int ia,int ib);

// Simple print visitor, display the blocks it traverses.

void Print_visitor(void * arg,int ia,int ib);

// Visitor checking the calling pattern verifies the constraints.

typedef struct sCheckPatternData
{
  int p;
  int * y;  // Counts occurrences for each block
  int invalid; // Counts invalid occurrences
} CheckPatternData;

void initCheckPatternData(CheckPatternData * data,int p);
void cleanupCheckPatternData(CheckPatternData * data);
void CheckPattern_visitor(void * arg,int ia,int ib);

// Visitor drawing the calling pattern in an EPS file.

typedef struct sTracePatternData
{
  int p;
  int length;  // length of a cell (PS units)
  int la,lb;  // Previous position, or (-1,-1)
  int index; // Current index
  int pass; // Current pass ; call pattern two times with pass=0 then pass=1
  FILE * f; // Output file
} TracePatternData;

void initTracePatternData(TracePatternData * data,int p,const char * filename);
void cleanupTracePatternData(TracePatternData * data);
void TracePattern_visitor(void * arg,int ia,int ib);

// Visitor computing blocks

typedef struct sComputeData
{
  int p; // Y width and height and leading dimension, and PA,PB width
  int s; // Slice dimension, PA,PB height
  const double * pa; // Packed A columns (S*P)
  const double * pb; // Packed B columns (S*P)
  double * y;  // Matrix receiving the result (P*P)
} ComputeData;

void initComputeData(ComputeData * data,int p,int s,const double * pa,const double * pb,double * y);
void cleanupComputeData(ComputeData * data);
void Compute_visitor(void * arg,int ia,int ib);

// Visitor updating the transposed block from the visited block

typedef struct sTransposeData
{
  int p;
  double * y;  // Matrix receiving the result, size P*P, leading dimension P.
} TransposeData;

void initTransposeData(TransposeData * data,int p,double * y);
void cleanupTransposeData(TransposeData * data);
void Transpose_visitor(void * arg,int ia,int ib);

#endif // #ifndef visitors_h
