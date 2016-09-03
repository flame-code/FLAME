// Blocks iterator patterns
//
// Copyright (c) 2010, Commissariat a l'Energie Atomique
// Eric Bainville, Mar 2010

#ifndef patterns_h
#define patterns_h

#include "visitors.h"

// Defines a block iterator.  Executing it should call the
// provided block visitor for all blocks to process.
// P is the number of columns of the two matrices.
// ARG are the pattern private parameters.
typedef void (* BlockPattern_2x4_Proc)(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg);

// Patterns

// Regular loop on column groups
typedef struct sColumnsPatternData
{
  int nb; // Block size (must be multiple of 4)
} ColumnsPatternData;
void columns_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg);

// Regular loop on row groups
typedef struct sRowsPatternData
{
  int na; // Block size (must be multiple of 2)
} RowsPatternData;
void rows_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg);

// Block loop
typedef struct sBlocksPatternData
{
  int na,nb; // Block size (must be multiple of 2 and 4 resp.)
  int order; // Loop order (0 or 1)
} BlocksPatternData;
void blocks_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg);

// Peano curve
void peano_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg);

// Identify pattern and its arguments into string S,LEN
void idendifyPattern(BlockPattern_2x4_Proc pattern,void * pattern_arg,char * s,int len);

#endif // #ifndef patterns_h
