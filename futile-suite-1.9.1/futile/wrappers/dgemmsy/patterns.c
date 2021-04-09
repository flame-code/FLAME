//! @file
//!  Blocks iterator patterns
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 
//! Eric Bainville, Mar 2010


#include "patterns.h"

void columns_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg)
{
  ColumnsPatternData * data = (ColumnsPatternData *)arg;
  int p2 = 4 + (p>>1); // vertical size of band
  int ia,ib0,ib1;

  // Loop on column blocks
  for (ib0=0;ib0<p;ib0+=data->nb)
  {
    // Then on rows
    for (ia=0;ia<p;ia+=2)
    {
      // And on columns inside each block
      for (ib1=0;ib1<data->nb;ib1+=4)
      {
	int ib = ib0 + ib1;
	if (ia >= p || ib >= p) continue;
	if ( ((p + ia - ib)%p) >= p2 ) continue;
	visitor(visitor_arg,ia,ib);
      }
    }
  }
}

void rows_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg)
{
  RowsPatternData * data = (RowsPatternData *)arg;
  int p2 = 4 + (p>>1); // vertical size of band
  int ia0,ia1,ib;

  // Loop on row blocks
  for (ia0=0;ia0<p;ia0+=data->na)
  {
    // Then on columns
    for (ib=0;ib<p;ib+=4)
    {
      // And on rows inside each block
      for (ia1=0;ia1<data->na;ia1+=2)
      {
	int ia = ia0 + ia1;
	if (ia >= p || ib >= p) continue;
	if ( ((p + ia - ib)%p) >= p2 ) continue;
	visitor(visitor_arg,ia,ib);
      }
    }
  }
}

void blocks_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg)
{
  BlocksPatternData * data = (BlocksPatternData *)arg;
  int p2 = 4 + (p>>1); // vertical size of band
  int ia0,ib0,ia1,ib1;

#define BLOCKS_ITERATION \
      int ia = ia0 + ia1; \
      int ib = ib0 + ib1; \
      if (ia >= p || ib >= p) continue; \
      if ( ((p + ia - ib)%p) >= p2 ) continue; \
      visitor(visitor_arg,ia,ib);

  switch (data->order)
  {
  case 0:
    for (ia0=0;ia0<p;ia0+=data->na) for (ib0=0;ib0<p;ib0+=data->nb)
    for (ia1=0;ia1<data->na;ia1+=2) for (ib1=0;ib1<data->nb;ib1+=4)
    { BLOCKS_ITERATION }
    break;
  case 1:
    for (ib0=0;ib0<p;ib0+=data->nb) for (ia0=0;ia0<p;ia0+=data->na)
    for (ib1=0;ib1<data->nb;ib1+=4) for (ia1=0;ia1<data->na;ia1+=2)
    { BLOCKS_ITERATION }
    break;
  }
}

// Put all call arguments into one struct for recursive calls
struct PeanoData
{
  int p;
  BlockVisitor_2x4_Proc visitor;
  void * visitor_arg;
};

// Corner bit 0: 0=left, 1=right
// Corner bit 1: 0=top, 1=bottom
void peano_rec(int ia,int ib,int sa,int sb,int in_corner,int out_corner,struct PeanoData * data)
{
  int sa2 = sa>>1;
  int sb2 = sb>>1;
  int c;
  // Final case
  if (sa == 2 && sb == 4)
  {
    (data->visitor)(data->visitor_arg,ia,ib);
    return;
  }
  // Subdivide
  c = in_corner;
  peano_rec(ia+((c&1)?sa2:0),ib+((c&2)?sb2:0),sa2,sb2,in_corner,out_corner^3,data);
  c = out_corner^3;
  peano_rec(ia+((c&1)?sa2:0),ib+((c&2)?sb2:0),sa2,sb2,in_corner,out_corner,data);
  c = in_corner^3;
  peano_rec(ia+((c&1)?sa2:0),ib+((c&2)?sb2:0),sa2,sb2,in_corner,out_corner,data);
  c = out_corner;
  peano_rec(ia+((c&1)?sa2:0),ib+((c&2)?sb2:0),sa2,sb2,in_corner^3,out_corner,data);
}

void peano_pattern(void * arg,int p,BlockVisitor_2x4_Proc visitor,void * visitor_arg)
{
  struct PeanoData data;
  data.p = p;
  data.visitor = visitor;
  data.visitor_arg = visitor_arg;

  peano_rec(0,0,p>>1,p,1,3,&data);
  peano_rec(p>>1,0,p>>1,p,2,0,&data);
}

void idendifyPattern(BlockPattern_2x4_Proc pattern,void * pattern_arg,char * s,int len)
{
  if (pattern == columns_pattern)
  {
    ColumnsPatternData * data = (ColumnsPatternData *)pattern_arg;
    snprintf(s,len,"COLUMNS(NCOLS=%d)",data->nb);
    return;
  }
  if (pattern == rows_pattern)
  {
    RowsPatternData * data = (RowsPatternData *)pattern_arg;
    snprintf(s,len,"ROWS(NROWS=%d)",data->na);
    return;
  }
  if (pattern == peano_pattern)
  {
    snprintf(s,len,"PEANO");
    return;
  }
  if (pattern == blocks_pattern)
  {
    BlocksPatternData * data = (BlocksPatternData *)pattern_arg;
    snprintf(s,len,"BLOCKS(NROWS=%d,NCOLS=%d,ORDER=%d)",data->na,data->nb,data->order);
    return;
  }

  // Default
  snprintf(s,len,"UNDEF");
}
