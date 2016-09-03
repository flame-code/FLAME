//! @file
//! Blocks visitors (dgemmsy)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 
//! Eric Bainville, Mar 2010


#include <stdlib.h>
#include <string.h>

#include "visitors.h"
#include "utils.h"
#include "gemm_block_c.h"

void Nop_visitor(void * arg,int ia,int ib)
{
}

void Print_visitor(void * arg,int ia,int ib)
{
  printf("%d %d\n",ia,ib);
}

void initCheckPatternData(CheckPatternData * data,int p)
{
  size_t sz = p*p*sizeof(int); // to allocate
  if (data == 0 || p <= 0) return; // Invalid
  data->p = p;
  data->y = (int *)malloc(sz);
  data->invalid = 0;
  memset(data->y,0,sz);
}

void cleanupCheckPatternData(CheckPatternData * data)
{
  int total = 0;
  int value = 0,count,p,i;
  double total_processed = 0;
  double p2;
  if (data == 0) return; // Invalid
  p = data->p;
  p2 = (double)p*p;

  // Display value counts
  printf("Pattern check:\n");
  printf("- P = %d\n",p);
  printf("- Invalid values = %d\n",data->invalid);
  while (total < p*p)
  {
    count = 0;
    for (i=0;i<p*p;i++) if (data->y[i] == value) count++;
    printf("- Values[%d] = %d\n",value,count);
    total += count;
    total_processed += value * count;
    value++;
  }
  printf("- Extra values processed = %.2f%%\n",100.0*(total_processed/p2-1.0));
  free(data->y);
}

void CheckPattern_visitor(void * arg,int ia,int ib)
{
  CheckPatternData * data = (CheckPatternData *)arg;
  int p,da,db,aa,bb;
  if (data == 0) return; // Invalid
  p = data->p;

  if (ia < 0 || ia >= p || (ia&1) != 0
    || ib < 0 || ib >= p || (ib&3) != 0)
  {
    data->invalid++;
    return; // Invalid block
  }

  // Update cells and transposed cells
  for (da=0;da<2;da++) for (db=0;db<4;db++)
  {
    aa = ia + da;
    bb = ib + db;
    data->y[aa+p*bb]++;
    if (aa != bb) data->y[bb+p*aa]++;
  }
}

void initTracePatternData(TracePatternData * data,int p,const char * filename)
{
  int side;
  FILE * f = 0;
  if (data == 0) return; // Invalid
  f = fopen(filename,"wt");
  data->p = p;
  data->f = f;
  data->la = data->lb = -1;
  data->index = 0;
  data->length = 550 / p;
  data->pass = 0;
  if (data->length < 1) data->length = 1;
  if (f == 0) return;

  side = (p+1) * data->length;
  // EPS prologue
  fprintf(f,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(f,"%%%%BoundingBox: 0 0 %d %d\n",side,side);
  fprintf(f,"%%%%EndComments\n");
  fprintf(f,"/l %d def\n",data->length);
  fprintf(f,"/p %d def\n",p);
  fprintf(f,"/r l 2 div def\n");
  fprintf(f,"/lp l p mul def\n");
  fprintf(f,"/link_style { 1 0 0 setrgbcolor 3.1 setlinewidth } bind def\n");
  fprintf(f,"/tr { l mul exch p sub neg l mul } bind def\n");
  // d_block IA IB: draw block (IA is row and IB is column)
  fprintf(f,"/d_block { newpath tr moveto l 4 mul 0 rlineto 0 l -2 mul rlineto l -4 mul 0 rlineto closepath gsave 1 setgray fill grestore 0 setgray 1 setlinewidth stroke } bind def\n");
  fprintf(f,"/d_label { newpath 3 1 roll tr moveto l l -1.5 mul rmoveto show } bind def\n");
  // Top, Bottom, Left, Right links
  /*
  fprintf(f,"/d_blink { gsave tr translate newpath 0 l -2 mul moveto l 4 mul 0 rlineto link_style stroke grestore } bind def\n");
  fprintf(f,"/d_tlink { gsave tr translate newpath 0 0 moveto l 4 mul 0 rlineto link_style stroke grestore } bind def\n");
  */
  fprintf(f,"/d_llink { gsave tr translate newpath l 2 mul l -1 mul moveto l -2 mul 0 rlineto link_style stroke grestore } bind def\n");
  fprintf(f,"/d_rlink { gsave tr translate newpath l 2 mul l -1 mul moveto l 2 mul 0 rlineto link_style stroke grestore } bind def\n");
  fprintf(f,"/d_blink { gsave tr translate newpath l 2 mul l -1 mul moveto 0 l neg rlineto link_style stroke grestore } bind def\n");
  fprintf(f,"/d_tlink { gsave tr translate newpath l 2 mul l -1 mul moveto 0 l rlineto link_style stroke grestore } bind def\n");
  fprintf(f,"/Courier-Bold 10 selectfont\n");
  fprintf(f,"gsave\n");
  fprintf(f,"l 2 div dup translate\n");
  fprintf(f,"newpath 0 0 moveto lp 0 rlineto 0 lp rlineto lp neg 0 rlineto closepath 0.95 setgray fill\n");
}

void cleanupTracePatternData(TracePatternData * data)
{
  FILE * f = 0;
  int i;
  if (data == 0 || data->f == 0) return; // Invalid
  f = data->f;

  // Add a 1x1 grid
  fprintf(f,"0.1 setlinewidth 0 0 1 setrgbcolor newpath\n");
  for (i=0;i<=data->p;i++)
  {
    fprintf(f,"0 %d l mul moveto lp 0 rlineto\n",i);
    fprintf(f,"%d l mul 0 moveto 0 lp rlineto\n",i);
  }
  fprintf(f,"stroke\n");

  fprintf(f,"grestore\n");
  fclose(f);
}

void TracePattern_visitor(void * arg,int ia,int ib)
{
  FILE * f;
  int aa,bb,p;
  TracePatternData * data = (TracePatternData *)arg;
  if (data == 0 || data->f == 0) return; // Invalid
  f = data->f;
  p = data->p;

  if (data->pass == 0)
  {
    fprintf(f,"%d %d d_block\n",ia,ib);

    aa = (ia+p-2)%p; bb = ib;
    if (aa == data->la && bb == data->lb)
    {
      fprintf(f,"%d %d d_blink\n",aa,bb);
      fprintf(f,"%d %d d_tlink\n",ia,ib);
    }
    aa = (ia+2)%p; bb = ib;
    if (aa == data->la && bb == data->lb)
    {
      fprintf(f,"%d %d d_tlink\n",aa,bb);
      fprintf(f,"%d %d d_blink\n",ia,ib);
    }
    aa = ia; bb = (ib+p-4)%p;
    if (aa == data->la && bb == data->lb)
    {
      fprintf(f,"%d %d d_rlink\n",aa,bb);
      fprintf(f,"%d %d d_llink\n",ia,ib);
    }
    aa = ia; bb = (ib+4)%p;
    if (aa == data->la && bb == data->lb)
    {
      fprintf(f,"%d %d d_llink\n",aa,bb);
      fprintf(f,"%d %d d_rlink\n",ia,ib);
    }
  }
  if (data->pass == 1)
  {
    data->index++;
    fprintf(f,"%d %d (%3d) d_label\n",ia,ib,data->index);
  }
  data->la = ia;
  data->lb = ib;
}

void initComputeData(ComputeData * data,int p,int s,const double * pa,const double * pb,double * y)
{
  if (data == 0) return; // Invalid
  data->p = p;
  data->s = s;
  data->pa = pa;
  data->pb = pb;
  data->y = y;
}

void cleanupComputeData(ComputeData * data)
{
  if (data == 0) return; // Invalid
}

void Compute_visitor(void * arg,int ia,int ib)
{
  ComputeData * data = (ComputeData *)arg;
  double *y;
  int p;
  if (data == 0) return; // Invalid

  p = data->p;
  y = data->y + ia + p * ib;
  gemm_block_2x4_c(data->pa+ia*data->s,data->pb+ib*data->s,data->s,y,p);
}

void initTransposeData(TransposeData * data,int p,double * y)
{
  if (data == 0) return; // Invalid
  data->p = p;
  data->y = y;
}

void cleanupTransposeData(TransposeData * data)
{
}

void Transpose_visitor(void * arg,int ia,int ib)
{
  double *y,*ty;
  int p;
  TransposeData * data = (TransposeData *)arg;
  if (data == 0) return; // Invalid

  // Update Y^T
  p = data->p;
  y = data->y + ia + p * ib;
  ty = data->y + ib + p * ia;

  ty[0] = y[0];
  ty[1] = y[p];
  ty[2] = y[2*p];
  ty[3] = y[3*p];
  ty[p] = y[1];
  ty[1+p] = y[1+p];
  ty[2+p] = y[1+2*p];
  ty[3+p] = y[1+3*p];
}
