!> @file
!!  Test of the buffer allocation informations examples of override
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program f_buffer_allocations
  use futile
  implicit none
  !declare a number of allocatable arrays an pointers of different implicit type and ranks
  logical, dimension(:), allocatable ::             l1
  logical, dimension(:,:), allocatable ::           l2
  logical, dimension(:,:,:), allocatable ::         l3
  !logical, dimension(:,:,:,:), allocatable ::       l4
  !logical, dimension(:,:,:,:,:), allocatable ::     l5
  !logical, dimension(:,:,:,:,:,:), allocatable ::   l6
  !logical, dimension(:,:,:,:,:,:,:), allocatable :: l7
  logical(f_byte), dimension(:), allocatable ::             b1
  logical(f_byte), dimension(:,:), allocatable ::           b2
  logical(f_byte), dimension(:,:,:), allocatable ::         b3
  !logical(f_byte), dimension(:,:,:,:), allocatable ::       b4
  !logical(f_byte), dimension(:,:,:,:,:), allocatable ::     b5
  !logical(f_byte), dimension(:,:,:,:,:,:), allocatable ::   b6
  !logical(f_byte), dimension(:,:,:,:,:,:,:), allocatable :: b7
  integer(f_integer), dimension(:), allocatable ::             i1
  integer(f_integer), dimension(:,:), allocatable ::           i2
  integer(f_integer), dimension(:,:,:), allocatable ::         i3
  integer(f_integer), dimension(:,:,:,:), allocatable ::       i4
  !integer(f_integer), dimension(:,:,:,:,:), allocatable ::     i5
  !integer(f_integer), dimension(:,:,:,:,:,:), allocatable ::   i6
  !integer(f_integer), dimension(:,:,:,:,:,:,:), allocatable :: i7
  !integer(f_short), dimension(:), allocatable ::             s1
  !integer(f_short), dimension(:,:), allocatable ::           s2
  !integer(f_short), dimension(:,:,:), allocatable ::         s3
  !integer(f_short), dimension(:,:,:,:), allocatable ::       s4
  !integer(f_short), dimension(:,:,:,:,:), allocatable ::     s5
  !integer(f_short), dimension(:,:,:,:,:,:), allocatable ::   s6
  !integer(f_short), dimension(:,:,:,:,:,:,:), allocatable :: s7
  integer(f_long), dimension(:), allocatable ::             li1
  integer(f_long), dimension(:,:), allocatable ::           li2
  integer(f_long), dimension(:,:,:), allocatable ::         li3
  integer(f_long), dimension(:,:,:,:), allocatable ::       li4
  !integer(f_long), dimension(:,:,:,:,:), allocatable ::     li5
  !integer(f_long), dimension(:,:,:,:,:,:), allocatable ::   li6
  !integer(f_long), dimension(:,:,:,:,:,:,:), allocatable :: li7
  real(f_simple), dimension(:), allocatable ::             r1
  real(f_simple), dimension(:,:), allocatable ::           r2
  real(f_simple), dimension(:,:,:), allocatable ::         r3
  real(f_simple), dimension(:,:,:,:), allocatable ::       r4
  !real(f_simple), dimension(:,:,:,:,:), allocatable ::     r5
  !real(f_simple), dimension(:,:,:,:,:,:), allocatable ::   r6
  !real(f_simple), dimension(:,:,:,:,:,:,:), allocatable :: r7
  real(f_double), dimension(:), allocatable ::             d1
  real(f_double), dimension(:,:), allocatable ::           d2
  real(f_double), dimension(:,:,:), allocatable ::         d3
  real(f_double), dimension(:,:,:,:), allocatable ::       d4
  real(f_double), dimension(:,:,:,:,:), allocatable ::     d5
  real(f_double), dimension(:,:,:,:,:,:), allocatable ::   d6
  real(f_double), dimension(:,:,:,:,:,:,:), allocatable :: d7
  !real(f_quadruple), dimension(:), allocatable ::             q1
  !real(f_quadruple), dimension(:,:), allocatable ::           q2
  !real(f_quadruple), dimension(:,:,:), allocatable ::         q3
  !real(f_quadruple), dimension(:,:,:,:), allocatable ::       q4
  !real(f_quadruple), dimension(:,:,:,:,:), allocatable ::     q5
  !real(f_quadruple), dimension(:,:,:,:,:,:), allocatable ::   q6
  !real(f_quadruple), dimension(:,:,:,:,:,:,:), allocatable :: q7
  !complex(f_simple), dimension(:), allocatable ::             c1
  !complex(f_simple), dimension(:,:), allocatable ::           c2
  !complex(f_simple), dimension(:,:,:), allocatable ::         c3
  !complex(f_simple), dimension(:,:,:,:), allocatable ::       c4
  !complex(f_simple), dimension(:,:,:,:,:), allocatable ::     c5
  !complex(f_simple), dimension(:,:,:,:,:,:), allocatable ::   c6
  !complex(f_simple), dimension(:,:,:,:,:,:,:), allocatable :: c7
  !complex(f_double), dimension(:), allocatable ::             z1
  complex(f_double), dimension(:,:), allocatable ::           z2
  complex(f_double), dimension(:,:,:), allocatable ::         z3
  !complex(f_double), dimension(:,:,:,:), allocatable ::       z4
  !complex(f_double), dimension(:,:,:,:,:), allocatable ::     z5
  !complex(f_double), dimension(:,:,:,:,:,:), allocatable ::   z6
  !complex(f_double), dimension(:,:,:,:,:,:,:), allocatable :: z7
  logical, dimension(:), pointer ::             l1_ptr
  logical, dimension(:,:), pointer ::           l2_ptr
  logical, dimension(:,:,:), pointer ::         l3_ptr
  !logical, dimension(:,:,:,:), pointer ::       l4_ptr
  !logical, dimension(:,:,:,:,:), pointer ::     l5_ptr
  !logical, dimension(:,:,:,:,:,:), pointer ::   l6_ptr
  !logical, dimension(:,:,:,:,:,:,:), pointer :: l7_ptr
  logical(f_byte), dimension(:), pointer ::             b1_ptr
  logical(f_byte), dimension(:,:), pointer ::           b2_ptr
  logical(f_byte), dimension(:,:,:), pointer ::         b3_ptr
  !logical(f_byte), dimension(:,:,:,:), pointer ::       b4_ptr
  !logical(f_byte), dimension(:,:,:,:,:), pointer ::     b5_ptr
  !logical(f_byte), dimension(:,:,:,:,:,:), pointer ::   b6_ptr
  !logical(f_byte), dimension(:,:,:,:,:,:,:), pointer :: b7_ptr
  integer(f_integer), dimension(:), pointer ::             i1_ptr
  integer(f_integer), dimension(:,:), pointer ::           i2_ptr
  integer(f_integer), dimension(:,:,:), pointer ::         i3_ptr
  integer(f_integer), dimension(:,:,:,:), pointer ::       i4_ptr
  !integer(f_integer), dimension(:,:,:,:,:), pointer ::     i5_ptr
  !integer(f_integer), dimension(:,:,:,:,:,:), pointer ::   i6_ptr
  !integer(f_integer), dimension(:,:,:,:,:,:,:), pointer :: i7_ptr
  !integer(f_short), dimension(:), pointer ::             s1_ptr
  !integer(f_short), dimension(:,:), pointer ::           s2_ptr
  !integer(f_short), dimension(:,:,:), pointer ::         s3_ptr
  !integer(f_short), dimension(:,:,:,:), pointer ::       s4_ptr
  !integer(f_short), dimension(:,:,:,:,:), pointer ::     s5_ptr
  !integer(f_short), dimension(:,:,:,:,:,:), pointer ::   s6_ptr
  !integer(f_short), dimension(:,:,:,:,:,:,:), pointer :: s7_ptr
  integer(f_long), dimension(:), pointer ::             li1_ptr
  !integer(f_long), dimension(:,:), pointer ::           li2_ptr
  !integer(f_long), dimension(:,:,:), pointer ::         li3_ptr
  !integer(f_long), dimension(:,:,:,:), pointer ::       li4_ptr
  !integer(f_long), dimension(:,:,:,:,:), pointer ::     li5_ptr
  !integer(f_long), dimension(:,:,:,:,:,:), pointer ::   li6_ptr
  !integer(f_long), dimension(:,:,:,:,:,:,:), pointer :: li7_ptr
  !real(f_simple), dimension(:), pointer ::             r1_ptr
  !real(f_simple), dimension(:,:), pointer ::           r2_ptr
  !real(f_simple), dimension(:,:,:), pointer ::         r3_ptr
  !real(f_simple), dimension(:,:,:,:), pointer ::       r4_ptr
  !real(f_simple), dimension(:,:,:,:,:), pointer ::     r5_ptr
  !real(f_simple), dimension(:,:,:,:,:,:), pointer ::   r6_ptr
  !real(f_simple), dimension(:,:,:,:,:,:,:), pointer :: r7_ptr
  real(f_double), dimension(:), pointer ::             d1_ptr
  real(f_double), dimension(:,:), pointer ::           d2_ptr
  real(f_double), dimension(:,:,:), pointer ::         d3_ptr
  real(f_double), dimension(:,:,:,:), pointer ::       d4_ptr
  real(f_double), dimension(:,:,:,:,:), pointer ::     d5_ptr
  real(f_double), dimension(:,:,:,:,:,:), pointer ::   d6_ptr
  !real(f_double), dimension(:,:,:,:,:,:,:), pointer :: d7_ptr
  !real(f_quadruple), dimension(:), pointer ::             q1_ptr
  !real(f_quadruple), dimension(:,:), pointer ::           q2_ptr
  !real(f_quadruple), dimension(:,:,:), pointer ::         q3_ptr
  !real(f_quadruple), dimension(:,:,:,:), pointer ::       q4_ptr
  !real(f_quadruple), dimension(:,:,:,:,:), pointer ::     q5_ptr
  !real(f_quadruple), dimension(:,:,:,:,:,:), pointer ::   q6_ptr
  !real(f_quadruple), dimension(:,:,:,:,:,:,:), pointer :: q7_ptr
  !complex(f_simple), dimension(:), pointer ::             c1_ptr
  !complex(f_simple), dimension(:,:), pointer ::           c2_ptr
  !complex(f_simple), dimension(:,:,:), pointer ::         c3_ptr
  !complex(f_simple), dimension(:,:,:,:), pointer ::       c4_ptr
  !complex(f_simple), dimension(:,:,:,:,:), pointer ::     c5_ptr
  !complex(f_simple), dimension(:,:,:,:,:,:), pointer ::   c6_ptr
  !complex(f_simple), dimension(:,:,:,:,:,:,:), pointer :: c7_ptr
  !complex(f_double), dimension(:), pointer ::             z1_ptr
  !complex(f_double), dimension(:,:), pointer ::           z2_ptr
  !complex(f_double), dimension(:,:,:), pointer ::         z3_ptr
  !complex(f_double), dimension(:,:,:,:), pointer ::       z4_ptr
  !complex(f_double), dimension(:,:,:,:,:), pointer ::     z5_ptr
  !complex(f_double), dimension(:,:,:,:,:,:), pointer ::   z6_ptr
  !complex(f_double), dimension(:,:,:,:,:,:,:), pointer :: z7_ptr
  !allocation ranges
  integer, dimension(1) :: n1
  integer, dimension(2) :: n2
  integer, dimension(3) :: n3
  integer, dimension(4) :: n4
  integer, dimension(5) :: n5
  integer, dimension(6) :: n6
  integer, dimension(7) :: n7

  !other tests
  real(f_double), dimension(:), pointer :: d1_ptr_exotic

  call f_lib_initialize()

  call yaml_new_document()

  !pick a typical size for allocation
  n1=5
  n2=5
  n3=5
  n4=5
  n5=5
  n6=5
  n7=5

  !here all the allocations with intrinsics, I comment out the ones which are not yet implemented

  l1=f_malloc(n1,id='l1'); call buffer_info(shape(l1),lbound(l1),ubound(l1),kind(l1),'Logical'); call f_free(l1)
  l2=f_malloc(n2,id='l2'); call buffer_info(shape(l2),lbound(l2),ubound(l2),kind(l2),'Logical'); call f_free(l2)
  l3=f_malloc(n3,id='l3'); call buffer_info(shape(l3),lbound(l3),ubound(l3),kind(l3),'Logical'); call f_free(l3)
  !l4=f_malloc(n4,id='l4'); call buffer_info(shape(l4),lbound(l4),ubound(l4),kind(l4),'Logical'); call f_free(l4)
  !l5=f_malloc(n5,id='l5'); call buffer_info(shape(l5),lbound(l5),ubound(l5),kind(l5),'Logical'); call f_free(l5)
  !l6=f_malloc(n6,id='l6'); call buffer_info(shape(l6),lbound(l6),ubound(l6),kind(l6),'Logical'); call f_free(l6)
  !l7=f_malloc(n7,id='l7'); call buffer_info(shape(l7),lbound(l7),ubound(l7),kind(l7),'Logical'); call f_free(l7)

  b1=f_malloc(n1,id='b1'); call buffer_info(shape(b1),lbound(b1),ubound(b1),kind(b1),'Byte'); call f_free(b1)
  b2=f_malloc(n2,id='b2'); call buffer_info(shape(b2),lbound(b2),ubound(b2),kind(b2),'Byte'); call f_free(b2)
  b3=f_malloc(n3,id='b3'); call buffer_info(shape(b3),lbound(b3),ubound(b3),kind(b3),'Byte'); call f_free(b3)
  !b4=f_malloc(n4,id='b4'); call buffer_info(shape(b4),lbound(b4),ubound(b4),kind(b4),'Byte'); call f_free(b4)
  !b5=f_malloc(n5,id='b5'); call buffer_info(shape(b5),lbound(b5),ubound(b5),kind(b5),'Byte'); call f_free(b5)
  !b6=f_malloc(n6,id='b6'); call buffer_info(shape(b6),lbound(b6),ubound(b6),kind(b6),'Byte'); call f_free(b6)
  !b7=f_malloc(n7,id='b7'); call buffer_info(shape(b7),lbound(b7),ubound(b7),kind(b7),'Byte'); call f_free(b7)

  !s1=f_malloc(n1,id='s1'); call buffer_info(shape(s1),lbound(s1),ubound(s1),kind(s1),'Short'); call f_free(s1)
  !s2=f_malloc(n2,id='s2'); call buffer_info(shape(s2),lbound(s2),ubound(s2),kind(s2),'Short'); call f_free(s2)
  !s3=f_malloc(n3,id='s3'); call buffer_info(shape(s3),lbound(s3),ubound(s3),kind(s3),'Short'); call f_free(s3)
  !s4=f_malloc(n4,id='s4'); call buffer_info(shape(s4),lbound(s4),ubound(s4),kind(s4),'Short'); call f_free(s4)
  !s5=f_malloc(n5,id='s5'); call buffer_info(shape(s5),lbound(s5),ubound(s5),kind(s5),'Short'); call f_free(s5)
  !s6=f_malloc(n6,id='s6'); call buffer_info(shape(s6),lbound(s6),ubound(s6),kind(s6),'Short'); call f_free(s6)
  !s7=f_malloc(n7,id='s7'); call buffer_info(shape(s7),lbound(s7),ubound(s7),kind(s7),'Short'); call f_free(s7)

  i1=f_malloc(n1,id='i1'); call buffer_info(shape(i1),lbound(i1),ubound(i1),kind(i1),'Int'); call f_free(i1)
  i2=f_malloc(n2,id='i2'); call buffer_info(shape(i2),lbound(i2),ubound(i2),kind(i2),'Int'); call f_free(i2)
  i3=f_malloc(n3,id='i3'); call buffer_info(shape(i3),lbound(i3),ubound(i3),kind(i3),'Int'); call f_free(i3)
  i4=f_malloc(n4,id='i4'); call buffer_info(shape(i4),lbound(i4),ubound(i4),kind(i4),'Int'); call f_free(i4)
  !i5=f_malloc(n5,id='i5'); call buffer_info(shape(i5),lbound(i5),ubound(i5),kind(i5),'Int'); call f_free(i5)
  !i6=f_malloc(n6,id='i6'); call buffer_info(shape(i6),lbound(i6),ubound(i6),kind(i6),'Int'); call f_free(i6)
  !i7=f_malloc(n7,id='i7'); call buffer_info(shape(i7),lbound(i7),ubound(i7),kind(i7),'Int'); call f_free(i7)

  li1=f_malloc(n1,id='li1'); call buffer_info(shape(li1),lbound(li1),ubound(li1),kind(li1),'Long'); call f_free(li1)
  li2=f_malloc(n2,id='li2'); call buffer_info(shape(li2),lbound(li2),ubound(li2),kind(li2),'Long'); call f_free(li2)
  li3=f_malloc(n3,id='li3'); call buffer_info(shape(li3),lbound(li3),ubound(li3),kind(li3),'Long'); call f_free(li3)
  li4=f_malloc(n4,id='li4'); call buffer_info(shape(li4),lbound(li4),ubound(li4),kind(li4),'Long'); call f_free(li4)
  !li5=f_malloc(n5,id='li5'); call buffer_info(shape(li5),lbound(li5),ubound(li5),kind(li5),'Long'); call f_free(li5)
  !li6=f_malloc(n6,id='li6'); call buffer_info(shape(li6),lbound(li6),ubound(li6),kind(li6),'Long'); call f_free(li6)
  !li7=f_malloc(n7,id='li7'); call buffer_info(shape(li7),lbound(li7),ubound(li7),kind(li7),'Long'); call f_free(li7)

  r1=f_malloc(n1,id='r1'); call buffer_info(shape(r1),lbound(r1),ubound(r1),kind(r1),'Float');
  call f_zero(r1)
  call detect_rnan(r1,shape(r1))
  call f_free(r1)
  r2=f_malloc(n2,id='r2'); call buffer_info(shape(r2),lbound(r2),ubound(r2),kind(r2),'Float'); call f_free(r2)
  r3=f_malloc(n3,id='r3'); call buffer_info(shape(r3),lbound(r3),ubound(r3),kind(r3),'Float'); call f_free(r3)
  r4=f_malloc(n4,id='r4'); call buffer_info(shape(r4),lbound(r4),ubound(r4),kind(r4),'Float'); call f_free(r4)
  !r5=f_malloc(n5,id='r5'); call buffer_info(shape(r5),lbound(r5),ubound(r5),kind(r5),'Float'); call f_free(r5)
  !r6=f_malloc(n6,id='r6'); call buffer_info(shape(r6),lbound(r6),ubound(r6),kind(r6),'Float'); call f_free(r6)
  !r7=f_malloc(n7,id='r7'); call buffer_info(shape(r7),lbound(r7),ubound(r7),kind(r7),'Float'); call f_free(r7)

  d1=f_malloc(n1,id='d1'); call buffer_info(shape(d1),lbound(d1),ubound(d1),kind(d1),'Double')
  call f_zero(d1)
  call detect_dnan(d1,shape(d1))
  call f_free(d1)
  d2=f_malloc(n2,id='d2'); call buffer_info(shape(d2),lbound(d2),ubound(d2),kind(d2),'Double')
  call f_zero(d2)
  call detect_dnan(d2,shape(d2))
  call f_free(d2)
  d3=f_malloc(n3,id='d3'); call buffer_info(shape(d3),lbound(d3),ubound(d3),kind(d3),'Double')
  call f_zero(d3)
  call detect_dnan(d3,shape(d3))
  call f_free(d3)
  d4=f_malloc(n4,id='d4'); call buffer_info(shape(d4),lbound(d4),ubound(d4),kind(d4),'Double');
  call f_zero(d4)
  call detect_dnan(d4,shape(d4))
  call f_free(d4)
  d5=f_malloc(n5,id='d5'); call buffer_info(shape(d5),lbound(d5),ubound(d5),kind(d5),'Double')
  call f_zero(d5)
  call detect_dnan(d5,shape(d5))
  call f_free(d5)
  d6=f_malloc(n6,id='d6'); call buffer_info(shape(d6),lbound(d6),ubound(d6),kind(d6),'Double')
  call f_zero(d6)
  call detect_dnan(d6,shape(d6))
  call f_free(d6)
  d7=f_malloc(n7,id='d7'); call buffer_info(shape(d7),lbound(d7),ubound(d7),kind(d7),'Double')
  call f_zero(d7)
  call detect_dnan(d7,shape(d7))
 call f_free(d7)

  !q1=f_malloc(n1,id='q1'); call buffer_info(shape(q1),lbound(q1),ubound(q1),kind(q1),'Quadruple'); call f_free(q1)
  !q2=f_malloc(n2,id='q2'); call buffer_info(shape(q2),lbound(q2),ubound(q2),kind(q2),'Quadruple'); call f_free(q2)
  !q3=f_malloc(n3,id='q3'); call buffer_info(shape(q3),lbound(q3),ubound(q3),kind(q3),'Quadruple'); call f_free(q3)
  !q4=f_malloc(n4,id='q4'); call buffer_info(shape(q4),lbound(q4),ubound(q4),kind(q4),'Quadruple'); call f_free(q4)
  !q5=f_malloc(n5,id='q5'); call buffer_info(shape(q5),lbound(q5),ubound(q5),kind(q5),'Quadruple'); call f_free(q5)
  !q6=f_malloc(n6,id='q6'); call buffer_info(shape(q6),lbound(q6),ubound(q6),kind(q6),'Quadruple'); call f_free(q6)
  !q7=f_malloc(n7,id='q7'); call buffer_info(shape(q7),lbound(q7),ubound(q7),kind(q7),'Quadruple'); call f_free(q7)

  !c1=f_malloc(n1,id='c1'); call buffer_info(shape(c1),lbound(c1),ubound(c1),kind(c1),'Complex'); call f_free(c1)
  !c2=f_malloc(n2,id='c2'); call buffer_info(shape(c2),lbound(c2),ubound(c2),kind(c2),'Complex'); call f_free(c2)
  !c3=f_malloc(n3,id='c3'); call buffer_info(shape(c3),lbound(c3),ubound(c3),kind(c3),'Complex'); call f_free(c3)
  !c4=f_malloc(n4,id='c4'); call buffer_info(shape(c4),lbound(c4),ubound(c4),kind(c4),'Complex'); call f_free(c4)
  !c5=f_malloc(n5,id='c5'); call buffer_info(shape(c5),lbound(c5),ubound(c5),kind(c5),'Complex'); call f_free(c5)
  !c6=f_malloc(n6,id='c6'); call buffer_info(shape(c6),lbound(c6),ubound(c6),kind(c6),'Complex'); call f_free(c6)
  !c7=f_malloc(n7,id='c7'); call buffer_info(shape(c7),lbound(c7),ubound(c7),kind(c7),'Complex'); call f_free(c7)

  !z1=f_malloc(n1,id='z1'); call buffer_info(shape(z1),lbound(z1),ubound(z1),kind(z1),'Double Complex'); call f_free(z1)
  z2=f_malloc(n2,id='z2'); call buffer_info(shape(z2),lbound(z2),ubound(z2),kind(z2),'Double Complex'); call f_free(z2)
  z3=f_malloc(n3,id='z3'); call buffer_info(shape(z3),lbound(z3),ubound(z3),kind(z3),'Double Complex'); call f_free(z3)
  !z4=f_malloc(n4,id='z4'); call buffer_info(shape(z4),lbound(z4),ubound(z4),kind(z4),'Double Complex'); call f_free(z4)
  !z5=f_malloc(n5,id='z5'); call buffer_info(shape(z5),lbound(z5),ubound(z5),kind(z5),'Double Complex'); call f_free(z5)
  !z6=f_malloc(n6,id='z6'); call buffer_info(shape(z6),lbound(z6),ubound(z6),kind(z6),'Double Complex'); call f_free(z6)
  !z7=f_malloc(n7,id='z7'); call buffer_info(shape(z7),lbound(z7),ubound(z7),kind(z7),'Double Complex'); call f_free(z7)

  !pointers

  l1_ptr=f_malloc_ptr(n1,id='l1_ptr')
  call buffer_info(shape(l1_ptr),lbound(l1_ptr),ubound(l1_ptr),kind(l1_ptr),'Logical_ptr')
  call f_free_ptr(l1_ptr)
  l2_ptr=f_malloc_ptr(n2,id='l2_ptr')
  call buffer_info(shape(l2_ptr),lbound(l2_ptr),ubound(l2_ptr),kind(l2_ptr),'Logical_ptr')
  call f_free_ptr(l2_ptr)
  l3_ptr=f_malloc_ptr(n3,id='l3_ptr')
  call buffer_info(shape(l3_ptr),lbound(l3_ptr),ubound(l3_ptr),kind(l3_ptr),'Logical_ptr')
  call f_free_ptr(l3_ptr)
  !l4_ptr=f_malloc_ptr(n4,id='l4_ptr')
  !call buffer_info(shape(l4_ptr),lbound(l4_ptr),ubound(l4_ptr),kind(l4_ptr),'Logical_ptr')
  !call f_free_ptr(l4_ptr)
  !l5_ptr=f_malloc_ptr(n5,id='l5_ptr')
  !call buffer_info(shape(l5_ptr),lbound(l5_ptr),ubound(l5_ptr),kind(l5_ptr),'Logical_ptr')
  !call f_free_ptr(l5_ptr)
  !l6_ptr=f_malloc_ptr(n6,id='l6_ptr')
  !call buffer_info(shape(l6_ptr),lbound(l6_ptr),ubound(l6_ptr),kind(l6_ptr),'Logical_ptr')
  !call f_free_ptr(l6_ptr)
  !l7_ptr=f_malloc_ptr(n7,id='l7_ptr')
  !call buffer_info(shape(l7_ptr),lbound(l7_ptr),ubound(l7_ptr),kind(l7_ptr),'Logical_ptr')
  !call f_free_ptr(l7_ptr)

  b1_ptr=f_malloc_ptr(n1,id='b1_ptr')
  call buffer_info(shape(b1_ptr),lbound(b1_ptr),ubound(b1_ptr),kind(b1_ptr),'Byte_ptr')
  call f_free_ptr(b1_ptr)
  b2_ptr=f_malloc_ptr(n2,id='b2_ptr')
  call buffer_info(shape(b2_ptr),lbound(b2_ptr),ubound(b2_ptr),kind(b2_ptr),'Byte_ptr')
  call f_free_ptr(b2_ptr)
  b3_ptr=f_malloc_ptr(n3,id='b3_ptr')
  call buffer_info(shape(b3_ptr),lbound(b3_ptr),ubound(b3_ptr),kind(b3_ptr),'Byte_ptr')
  call f_free_ptr(b3_ptr)
  !b4_ptr=f_malloc_ptr(n4,id='b4_ptr')
  !call buffer_info(shape(b4_ptr),lbound(b4_ptr),ubound(b4_ptr),kind(b4_ptr),'Byte_ptr')
  !call f_free_ptr(b4_ptr)
  !b5_ptr=f_malloc_ptr(n5,id='b5_ptr')
  !call buffer_info(shape(b5_ptr),lbound(b5_ptr),ubound(b5_ptr),kind(b5_ptr),'Byte_ptr')
  !call f_free_ptr(b5_ptr)
  !b6_ptr=f_malloc_ptr(n6,id='b6_ptr')
  !call buffer_info(shape(b6_ptr),lbound(b6_ptr),ubound(b6_ptr),kind(b6_ptr),'Byte_ptr')
  !call f_free_ptr(b6_ptr)
  !b7_ptr=f_malloc_ptr(n7,id='b7_ptr')
  !call buffer_info(shape(b7_ptr),lbound(b7_ptr),ubound(b7_ptr),kind(b7_ptr),'Byte_ptr')
  !call f_free_ptr(b7_ptr)

  !s1_ptr=f_malloc_ptr(n1,id='s1_ptr')
  !call buffer_info(shape(s1_ptr),lbound(s1_ptr),ubound(s1_ptr),kind(s1_ptr),'Short_ptr')
  !call f_free_ptr(s1_ptr)
  !s2_ptr=f_malloc_ptr(n2,id='s2_ptr')
  !call buffer_info(shape(s2_ptr),lbound(s2_ptr),ubound(s2_ptr),kind(s2_ptr),'Short_ptr')
  !call f_free_ptr(s2_ptr)
  !s3_ptr=f_malloc_ptr(n3,id='s3_ptr')
  !call buffer_info(shape(s3_ptr),lbound(s3_ptr),ubound(s3_ptr),kind(s3_ptr),'Short_ptr')
  !call f_free_ptr(s3_ptr)
  !s4_ptr=f_malloc_ptr(n4,id='s4_ptr')
  !call buffer_info(shape(s4_ptr),lbound(s4_ptr),ubound(s4_ptr),kind(s4_ptr),'Short_ptr')
  !call f_free_ptr(s4_ptr)
  !s5_ptr=f_malloc_ptr(n5,id='s5_ptr')
  !call buffer_info(shape(s5_ptr),lbound(s5_ptr),ubound(s5_ptr),kind(s5_ptr),'Short_ptr')
  !call f_free_ptr(s5_ptr)
  !s6_ptr=f_malloc_ptr(n6,id='s6_ptr')
  !call buffer_info(shape(s6_ptr),lbound(s6_ptr),ubound(s6_ptr),kind(s6_ptr),'Short_ptr')
  !call f_free_ptr(s6_ptr)
  !s7_ptr=f_malloc_ptr(n7,id='s7_ptr')
  !call buffer_info(shape(s7_ptr),lbound(s7_ptr),ubound(s7_ptr),kind(s7_ptr),'Short_ptr')
  !call f_free_ptr(s7_ptr)

  i1_ptr=f_malloc_ptr(n1,id='i1_ptr')
  call buffer_info(shape(i1_ptr),lbound(i1_ptr),ubound(i1_ptr),kind(i1_ptr),'Int_ptr')
  call f_free_ptr(i1_ptr)
  i2_ptr=f_malloc_ptr(n2,id='i2_ptr')
  call buffer_info(shape(i2_ptr),lbound(i2_ptr),ubound(i2_ptr),kind(i2_ptr),'Int_ptr')
  call f_free_ptr(i2_ptr)
  i3_ptr=f_malloc_ptr(n3,id='i3_ptr')
  call buffer_info(shape(i3_ptr),lbound(i3_ptr),ubound(i3_ptr),kind(i3_ptr),'Int_ptr')
  call f_free_ptr(i3_ptr)
  i4_ptr=f_malloc_ptr(n4,id='i4_ptr')
  call buffer_info(shape(i4_ptr),lbound(i4_ptr),ubound(i4_ptr),kind(i4_ptr),'Int_ptr')
  call f_free_ptr(i4_ptr)
  !i5_ptr=f_malloc_ptr(n5,id='i5_ptr')
  !call buffer_info(shape(i5_ptr),lbound(i5_ptr),ubound(i5_ptr),kind(i5_ptr),'Int_ptr')
  !call f_free_ptr(i5_ptr)
  !i6_ptr=f_malloc_ptr(n6,id='i6_ptr')
  !call buffer_info(shape(i6_ptr),lbound(i6_ptr),ubound(i6_ptr),kind(i6_ptr),'Int_ptr')
  !call f_free_ptr(i6_ptr)
  !i7_ptr=f_malloc_ptr(n7,id='i7_ptr')
  !call buffer_info(shape(i7_ptr),lbound(i7_ptr),ubound(i7_ptr),kind(i7_ptr),'Int_ptr')
  !call f_free_ptr(i7_ptr)

  li1_ptr=f_malloc_ptr(n1,id='li1_ptr')
  call buffer_info(shape(li1_ptr),lbound(li1_ptr),ubound(li1_ptr),kind(li1_ptr),'Long_ptr')
  call f_free_ptr(li1_ptr)
  !li2_ptr=f_malloc_ptr(n2,id='li2_ptr')
  !call buffer_info(shape(li2_ptr),lbound(li2_ptr),ubound(li2_ptr),kind(li2_ptr),'Long_ptr')
  !call f_free_ptr(li2_ptr)
  !li3_ptr=f_malloc_ptr(n3,id='li3_ptr')
  !call buffer_info(shape(li3_ptr),lbound(li3_ptr),ubound(li3_ptr),kind(li3_ptr),'Long_ptr')
  !call f_free_ptr(li3_ptr)
  !li4_ptr=f_malloc_ptr(n4,id='li4_ptr')
  !call buffer_info(shape(li4_ptr),lbound(li4_ptr),ubound(li4_ptr),kind(li4_ptr),'Long_ptr')
  !call f_free_ptr(li4_ptr)
  !li5_ptr=f_malloc_ptr(n5,id='li5_ptr')
  !call buffer_info(shape(li5_ptr),lbound(li5_ptr),ubound(li5_ptr),kind(li5_ptr),'Long_ptr')
  !call f_free_ptr(li5_ptr)
  !li6_ptr=f_malloc_ptr(n6,id='li6_ptr')
  !call buffer_info(shape(li6_ptr),lbound(li6_ptr),ubound(li6_ptr),kind(li6_ptr),'Long_ptr')
  !call f_free_ptr(li6_ptr)
  !li7_ptr=f_malloc_ptr(n7,id='li7_ptr')
  !call buffer_info(shape(li7_ptr),lbound(li7_ptr),ubound(li7_ptr),kind(li7_ptr),'Long_ptr')
  !call f_free_ptr(li7_ptr)

  !r1_ptr=f_malloc_ptr(n1,id='r1_ptr')
  !call buffer_info(shape(r1_ptr),lbound(r1_ptr),ubound(r1_ptr),kind(r1_ptr),'Float_ptr')
  !call f_free_ptr(r1_ptr)
  !r2_ptr=f_malloc_ptr(n2,id='r2_ptr')
  !call buffer_info(shape(r2_ptr),lbound(r2_ptr),ubound(r2_ptr),kind(r2_ptr),'Float_ptr')
  !call f_free_ptr(r2_ptr)
  !r3_ptr=f_malloc_ptr(n3,id='r3_ptr')
  !call buffer_info(shape(r3_ptr),lbound(r3_ptr),ubound(r3_ptr),kind(r3_ptr),'Float_ptr')
  !call f_free_ptr(r3_ptr)
  !r4_ptr=f_malloc_ptr(n4,id='r4_ptr')
  !call buffer_info(shape(r4_ptr),lbound(r4_ptr),ubound(r4_ptr),kind(r4_ptr),'Float_ptr')
  !call f_free_ptr(r4_ptr)
  !r5_ptr=f_malloc_ptr(n5,id='r5_ptr')
  !call buffer_info(shape(r5_ptr),lbound(r5_ptr),ubound(r5_ptr),kind(r5_ptr),'Float_ptr')
  !call f_free_ptr(r5_ptr)
  !r6_ptr=f_malloc_ptr(n6,id='r6_ptr')
  !call buffer_info(shape(r6_ptr),lbound(r6_ptr),ubound(r6_ptr),kind(r6_ptr),'Float_ptr')
  !call f_free_ptr(r6_ptr)
  !r7_ptr=f_malloc_ptr(n7,id='r7_ptr')
  !call buffer_info(shape(r7_ptr),lbound(r7_ptr),ubound(r7_ptr),kind(r7_ptr),'Float_ptr')
  !call f_free_ptr(r7_ptr)

  d1_ptr=f_malloc_ptr(n1,id='d1_ptr')
  call buffer_info(shape(d1_ptr),lbound(d1_ptr),ubound(d1_ptr),kind(d1_ptr),'Double_ptr')
  !here we might detect some NaN
  call f_zero(d1_ptr)
  call detect_dnan(d1_ptr,shape(d1_ptr))
  call f_free_ptr(d1_ptr)

  d1_ptr=f_malloc_ptr(n1,id='d1_ptr',info="{Type: SHARED}")
  call buffer_info(shape(d1_ptr),lbound(d1_ptr),ubound(d1_ptr),kind(d1_ptr),'Double_ptr')
  !here we might detect some NaN
  call f_zero(d1_ptr)
  call detect_dnan(d1_ptr,shape(d1_ptr))
  call f_free_ptr(d1_ptr)

  d2_ptr=f_malloc_ptr(n2,id='d2_ptr',info='{alignment: 32}')
  call buffer_info(shape(d2_ptr),lbound(d2_ptr),ubound(d2_ptr),kind(d2_ptr),'Double_ptr')
  !here we might detect some NaN
  call f_zero(d2_ptr)
  call detect_dnan(d2_ptr,shape(d2_ptr))
  call f_free_ptr(d2_ptr)
  d3_ptr=f_malloc_ptr(n3,id='d3_ptr')
  call buffer_info(shape(d3_ptr),lbound(d3_ptr),ubound(d3_ptr),kind(d3_ptr),'Double_ptr')
  !here we might detect some NaN
  call f_zero(d3_ptr)
  call detect_dnan(d3_ptr,shape(d3_ptr))
  call f_free_ptr(d3_ptr)
  d4_ptr=f_malloc_ptr(n4,id='d4_ptr')
  call buffer_info(shape(d4_ptr),lbound(d4_ptr),ubound(d4_ptr),kind(d4_ptr),'Double_ptr')
  call f_zero(d4_ptr)
  call detect_dnan(d4_ptr,shape(d4_ptr))
  call f_free_ptr(d4_ptr)
  d5_ptr=f_malloc_ptr(n5,id='d5_ptr')
  call buffer_info(shape(d5_ptr),lbound(d5_ptr),ubound(d5_ptr),kind(d5_ptr),'Double_ptr')
  call f_zero(d5_ptr)
  call detect_dnan(d5_ptr,shape(d5_ptr))
  call f_free_ptr(d5_ptr)
  d6_ptr=f_malloc_ptr(n6,id='d6_ptr')
  call buffer_info(shape(d6_ptr),lbound(d6_ptr),ubound(d6_ptr),kind(d6_ptr),'Double_ptr')
  call f_zero(d6_ptr)
  call detect_dnan(d6_ptr,shape(d6_ptr))
  call f_free_ptr(d6_ptr)
  !d7_ptr=f_malloc_ptr(n7,id='d7_ptr')
  !call buffer_info(shape(d7_ptr),lbound(d7_ptr),ubound(d7_ptr),kind(d7_ptr),'Double_ptr')
  !call f_free_ptr(d7_ptr)

  !q1_ptr=f_malloc_ptr(n1,id='q1_ptr')
  !call buffer_info(shape(q1_ptr),lbound(q1_ptr),ubound(q1_ptr),kind(q1_ptr),'Quadruple_ptr')
  !call f_free_ptr(q1_ptr)
  !q2_ptr=f_malloc_ptr(n2,id='q2_ptr')
  !call buffer_info(shape(q2_ptr),lbound(q2_ptr),ubound(q2_ptr),kind(q2_ptr),'Quadruple_ptr')
  !call f_free_ptr(q2_ptr)
  !q3_ptr=f_malloc_ptr(n3,id='q3_ptr')
  !call buffer_info(shape(q3_ptr),lbound(q3_ptr),ubound(q3_ptr),kind(q3_ptr),'Quadruple_ptr')
  !call f_free_ptr(q3_ptr)
  !q4_ptr=f_malloc_ptr(n4,id='q4_ptr')
  !call buffer_info(shape(q4_ptr),lbound(q4_ptr),ubound(q4_ptr),kind(q4_ptr),'Quadruple_ptr')
  !call f_free_ptr(q4_ptr)
  !q5_ptr=f_malloc_ptr(n5,id='q5_ptr')
  !call buffer_info(shape(q5_ptr),lbound(q5_ptr),ubound(q5_ptr),kind(q5_ptr),'Quadruple_ptr')
  !call f_free_ptr(q5_ptr)
  !q6_ptr=f_malloc_ptr(n6,id='q6_ptr')
  !call buffer_info(shape(q6_ptr),lbound(q6_ptr),ubound(q6_ptr),kind(q6_ptr),'Quadruple_ptr')
  !call f_free_ptr(q6_ptr)
  !q7_ptr=f_malloc_ptr(n7,id='q7_ptr')
  !call buffer_info(shape(q7_ptr),lbound(q7_ptr),ubound(q7_ptr),kind(q7_ptr),'Quadruple_ptr')
  !call f_free_ptr(q7_ptr)

  !c1_ptr=f_malloc_ptr(n1,id='c1_ptr')
  !call buffer_info(shape(c1_ptr),lbound(c1_ptr),ubound(c1_ptr),kind(c1_ptr),'Complex_ptr')
  !call f_free_ptr(c1_ptr)
  !c2_ptr=f_malloc_ptr(n2,id='c2_ptr')
  !call buffer_info(shape(c2_ptr),lbound(c2_ptr),ubound(c2_ptr),kind(c2_ptr),'Complex_ptr')
  !call f_free_ptr(c2_ptr)
  !c3_ptr=f_malloc_ptr(n3,id='c3_ptr')
  !call buffer_info(shape(c3_ptr),lbound(c3_ptr),ubound(c3_ptr),kind(c3_ptr),'Complex_ptr')
  !call f_free_ptr(c3_ptr)
  !c4_ptr=f_malloc_ptr(n4,id='c4_ptr')
  !call buffer_info(shape(c4_ptr),lbound(c4_ptr),ubound(c4_ptr),kind(c4_ptr),'Complex_ptr')
  !call f_free_ptr(c4_ptr)
  !c5_ptr=f_malloc_ptr(n5,id='c5_ptr')
  !call buffer_info(shape(c5_ptr),lbound(c5_ptr),ubound(c5_ptr),kind(c5_ptr),'Complex_ptr')
  !call f_free_ptr(c5_ptr)
  !c6_ptr=f_malloc_ptr(n6,id='c6_ptr')
  !call buffer_info(shape(c6_ptr),lbound(c6_ptr),ubound(c6_ptr),kind(c6_ptr),'Complex_ptr')
  !call f_free_ptr(c6_ptr)
  !c7_ptr=f_malloc_ptr(n7,id='c7_ptr')
  !call buffer_info(shape(c7_ptr),lbound(c7_ptr),ubound(c7_ptr),kind(c7_ptr),'Complex_ptr')
  !call f_free_ptr(c7_ptr)

  !z1_ptr=f_malloc_ptr(n1,id='z1_ptr')
  !call buffer_info(shape(z1_ptr),lbound(z1_ptr),ubound(z1_ptr),kind(z1_ptr),'Double Complex_ptr')
  !call f_free_ptr(z1_ptr)
  !z2_ptr=f_malloc_ptr(n2,id='z2_ptr')
  !call buffer_info(shape(z2_ptr),lbound(z2_ptr),ubound(z2_ptr),kind(z2_ptr),'Double Complex_ptr')
  !call f_free_ptr(z2_ptr)
  !z3_ptr=f_malloc_ptr(n3,id='z3_ptr')
  !call buffer_info(shape(z3_ptr),lbound(z3_ptr),ubound(z3_ptr),kind(z3_ptr),'Double Complex_ptr')
  !call f_free_ptr(z3_ptr)
  !z4_ptr=f_malloc_ptr(n4,id='z4_ptr')
  !call buffer_info(shape(z4_ptr),lbound(z4_ptr),ubound(z4_ptr),kind(z4_ptr),'Double Complex_ptr')
  !call f_free_ptr(z4_ptr)
  !z5_ptr=f_malloc_ptr(n5,id='z5_ptr')
  !call buffer_info(shape(z5_ptr),lbound(z5_ptr),ubound(z5_ptr),kind(z5_ptr),'Double Complex_ptr')
  !call f_free_ptr(z5_ptr)
  !z6_ptr=f_malloc_ptr(n6,id='z6_ptr')
  !call buffer_info(shape(z6_ptr),lbound(z6_ptr),ubound(z6_ptr),kind(z6_ptr),'Double Complex_ptr')
  !call f_free_ptr(z6_ptr)
  !z7_ptr=f_malloc_ptr(n7,id='z7_ptr')
  !call buffer_info(shape(z7_ptr),lbound(z7_ptr),ubound(z7_ptr),kind(z7_ptr),'Double Complex_ptr')
  !call f_free_ptr(z7_ptr)


  !allocate the pointer in the normal case but with zero size
  d1_ptr=f_malloc_ptr(-1.to.-2,id='d1_ptr')


  !now insert some extra information in the buffer
  !allocate the pointer in the normal case
  d1_ptr_exotic=f_malloc_ptr(-1.to.34,id='d1_ptr_exotic',&
       info='{Type: SHARED}')

  call yaml_mapping_open('Status of the database',advance='no')
  call yaml_comment('Database status',hfill='~')
  call f_malloc_dump_status()
  call yaml_mapping_close(advance='no')
  call yaml_comment('',hfill='~')
  call f_free_ptr(d1_ptr)
  call f_free_ptr(d1_ptr_exotic)
  call f_lib_finalize()

  contains

    subroutine buffer_info(shp,lbnd,ubnd,sizeof,typeb)
      implicit none
      integer, intent(in) :: sizeof
      character(len=*), intent(in) :: typeb
      integer, dimension(:), intent(in) :: shp,lbnd,ubnd

      !call yaml_mapping_open(typeb//' buffer,kind '//trim(yaml_toa(sizeof)))
      call yaml_mapping_open(typeb//' buffer, kind '+sizeof+', rank'+yaml_toa(size(shp)))
        call yaml_map('Shape',shp)
        call yaml_map('Lbound',lbnd)
        call yaml_map('Ubound',ubnd)
      call yaml_mapping_close()
    end subroutine buffer_info

    subroutine detect_dnan(d1_ptr,shp)
      implicit none
      real(f_double), dimension(*), intent(in) :: d1_ptr
      integer, dimension(:), intent(in) :: shp
      !local variable
      integer(f_long) ::  isize

      isize=product(int(shp,f_long))

      call yaml_map('Initial value',d1_ptr(1))
      call yaml_map('Official Final value',d1_ptr(isize))
      if (f_nan_pad_size> 0) call yaml_map('Next value',d1_ptr(isize+1))

    end subroutine detect_dnan

    subroutine detect_rnan(d1_ptr,shp)
      implicit none
      real(f_simple), dimension(*), intent(in) :: d1_ptr
      integer, dimension(:), intent(in) :: shp
      !local variable
      integer(f_long) ::  isize

      isize=product(int(shp,f_long))

      call yaml_map('Initial value',d1_ptr(1))
      call yaml_map('Official Final value',d1_ptr(isize))
      if (f_nan_pad_size> 0) call yaml_map('Next value',d1_ptr(isize+1))

    end subroutine detect_rnan


    !> no way, seems impossible to modify the behaviour of an allocatable array
    subroutine allocate_test(n,origin)
      use iso_c_binding
      use smpi_shared
      implicit none
      integer, intent(in) :: n
      real(f_double), dimension(:), allocatable :: origin
      real(f_double), dimension(:), allocatable, target :: array

      call yaml_map('Original address',f_loc(origin))
      call yaml_map('Original shape',shape(origin))
      call yaml_map('Original values',origin)

      !allocate(array(n),source=origin)

      !call move_alloc(origin,array)

      call yaml_map('Allocated address',f_loc(array))
      call yaml_map('Allocated shape',shape(array))
      call yaml_map('Allocated values',array)

    end subroutine allocate_test

end program f_buffer_allocations
