!> @file
!!  Tool routine for FFT
!! @author
!!    Copyright (C) 2008-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Build the possible number to perform an FFT
!! using radix of 3,4,5,6,8
program fft_factors

   implicit none

   !Maximal combination of factorisation with m_radix numbers
   integer, parameter :: n_max=24000
   !Number of available radix
   integer, parameter :: n_radix_max = 8
   !Number of radix for the decomposition: m_radix
   integer, parameter :: n_decomp_max = 7

   integer, dimension(:,:), allocatable :: factor,selected
   integer, dimension(:), allocatable :: results,indexes
   integer, dimension(n_radix_max) :: radix
   integer, dimension(1) :: i_loc
   integer :: n_element,reference
   integer :: i,j,k,rad,ncount,end1,i_swap,j_swap

   !Allocations
   allocate(factor  (n_decomp_max,n_max))
   allocate(selected(n_decomp_max,n_max))
   allocate(results(n_max))

   !Possible radix
   radix(1)=1
   radix(2)=2
   radix(3)=3
   radix(4)=4
   radix(5)=5
   radix(6)=7
   radix(7)=8
   radix(n_radix_max)=6

   print *, "Build all factorisations"
   !Check for all combinations of n_decomp_max numbers between n_radix_max-1 radix except with 6
   !because 6 must be used only once and has to be the first one for the decomposition
   !We fill factor in rank i with the radix radix(rad)
   !The first factor factor(:,1) varies fast
   !The second factor factor(:,1) varies each 5 indexes (given by end1=5**(i-1)).
   !and so on
   do i=1,n_decomp_max-1
      ncount=1
      loop: do
         do rad=1,n_radix_max-1
            end1=(n_radix_max-1)**(i-1)
            if (i==1) then
               end1=1
            end if
            do j=1,end1
               factor(i,ncount)=radix(rad)
               ncount=ncount+1
               if (ncount>n_max) then
                  exit loop
               end if
            end do
         end do
      end do loop
   end do

   !Finally, we add all available radix (including 6) in factor(:,n_decomp_max)
   ncount=1
   loop6: do
      do rad=1,n_radix_max
         do j=1,n_max/n_radix_max
            factor(n_decomp_max,ncount)=radix(rad)
            ncount=ncount+1
            if (ncount>n_max) then
               exit loop6
            end if
         end do
      end do
   end do loop6

   !Now we calculate the numbers 1:ncount which are the multiplication of factors
   print *,'Doing the product'
   ncount=0
   do i=1,n_max
      reference=1
      do j=1,n_decomp_max
         reference=reference*factor(j,i)
      end do
      if (reference <= n_max) then
         ncount=ncount+1
         call sort_inplace(n_decomp_max,factor(:,i))
         selected(:,ncount)=factor(:,i)
         results(ncount)=reference
         !print *,ncount,'--',results(ncount),'selected',selected(:,ncount)
      end if
   end do
   n_element=ncount

   !Allocations
   deallocate(factor)
   allocate(indexes(n_element))


   print *,"Sort the results"
   !Sort the ncount results
   call sort_indexes(n_element,results,indexes)


   !Remove the same decompositions
   i=0
   loop_i: do
      i=i+1
      if (i>=n_element) then
         exit loop_i
      end if
      i_swap=indexes(i)
      if (results(i_swap) <= 0) then
         cycle loop_i
      end if
      loop_j: do j=i+1,n_element
         j_swap=indexes(j)
         !print *,i,j,results(i_swap), results(j_swap),'---',selected(:,i_swap),' ?? ',selected(:,j_swap)
         if (results(i_swap) == results(j_swap)) then
            do k=1,n_decomp_max
               if (selected(k,i_swap) /= selected(k,j_swap)) then
                  cycle loop_j
               end if
            end do
            !print *,j,results(j_swap), " == 0"
            results(j_swap) = -results(j_swap)
         else
            if (results(i_swap) /= abs(results(j_swap))) then
               cycle loop_i
            end if
         end if
      end do loop_j
   end do loop_i

   !Display
   do i=1,n_element
      i_swap=indexes(i)
      if (results(i_swap) > 0) then
         print *,results(i_swap),(selected(j,i_swap),j=n_decomp_max,1,-1)
      end if
   end do

   !De-Allocations
   deallocate(indexes)
   deallocate(results)
   deallocate(selected)

   contains


   !Algorithm not efficient but the number of indexes is small
   subroutine sort_inplace(n,array)
      implicit none
      !Arguements
      integer, intent(in) :: n
      integer, dimension(n), intent(inout) :: array
      !Local variables
      integer, dimension(1) :: i_loc
      integer :: i,i_swap,iadd
      do i=1,n
         i_loc = minloc(array(i:n))
         iadd = i - 1 + i_loc(1)
         if (i /= iadd) then
            i_swap = array(i)
            array(i) = array(iadd)
            array(iadd) = i_swap
         end if
      end do
   END SUBROUTINE sort_inplace


   !Algorithm not efficient but the number of indexes is small
   subroutine sort_indexes(n,array,indexes)
      implicit none
      !Arguments
      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: array
      integer, dimension(n), intent(out) :: indexes
      !Local variables
      integer, dimension(:), allocatable :: temp
      integer, dimension(1) :: i_loc
      integer :: i,i_max
      allocate(temp(n))
      temp(:)=array(:)
      i_max=huge(1)
      do i=1,n
         i_loc=minloc(temp)
         indexes(i)=i_loc(1)
         temp(i_loc(1))=i_max
      end do
      deallocate(temp)
   END SUBROUTINE sort_indexes

END PROGRAM fft_factors
