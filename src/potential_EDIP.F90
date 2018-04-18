module interface_edip
  use global
  use defs_basis

  implicit none

  private
  public :: &
    edip
 

contains


subroutine edip(parini,latvec_bohr,xred0,fxyz,strten,etot)
!      subroutine energyandforces(nat,latvec,rxyz0,fxyz,stress,pressure,etot,count1)
 


!   forces-edip.f
!   -------------

!   Version 1.0f

!   Force and Energy Calculation with the
!   Environment-Dependent Interatomic Potential

!   written by Martin Z. Bazant, 
!   Department of Physics, Harvard University
!   April - October 1997
!   (based on forces.c by Martin Z. Bazant, June 1994)

!   New address (2000):
!   Professor Martin Z. Bazant
!   Department of Mathematics 2-363B
!   Massachusetts Institute of Technology
!   Cambridge, MA 02139-4307

!   E-mail:
!   bazant@math.mit.edu

!   translated from c to FORTRAN 
!   by Noam Bernstein, noamb@cmt.harvard.edu, December 1997


!   COPYRIGHT NOTICE
!   ----------------

!   forces-edip, copyright 1997 by Martin Z. Bazant and Harvard University.
!   Permission is granted to use forces-edip for academic use only, at
!   no cost. Unauthorized sale or commerical use of this software
!   is prohibited by United States copyright law. Any publication describing
!   research involving this software should contain the following citations,
!   at once and in this order, to give proper credit to the theoretical work and
!   fitting that produced EDIP and this subroutine:

!     1.  M. Z. Bazant and E. Kaxiras, Phys. Rev. Lett. 77, 4370 (1996).
!     2.  M. Z. Bazant, E. Kaxiras, J. F. Justo, Phys. Rev. B 56, 8542 (1997).
!     3.  J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, and S. Yip, 
!	    Phys. Rev. B 58, 2539 (1998).

!   This software has been extensively tested for molecular dynamics simulations
!   on Sun, SGI and IBM architectures, but no guarantees are made.


!   WEBSITE
!   -------

!   Updated versions of this software are available at the EDIP distribution site,
!   http://pelion.eas.harvard.edu/software/EDIP/
!   Postscript files of related papers are available at the Kaxiras group web site
!   in the Department of Physics at Harvard University, http://pelion.eas.harvard.edu,
!   under 'Empirical Methods'.
!   A description of the algorithm used in this subroutine can be found
!   in the Ph.D. Thesis of M. Z. Bazant (1997), chapter 6, on the web at
!   http://pelion.eas.harvard.edu/~bazant/thesis/


!   INTERFACE
!   ---------

!     compute_forces_EDIP(N, p, lx, ly, lz, E, f)
!     N : number of particles
!     p : array (3,N) of positions in Angstroms
!     E : returned energy in eV
!     f : returned forces array (3,N) in eV/Angstroms

!     neighbors(p_nbrs(i)),...,neighbors(p_nbrs(i+1)) are the 
!     atoms which are "neighbors" of atom i, using a standard Verlet 
!     neighbor list. These are a global arrays, not passed, that are declared
!     in "edip_neighbors_include.h". This way of storing atomi! positions
!     is not unique, and will require a customized patch to the main MD program.

!     The parameters of the potential initialized in input_EDIP_params() are global
!     variables declared in "edip_pot_include.h".


!   PARAMETERS
!   ---------- 

!    par_cap_A,par_cap_B,par_rh,par_a,par_sig
!    par_lam,par_gam,par_b,par_c,par_delta
!    par_mu,par_Qo,par_palp,par_bet,par_alp

!    5.6714030     2.0002804     1.2085196     3.1213820     0.5774108
!    1.4533108     1.1247945     3.1213820     2.5609104    78.7590539
!    0.6966326   312.1341346     1.4074424     0.0070975     3.1083847

!    Connection between these parameters and those given in the paper,
!    Justo et al., Phys. Rev. B 58, 2539 (1998):

!    A((B/r)**rh-palp*exp(-bet*Z*Z)) = A'((B'/r)**rh-exp(-bet*Z*Z))

!    so in the paper (')
!    A' = A*palp
!    B' = B * palp**(-1/rh)
!    eta = detla/Qo

!    Non-adjustable parameters for tau(Z) from Ismail & Kaxiras, 1993,
!    also in Bazant, Kaxiras, Justo, PRB (1997):

!    u1 = -0.165799;
!    u2 = 32.557;
!    u3 = 0.286198;
!    u4 = 0.66;

!        subroutine compute_forces_EDIP(N_own, pos, L_x, L_y, L_z, 
!     &	    E_potential, f)
          use mod_parini, only: typ_parini
          implicit none
          type(typ_parini), intent(in):: parini
          
!  ------------------------- VARIABLE DECLARATIONS -------------------------
!          integer N_own, max_nbrs
          integer:: max_nbrs
!          double precision pos(3,N_own)
          real(8):: pos(3,parini%nat)
          real(8):: E_potential
!          double precision f(3,N_own)
          real(8):: f(3,parini%nat)
          real(8):: L_x, L_y, L_z

          integer i,j,k,l,n
          real(8):: dx,dy,dz,r,rsqr,asqr
          real(8):: rinv,rmainv,xinv,xinv3,den,Z,fZ
          real(8):: dV2j,dV2ijx,dV2ijy,dV2ijz,pZ,dp
          real(8):: temp0,temp1,temp2
          real(8):: Qort,muhalf,u5
          real(8):: rmbinv,winv,dwinv,tau,dtau,lcos,x,H,dHdx,dhdl
          real(8):: dV3rij,dV3rijx,dV3rijy,dV3rijz
          real(8):: dV3rik,dV3rikx,dV3riky,dV3rikz
          real(8):: dV3l,dV3ljx,dV3ljy,dV3ljz,dV3lkx,dV3lky,dV3lkz
          real(8):: dV2dZ,dxdZ,dV3dZ
          real(8):: dEdrl,dEdrlx,dEdrly,dEdrlz
          real(8):: bmc,cmbinv
          real(8):: fjx,fjy,fjz,fkx,fky,fkz
        
!          double precision s2_t0(MAX_NBRS_1)
!          double precision s2_t1(MAX_NBRS_1)
!          double precision s2_t2(MAX_NBRS_1)
!          double precision s2_t3(MAX_NBRS_1)
!          double precision s2_dx(MAX_NBRS_1)
!          double precision s2_dy(MAX_NBRS_1)
!          double precision s2_dz(MAX_NBRS_1)
!          double precision s2_r(MAX_NBRS_1)
          real(8),allocatable:: s2_t0(:)
          real(8),allocatable:: s2_t1(:)
          real(8),allocatable:: s2_t2(:)
          real(8),allocatable:: s2_t3(:)
          real(8),allocatable:: s2_dx(:)
          real(8),allocatable:: s2_dy(:)
          real(8),allocatable:: s2_dz(:)
          real(8),allocatable::  s2_r(:)

          integer n2                
!   size of s2[]
!          integer num2(MAX_NBRS_1)  
          integer,allocatable:: num2(:)  
!   atom ID numbers for s2[]
       
!          double precision s3_g(MAX_NBRS_1)
!          double precision s3_dg(MAX_NBRS_1)
!          double precision s3_rinv(MAX_NBRS_1)
!          double precision s3_dx(MAX_NBRS_1)
!          double precision s3_dy(MAX_NBRS_1)
!          double precision s3_dz(MAX_NBRS_1)
!          double precision s3_r(MAX_NBRS_1)
          real(8),allocatable::   s3_g(:)
          real(8),allocatable::  s3_dg(:)
          real(8),allocatable::s3_rinv(:)
          real(8),allocatable::  s3_dx(:)
          real(8),allocatable::  s3_dy(:)
          real(8),allocatable::  s3_dz(:)
          real(8),allocatable::   s3_r(:)

          integer n3                
!   size of s3[]
!          integer num3(MAX_NBRS_1)  
          integer,allocatable:: num3(:)  
!   atom ID numbers for s3[]
        
!          double precision sz_df(MAX_NBRS_1)
!          double precision sz_sum(MAX_NBRS_1)
!          double precision sz_dx(MAX_NBRS_1)
!          double precision sz_dy(MAX_NBRS_1)
!          double precision sz_dz(MAX_NBRS_1)
!          double precision sz_r(MAX_NBRS_1)

          real(8),allocatable::  sz_df(:)
          real(8),allocatable:: sz_sum(:)
          real(8),allocatable::  sz_dx(:)
          real(8),allocatable::  sz_dy(:)
          real(8),allocatable::  sz_dz(:)
          real(8),allocatable::   sz_r(:)
          integer nz                
!   size of sz[]
!          integer numz(MAX_NBRS_1)  
          integer,allocatable:: numz(:)  
!   atom ID numbers for sz[]
        
          integer nj,nk,nl         
!   indices for the store arrays
          real(8)::V2, V3, virial, virial_xyz(3)
          real(8)::L_x_div_2, L_y_div_2, L_z_div_2
          real(8)::coord_total
        
          integer fixZ, tricks_Zfix
          parameter (fixZ = 0, tricks_Zfix = 5)
          real(8)::    par_cap_A,par_cap_B,par_rh,par_a,par_sig
          real(8)::    par_lam,par_gam,par_b,par_c,par_delta
          real(8)::    par_mu,par_Qo,par_palp,par_bet,par_alp
          real(8)::    u1,u2,u3,u4
          real(8):: par_bg
          real(8):: par_eta
          real(8):: pot_cutoff
          real(8):: delta_safe

! My variables
          integer, allocatable:: lsta(:,:), lstb(:)
          real(8), allocatable:: rel(:,:)
          real(8):: etot,cut2,coord_iat,ener_iat,xarg,coord,coord2,ener,ener2
          real(8):: latvec(3,3),latvec_bohr(3,3),latvecinv(3,3),stress(3,3),gradvol(3,3),vol,rxyz(3,parini%nat),fxyz(3,parini%nat),trans(3,3),transinv(3,3),a(3,3),stressvol(3,3),xred0(3,parini%nat),xred(3,parini%nat),tmplat(3,3)
          integer:: nnbrx,istop,kk,ll
          real(8):: sil_sjl,si1_sj1,si2_sj2,si3_sj3,tkmsum,tt,strten(6)
          xred=xred0
!         End my variables
! EDIP parameters
!         taken from Justo et al., Phys. Rev. B 58, 2539 (1998).
          par_cap_A = 5.6714030d0
          par_cap_B = 2.0002804d0
          par_rh = 1.2085196d0
          par_a = 3.1213820d0
          par_sig = 0.5774108d0
          par_lam = 1.4533108d0
          par_gam = 1.1247945d0
          par_b = 3.1213820d0
          par_c = 2.5609104d0
          par_delta = 78.7590539d0
          par_mu = 0.6966326d0
          par_Qo = 312.1341346d0
          par_palp = 1.4074424d0
          par_bet = 0.0070975d0
          par_alp = 3.1083847d0

          par_bg=par_a
          par_eta = par_delta/par_Qo
          pot_cutoff = par_a
          delta_safe = 0.2d0

          u1 = -0.165799d0
          u2 = 32.557d0
          u3 = 0.286198d0
          u4 = 0.66d0
!end parameters

!Do some preparation including the construction of the pair list 
          nnbrx=50 !number of geighbors
          max_nbrs=nnbrx
          allocate(lsta(2,parini%nat),lstb(nnbrx*parini%nat),rel(5,nnbrx*parini%nat))
          latvec=latvec_bohr*Bohr_Ang
          trans=latvec
          cut2=par_a!-1.d-8
          call backtocell(parini%nat,latvec,xred)
          call rxyz_int2cart(latvec,xred,rxyz,parini%nat)
          call invertmat(trans,transinv,3)
          call makepairlist(lsta,lstb,rel,parini%nat,rxyz,latvec,cut2,nnbrx)


!Allocation of temporary arrays
          allocate(  s2_t0(nnbrx)) 
          allocate(  s2_t1(nnbrx))
          allocate(  s2_t2(nnbrx))
          allocate(  s2_t3(nnbrx))
          allocate(  s2_dx(nnbrx))
          allocate(  s2_dy(nnbrx))
          allocate(  s2_dz(nnbrx))
          allocate(   s2_r(nnbrx))
          allocate(   num2(nnbrx))
          allocate(   s3_g(nnbrx))
          allocate(  s3_dg(nnbrx))
          allocate(s3_rinv(nnbrx))
          allocate(  s3_dx(nnbrx))
          allocate(  s3_dy(nnbrx))
          allocate(  s3_dz(nnbrx))
          allocate(   s3_r(nnbrx))
          allocate(   num3(nnbrx))
          allocate(  sz_df(nnbrx))
          allocate( sz_sum(nnbrx))
          allocate(  sz_dx(nnbrx))
          allocate(  sz_dy(nnbrx))
          allocate(  sz_dz(nnbrx))
          allocate(   sz_r(nnbrx))
          allocate(   numz(nnbrx))
 

!          L_x_div_2 = L_x/2.0D0
!          L_y_div_2 = L_y/2.0D0
!          L_z_div_2 = L_z/2.0D0
        
!          do i=1, N_own
          do i=1, parini%nat
            f(1,i) = 0.0d0
            f(2,i) = 0.0d0
            f(3,i) = 0.0d0
          end do
          coord=0.d0 
          coord2=0.d0
          coord_total=0.0d0
          ener=0.d0
          ener2=0.d0


          virial=0.0d0
          virial_xyz(1)= 0.0d0
          virial_xyz(2)= 0.0d0
          virial_xyz(3)= 0.0d0
          V2=0.0d0
          V3=0.0d0
          stress(:,:)=0.d0      
          
!   COMBINE COEFFICIENTS
        
          asqr = par_a*par_a
          Qort = sqrt(par_Qo)
          muhalf = par_mu*0.5D0
          u5 = u2*u4
          bmc = par_b-par_c
          cmbinv = 1.0D0/(par_c-par_b)
        
        
          
!  --- LEVEL 1: OUTER LOOP OVER ATOMS ---
        
          do i=1, parini%nat 
            
!   RESET COORDINATION AND NEIGHBOR NUMBERS
            coord_iat=0.d0
            ener_iat=0.d0
            Z = 0.0d0
            n2 = 1
            n3 = 1
            nz = 1
        
            
!  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---
        
!            do n=p_nbrs(i), p_nbrs(i+1)-1
!              j = neighbors(n)
             do n=lsta(1,i),lsta(2,i)
              j=lstb(n)
        
              
!!   TEST IF WITHIN OUTER CUTOFF
!        
!              dx = pos(1,j)-pos(1,i)
!              if (dx .gt. L_x_div_2) then
!        	dx = dx - L_x
!              else if (dx .lt. -L_x_div_2) then
!        	dx = dx + L_x
!              end if 
!!  dx periodic
!              if(dabs(dx) < par_a) then
!              dy = pos(2,j)-pos(2,i)
!              if (dy .gt. L_y_div_2) then
!        	dy = dy - L_y
!              else if (dy .lt. -L_y_div_2) then
!        	dy = dy + L_y
!              end if 
!!  dy periodic
!              if(dabs(dy) < par_a) then
!              dz = pos(3,j)-pos(3,i)
!              if (dz .gt. L_z_div_2) then
!        	dz = dz - L_z
!              else if (dz .lt. -L_z_div_2) then
!        	dz = dz + L_z
!              end if 
!!  dy periodic
!              if(dabs(dz) < par_a) then
!              rsqr = dx*dx + dy*dy + dz*dz
!              if(rsqr < asqr) then
!                r = sqrt(rsqr)
        
                
!   PARTS OF TWO-BODY INTERACTION r<par_a
        
                num2(n2) = j
!                rinv = 1.0/r
!                dx = dx * rinv
!                dy = dy * rinv
!                dz = dz * rinv
                dx = -rel(1,n)
                dy = -rel(2,n)
                dz = -rel(3,n)
                r=rel(4,n)
                rinv=rel(5,n)

                rmainv = 1.0/(r-par_a)
                s2_t0(n2) = par_cap_A*dexp(par_sig*rmainv)
                s2_t1(n2) = (par_cap_B*rinv)**par_rh
                s2_t2(n2) = par_rh*rinv
                s2_t3(n2) = par_sig*rmainv*rmainv
                s2_dx(n2) = dx
                s2_dy(n2) = dy
                s2_dz(n2) = dz
                s2_r(n2) = r
                n2 = n2 + 1

                if (n2.gt.max_nbrs) then
                write(*,*) 'WARNING enlarge max_nbrs'
                istop=1
                return
                endif

!Additional part from stefan
! coordination number calculated with soft cutoff between first and
! second nearest neighbor
        if (r.le.2.36d0) then
        coord_iat=coord_iat+1.d0
        else if (r.ge.3.83d0) then
        else
        xarg=(r-2.36d0)*(1.d0/(3.83d0-2.36d0))
        coord_iat=coord_iat+(2*xarg+1.d0)*(xarg-1.d0)**2
        endif
!-----------------------------
                
!   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b
        
                if(r < par_bg)  then
        
                  num3(n3) = j
                  rmbinv = 1.0d0/(r-par_bg)
                  temp1 = par_gam*rmbinv
                  temp0 = dexp(temp1)
                  s3_g(n3) = temp0
                  s3_dg(n3) = -rmbinv*temp1*temp0
                  s3_dx(n3) = dx
                  s3_dy(n3) = dy
                  s3_dz(n3) = dz
                  s3_rinv(n3) = rinv
                  s3_r(n3) = r
                  n3 = n3 + 1
        
!        	  if(fixZ .eq. 0) then
!Additional part from Stefan
                  if (n3.gt.max_nbrs) then
                  write(*,*) 'WARNING enlarge max_nbrs'
                  istop=1
                  return
                  endif
!--------------------------      
                  
!   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b
        
                  if(r .lt. par_b) then
                    if(r .lt. par_c) then
                    Z = Z + 1.0d0
                   else
                    xinv = bmc/(r-par_c)
                    xinv3 = xinv*xinv*xinv
                    den = 1.0d0/(1.d0 - xinv3)
                    temp1 = par_alp*den
                    fZ = dexp(temp1)
                    Z = Z + fZ
                    numz(nz) = j
                    sz_df(nz) = fZ*temp1*den*3.0d0*xinv3*xinv*cmbinv   
!   df/dr
                    sz_dx(nz) = dx
                    sz_dy(nz) = dy
                    sz_dz(nz) = dz
                    sz_r(nz) = r
                    nz = nz + 1
!Additional part from Stefan
                    if (nz.gt.max_nbrs) then
                    write(*,*) 'WARNING enlarge max_nbrs'
                    istop=1
                    return
                    endif
!---------------------
                   end if 
!  r < par_C
                  end if 
!  r < par_b
                  end if 
!  fixZ .eq. 0
!                end if 
!  r < par_bg
!              end if 
!  rsqr < asqr
!              end if 
!  dz < par_a
!              end if 
!  dy < par_a
!              end if 
!  dz < par_a
            end do
        
!            if(fixZ .ne. 0) then
!        
!              Z = tricks_Zfix
!              coord_total = coord_total + Z
!              pZ = par_palp*dexp(-par_bet*Z*Z)
!              dp = 0.0
        
!            else
        
!              coord_total = coord_total + Z
        
              
!   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES
        
              do nl=1, nz-1
               sz_sum(nl)=0.0d0
              end do
        
              
!   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION
        
              temp0 = par_bet*Z
              pZ = par_palp*dexp(-temp0*Z)         
!   bond order
              dp = -2.0d0*temp0*pZ         
!   derivative of bond order
        
!            end if
        
            
!  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---
        
            do nj=1, n2-1
        
              temp0 = s2_t1(nj) - pZ
        
              
!   two-body energy V2(rij,Z)
        
!              V2 = V2 + temp0*s2_t0(nj)
              ener_iat = ener_iat + temp0*s2_t0(nj)
              
!   two-body forces
        
              dV2j = - s2_t0(nj) * (s2_t1(nj)*s2_t2(nj)  + temp0 * s2_t3(nj))   
!   dV2/dr
              dV2ijx = dV2j * s2_dx(nj)
              dV2ijy = dV2j * s2_dy(nj)
              dV2ijz = dV2j * s2_dz(nj)
              f(1,i) = f(1,i) + dV2ijx
              f(2,i) = f(2,i) + dV2ijy
              f(3,i) = f(3,i) + dV2ijz
              j = num2(nj)
              f(1,j) = f(1,j) - dV2ijx
              f(2,j) = f(2,j) - dV2ijy
              f(3,j) = f(3,j) - dV2ijz
        
              
!   dV2/dr contribution to virial
        
              virial_xyz(1) = virial_xyz(1) - s2_r(nj)*(dV2ijx*s2_dx(nj))
              virial_xyz(2) = virial_xyz(2) - s2_r(nj)*(dV2ijy*s2_dy(nj))
              virial_xyz(3) = virial_xyz(3) - s2_r(nj)*(dV2ijz*s2_dz(nj))
              virial = virial - s2_r(nj) * (dV2ijx*s2_dx(nj)+ dV2ijy*s2_dy(nj) + dV2ijz*s2_dz(nj))

                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*s2_dx(nj)+transinv(1,2)*s2_dy(nj)+transinv(1,3)*s2_dz(nj)
                   si2_sj2=transinv(2,1)*s2_dx(nj)+transinv(2,2)*s2_dy(nj)+transinv(2,3)*s2_dz(nj)
                   si3_sj3=transinv(3,1)*s2_dx(nj)+transinv(3,2)*s2_dy(nj)+transinv(3,3)*s2_dz(nj)
                   stress(1,1)=stress(1,1)-s2_r(nj)*dV2ijx*si1_sj1
                   stress(1,2)=stress(1,2)-s2_r(nj)*dV2ijx*si2_sj2
                   stress(1,3)=stress(1,3)-s2_r(nj)*dV2ijx*si3_sj3
                   stress(2,1)=stress(2,1)-s2_r(nj)*dV2ijy*si1_sj1
                   stress(2,2)=stress(2,2)-s2_r(nj)*dV2ijy*si2_sj2
                   stress(2,3)=stress(2,3)-s2_r(nj)*dV2ijy*si3_sj3
                   stress(3,1)=stress(3,1)-s2_r(nj)*dV2ijz*si1_sj1
                   stress(3,2)=stress(3,2)-s2_r(nj)*dV2ijz*si2_sj2
                   stress(3,3)=stress(3,3)-s2_r(nj)*dV2ijz*si3_sj3

        
!              if(fixZ .eq. 0) then
        
              
!  --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---
        
              dV2dZ = - dp * s2_t0(nj)
              do nl=1, nz-1
                 sz_sum(nl) =  sz_sum(nl) + dV2dZ
              end do
      
!              end if 
!  fixZ
            end do
!Commented out by Stefan        
!            if(fixZ .ne. 0) then
!              winv = Qort*dexp(-muhalf*Z)
!              dwinv = 0.0
!              temp0 = dexp(-u4*Z)
!              tau = u1+u2*temp0*(u3-temp0)
!              dtau = 0.0
!            else
!----------------------------------
              
!   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION
        
              winv = Qort*exp(-muhalf*Z) 
!   inverse width of angular function
              dwinv = -muhalf*winv       
!   its derivative
              temp0 = exp(-u4*Z)
              tau = u1+u2*temp0*(u3-temp0) 
!   -cosine of angular minimum
              dtau = u5*temp0*(2.d0*temp0-u3) 
!   its derivative
!            end if
        
            
!  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---
        
            do nj=1, n3-2
        
              j=num3(nj)
        
              
!  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---
        
              do nk=nj+1, n3-1
        
                k=num3(nk)
        
       
!   angular function h(l,Z)
        
                lcos = s3_dx(nj) * s3_dx(nk) + s3_dy(nj) * s3_dy(nk) + s3_dz(nj) * s3_dz(nk)
                x = (lcos + tau)*winv
                temp0 = exp(-x*x)
        
                H = par_lam*(1.d0 - temp0 + par_eta*x*x)
                dHdx = 2.d0*par_lam*x*(temp0 + par_eta)
        
                dhdl = dHdx*winv
        
        
!   three-body energy
        
                temp1 = s3_g(nj) * s3_g(nk)
!                V3 = V3 + temp1*H
                ener_iat = ener_iat + temp1*H
        
        
!   (-) radial force on atom j
        
                dV3rij = s3_dg(nj) * s3_g(nk) * H
                dV3rijx = dV3rij * s3_dx(nj)
                dV3rijy = dV3rij * s3_dy(nj)
                dV3rijz = dV3rij * s3_dz(nj)
                fjx = dV3rijx
                fjy = dV3rijy
                fjz = dV3rijz
                
                
!   (-) radial force on atom k
        
                dV3rik = s3_g(nj) * s3_dg(nk) * H
                dV3rikx = dV3rik * s3_dx(nk)
                dV3riky = dV3rik * s3_dy(nk)
                dV3rikz = dV3rik * s3_dz(nk)
                fkx = dV3rikx
                fky = dV3riky
                fkz = dV3rikz
                
                
!   (-) angular force on j
        
                dV3l = temp1*dhdl
                dV3ljx = dV3l * (s3_dx(nk) - lcos * s3_dx(nj)) * s3_rinv(nj)
                dV3ljy = dV3l * (s3_dy(nk) - lcos * s3_dy(nj)) * s3_rinv(nj)
                dV3ljz = dV3l * (s3_dz(nk) - lcos * s3_dz(nj)) * s3_rinv(nj)
                fjx = fjx + dV3ljx
                fjy = fjy + dV3ljy
                fjz = fjz + dV3ljz
                
       
!   (-) angular force on k
        
                dV3lkx = dV3l * (s3_dx(nj) - lcos * s3_dx(nk)) * s3_rinv(nk)
                dV3lky = dV3l * (s3_dy(nj) - lcos * s3_dy(nk)) * s3_rinv(nk)
                dV3lkz = dV3l * (s3_dz(nj) - lcos * s3_dz(nk)) * s3_rinv(nk)
                fkx = fkx + dV3lkx
                fky = fky + dV3lky
                fkz = fkz + dV3lkz
       
       
!   apply radial + angular forces to i, j, k
        
                f(1,j) = f(1,j) - fjx
                f(2,j) = f(2,j) - fjy
                f(3,j) = f(3,j) - fjz
                f(1,k) = f(1,k) - fkx
                f(2,k) = f(2,k) - fky
                f(3,k) = f(3,k) - fkz
                f(1,i) = f(1,i) + fjx + fkx
                f(2,i) = f(2,i) + fjy + fky
                f(3,i) = f(3,i) + fjz + fkz
       
!   dV3/dR contributions to virial
        
                virial = virial - s3_r(nj) * (fjx*s3_dx(nj) + fjy*s3_dy(nj) + fjz*s3_dz(nj))
                virial = virial - s3_r(nk) * (fkx*s3_dx(nk) + fky*s3_dy(nk) + fkz*s3_dz(nk))
                virial_xyz(1) = virial_xyz(1) - s3_r(nj)*(fjx*s3_dx(nj))
                virial_xyz(2) = virial_xyz(2) - s3_r(nj)*(fjy*s3_dy(nj))
                virial_xyz(3) = virial_xyz(3) - s3_r(nj)*(fjz*s3_dz(nj))
                virial_xyz(1) = virial_xyz(1) - s3_r(nk)*(fkx*s3_dx(nk))
                virial_xyz(2) = virial_xyz(2) - s3_r(nk)*(fky*s3_dy(nk))
                virial_xyz(3) = virial_xyz(3) - s3_r(nk)*(fkz*s3_dz(nk))


                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*s3_dx(nj)+transinv(1,2)*s3_dy(nj)+transinv(1,3)*s3_dz(nj)
                   si2_sj2=transinv(2,1)*s3_dx(nj)+transinv(2,2)*s3_dy(nj)+transinv(2,3)*s3_dz(nj)
                   si3_sj3=transinv(3,1)*s3_dx(nj)+transinv(3,2)*s3_dy(nj)+transinv(3,3)*s3_dz(nj)
                   stress(1,1)=stress(1,1)-s3_r(nj)*fjx*si1_sj1
                   stress(1,2)=stress(1,2)-s3_r(nj)*fjx*si2_sj2
                   stress(1,3)=stress(1,3)-s3_r(nj)*fjx*si3_sj3
                   stress(2,1)=stress(2,1)-s3_r(nj)*fjy*si1_sj1
                   stress(2,2)=stress(2,2)-s3_r(nj)*fjy*si2_sj2
                   stress(2,3)=stress(2,3)-s3_r(nj)*fjy*si3_sj3
                   stress(3,1)=stress(3,1)-s3_r(nj)*fjz*si1_sj1
                   stress(3,2)=stress(3,2)-s3_r(nj)*fjz*si2_sj2
                   stress(3,3)=stress(3,3)-s3_r(nj)*fjz*si3_sj3


                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*s3_dx(nk)+transinv(1,2)*s3_dy(nk)+transinv(1,3)*s3_dz(nk)
                   si2_sj2=transinv(2,1)*s3_dx(nk)+transinv(2,2)*s3_dy(nk)+transinv(2,3)*s3_dz(nk)
                   si3_sj3=transinv(3,1)*s3_dx(nk)+transinv(3,2)*s3_dy(nk)+transinv(3,3)*s3_dz(nk)
                   stress(1,1)=stress(1,1)-s3_r(nk)*fkx*si1_sj1
                   stress(1,2)=stress(1,2)-s3_r(nk)*fkx*si2_sj2
                   stress(1,3)=stress(1,3)-s3_r(nk)*fkx*si3_sj3
                   stress(2,1)=stress(2,1)-s3_r(nk)*fky*si1_sj1
                   stress(2,2)=stress(2,2)-s3_r(nk)*fky*si2_sj2
                   stress(2,3)=stress(2,3)-s3_r(nk)*fky*si3_sj3
                   stress(3,1)=stress(3,1)-s3_r(nk)*fkz*si1_sj1
                   stress(3,2)=stress(3,2)-s3_r(nk)*fkz*si2_sj2
                   stress(3,3)=stress(3,3)-s3_r(nk)*fkz*si3_sj3






        
!                if(fixZ .eq. 0) then
        
        
!   prefactor for 4-body forces from coordination
                  dxdZ = dwinv*(lcos + tau) + winv*dtau
                  dV3dZ = temp1*dHdx*dxdZ
       
       
!  --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---
        
                  do nl=1, nz-1
                    sz_sum(nl) = sz_sum(nl) + dV3dZ
                  end do
!                end if
              end do
            end do
        
!            if(fixZ .eq. 0) then
        
            
!  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---
        
            do nl=1, nz-1
        
                dEdrl = sz_sum(nl) * sz_df(nl)
                dEdrlx = dEdrl * sz_dx(nl)
                dEdrly = dEdrl * sz_dy(nl)
                dEdrlz = dEdrl * sz_dz(nl)
                f(1,i) = f(1,i) + dEdrlx
                f(2,i) = f(2,i) + dEdrly
                f(3,i) = f(3,i) + dEdrlz
                l = numz(nl)
                f(1,l) = f(1,l) - dEdrlx
                f(2,l) = f(2,l) - dEdrly
                f(3,l) = f(3,l) - dEdrlz
        
                
!   dE/dZ*dZ/dr contribution to virial
        
                virial = virial - sz_r(nl) * (dEdrlx*sz_dx(nl)+ dEdrly*sz_dy(nl) + dEdrlz*sz_dz(nl))
                virial_xyz(1) = virial_xyz(1) - sz_r(nl)*(dEdrlx*sz_dx(nl))
                virial_xyz(2) = virial_xyz(2) - sz_r(nl)*(dEdrly*sz_dy(nl))
                virial_xyz(3) = virial_xyz(3) - sz_r(nl)*(dEdrlz*sz_dz(nl))


                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*sz_dx(nl)+transinv(1,2)*sz_dy(nl)+transinv(1,3)*sz_dz(nl)
                   si2_sj2=transinv(2,1)*sz_dx(nl)+transinv(2,2)*sz_dy(nl)+transinv(2,3)*sz_dz(nl)
                   si3_sj3=transinv(3,1)*sz_dx(nl)+transinv(3,2)*sz_dy(nl)+transinv(3,3)*sz_dz(nl)
                   stress(1,1)=stress(1,1)-sz_r(nl)*dEdrlx*si1_sj1
                   stress(1,2)=stress(1,2)-sz_r(nl)*dEdrlx*si2_sj2
                   stress(1,3)=stress(1,3)-sz_r(nl)*dEdrlx*si3_sj3
                   stress(2,1)=stress(2,1)-sz_r(nl)*dEdrly*si1_sj1
                   stress(2,2)=stress(2,2)-sz_r(nl)*dEdrly*si2_sj2
                   stress(2,3)=stress(2,3)-sz_r(nl)*dEdrly*si3_sj3
                   stress(3,1)=stress(3,1)-sz_r(nl)*dEdrlz*si1_sj1
                   stress(3,2)=stress(3,2)-sz_r(nl)*dEdrlz*si2_sj2
                   stress(3,3)=stress(3,3)-sz_r(nl)*dEdrlz*si3_sj3



            end do
        
!           end if
        coord=coord+coord_iat
        coord2=coord2+coord_iat**2
        ener = ener + ener_iat
        ener2 = ener2 + ener_iat**2
          end do

        call getvol(latvec,vol)
          if(vol.lt.0.d0) then
            write(77,*) parini%nat
            write(77,*) latvec(:,1)
            write(77,*) latvec(:,2)
            write(77,*) latvec(:,3)
            do i=1,parini%nat
              write(77,*)  rxyz(:,i),"Si"
            enddo

          endif
         
          etot=ener
          fxyz=f

          etot=etot/Ha_eV
          fxyz=fxyz/Ha_eV*Bohr_Ang
!Stress
!Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
        do i=1,3
           tmplat(:,i)=latvec(i,:)
        enddo
        stress=-matmul(stress,tmplat)/vol
        strten(1) = stress(1,1)
        strten(2) = stress(2,2)
        strten(3) = stress(3,3)
        strten(6) = stress(2,1)
        strten(5) = stress(3,1)
        strten(4) = stress(3,2)
!This is not very clear yet...
        strten=strten/Ha_eV*Bohr_Ang**3
        end subroutine




 subroutine findinter(iat,nat,lstbiat,cut,rxyz,rxyzexp,transvecall,latvec,ninter,nmax,relinteriat,infocode) 
 !This subroutine will find all neighboring atoms of atom iat including those in periodic
 !images within the range of cut. The cell is described by latvec[v1,v2,v3]. rxyz holds the
 !coordinates of all atoms, nat is the number of atoms. The output ist ninter=the number of
 !atoms in the list of neighbors. Relative coordinates are stored in relinter,
 !the index of the interactiong atoms in lstbiat. rxyz holds the pos of all periodic images
 !infocode vill return 1 if an error occured and the number of neighbors exceeds nmax, else 0
 real*8,  intent(in)  :: rxyz(3,nat), latvec(3,3), cut, rxyzexp(3,nat,3,3,3),transvecall(3,3,3,3)
 real*8, intent(inout):: relinteriat(5,nmax)
 integer, intent(in)  :: iat,nat,nmax
 integer,intent(inout):: ninter,infocode,lstbiat(nmax)
 real(8)              :: a(3),b(3),crossp(3),d2plane(6),d2corner(8),d2edge(12),d2all(3,3,3),ppoint0(3),tempvec(3)
 real(8)              :: nvectmp(3),pointtmp(3),ppointtmp(3),latvectmp(3),dist,nvec(3,3),cutplus
 real(8)              :: tempvece(3,3,3,3,2),tempvecp(3,3,3,3,2)
 integer              :: k,l,m,summe,kt,lt,mt
 real(8)              :: rxyzp(3,nat)
 !For second nearest images
 real(8),allocatable  :: d2all_2(:,:,:),rxyzexp_2(:,:,:,:,:),transvecall_2(:,:,:,:)
 !For third nearest images
 real(8),allocatable  :: d2all_3(:,:,:),rxyzexp_3(:,:,:,:,:),transvecall_3(:,:,:,:)
 !For the fourth nearest images
 real(8),allocatable  :: d2all_4(:,:,:),rxyzexp_4(:,:,:,:,:),transvecall_4(:,:,:,:)
 real(8)              :: dx2,dy2,dz2,dd2,cut2,aa,bb,thrd
 integer              :: k1,l1,m1
 real(8)              :: dproj(6),rotmat(3,3)
 logical              :: nec


 ninter=0 
 ppoint0=(/0.d0,0.d0,0.d0/)
 cutplus=cut+cut*0.01d0
 d2all=cutplus
 tempvece=0.d0

 call nveclatvec(latvec,nvec)
    !Bottom
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall(:,1,1,2),dist)
    d2all(:,:,1)=dist
    !Top
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall(:,1,1,3),dist)
    d2all(:,:,3)=dist

    !Back
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall(:,2,1,1),dist)
    d2all(1,:,:)=dist
    !Front
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall(:,3,1,1),dist)
    d2all(3,:,:)=dist

    !Left
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall(:,1,2,1),dist)
    d2all(:,1,:)=dist
    !Right
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall(:,1,3,1),dist)
    d2all(:,3,:)=dist


 d2all(2,2,2)=0.d0 !Always calculate in original cell
 !Calculate all information of relinter for atoms within the cell
 ninter=0
 do k=1,3
  do l=1,3
   do m=1,3
      if (d2all(k,l,m).lt.cut) then
      
      rxyzp(:,:)=rxyzexp(:,:,k,l,m)
!      write(*,*) 'calling reliatp' ,iat,k,l,m
!      write(*,*) "ninter1",ninter
      call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
!      write(*,*) "ninter2",ninter
      endif  
   enddo
  enddo
 enddo
 infocode=0
! write(*,*) 'all well'

! return

 !If second layer neccesary
 !Call subroutine to check here...
 call check_rep(latvec,cut,nec)
 if(nec) then
    return
 else
!    write(*,*) "entering layer 2"
    allocate(d2all_2(5,5,5),rxyzexp_2(3,nat,5,5,5),transvecall_2(3,5,5,5))
    call expand_2(rxyz,rxyzexp_2,transvecall_2,latvec,nat)
! call latvec2dproj(dproj,latvec,rotmat,rxyz,nat) 
! open(unit=12,file="secondexp.ascii")
! write(12,*) nat
! write(12,*) dproj(1),dproj(2),dproj(3)
! write(12,*) dproj(4),dproj(5),dproj(6)
! do k=1,5
! do l=1,5
! do m=1,5
! do jat=1,nat
! write(12,*) rxyzexp_2(:,jat,k,l,m),"au"
! enddo
! enddo
! enddo
! enddo
! close(12)
!!! stop

    d2all_2=cutplus
    !Calculate interaction with second nearest cells.
    !Check planes
    !Bottom
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_2(:,1,1,2),dist)
    d2all_2(:,:,1)=dist
    !Top
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_2(:,1,1,5),dist)
    d2all_2(:,:,5)=dist

    !Back
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_2(:,2,1,1),dist)
    d2all_2(1,:,:)=dist
    !Front
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_2(:,5,1,1),dist)
    d2all_2(5,:,:)=dist

    !Left
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_2(:,1,2,1),dist)
    d2all_2(:,1,:)=dist
    !Right
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_2(:,1,5,1),dist)
    d2all_2(:,5,:)=dist

 
!    !Substract corners
!    bb=5.d0*(1.d0-4.d0/3.d0)
!    aa=4.d0/3.d0
!    cut2=cut*cut
!    do k=2,5,3
!       do l=2,5,3
!          do m=2,5,3
!          dx2=rxyz(1,iat)-transvecall_2(1,k,l,m) 
!          dx2=dx2*dx2
!          dy2=rxyz(2,iat)-transvecall_2(2,k,l,m)
!          dy2=dy2*dy2
!          dz2=rxyz(3,iat)-transvecall_2(3,k,l,m)
!          dz2=dz2*dz2
!          dd2=dx2+dy2+dz2
!             if(dd2.gt.cut2) then
!              k1=int(bb+aa*real(k,8))
!              l1=int(bb+aa*real(l,8))
!              m1=int(bb+aa*real(m,8))
!              d2all_2(k1,l1,m1)=cutplus
!             endif   
!          enddo
!       enddo
!    enddo 
   
    do k=1,5
     do l=1,5
      do m=1,5
         if (d2all_2(k,l,m).ge.cut) then
         goto 1010
         else                        
         rxyzp(:,:)=rxyzexp_2(:,:,k,l,m)
   !      write(*,*) 'calling reliatp' ,iat,k,l,m
         call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
         endif
         1010 continue
      enddo
     enddo
    enddo
    infocode=0
    deallocate(d2all_2,rxyzexp_2,transvecall_2)
    nec=.false.
    call check_rep(latvec,1.d0/2.d0*cut,nec)

    if(nec) then
       return
    else

!    write(*,*) "entering layer 3"
       !Calculate third layer of periodic images if neccesary 
       allocate(d2all_3(7,7,7),rxyzexp_3(3,nat,7,7,7),transvecall_3(3,7,7,7))
       call expand_3(rxyz,rxyzexp_3,transvecall_3,latvec,nat)
! call latvec2dproj(dproj,latvec,rotmat,rxyz,nat) 
! open(unit=12,file="secondexp3.ascii")
! write(12,*) nat
! write(12,*) dproj(1),dproj(2),dproj(3)
! write(12,*) dproj(4),dproj(5),dproj(6)
! do k=1,7
! do l=1,7
! do m=1,7
! do jat=1,nat
! write(12,*) rxyzexp_3(:,jat,k,l,m),"au"
! enddo
! enddo
! enddo
! enddo
! close(12)
! stop


!       write(*,*) "after expand 3"
       d2all_3=cutplus
       !Calculate interaction with third nearest cells.
       !Check planes
       !Bottom
       call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_3(:,1,1,2),dist)
       d2all_3(:,:,1)=dist
       !Top
       call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_3(:,1,1,7),dist)
       d2all_3(:,:,7)=dist
   
       !Back
       call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_3(:,2,1,1),dist)
       d2all_3(1,:,:)=dist
       !Front
       call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_3(:,7,1,1),dist)
       d2all_3(7,:,:)=dist
   
       !Left
       call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_3(:,1,2,1),dist)
       d2all_3(:,1,:)=dist
       !Right
       call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_3(:,1,7,1),dist)
       d2all_3(:,7,:)=dist
   
    
!       !Substract corners
!       bb=7.d0*(1.d0-6.d0/5.d0)
!       aa=6.d0/5.d0
!       cut2=cut*cut
!       do k=2,7,5
!          do l=2,7,5
!             do m=2,7,5
!             dx2=rxyz(1,iat)-transvecall_3(1,k,l,m) 
!             dx2=dx2*dx2
!             dy2=rxyz(2,iat)-transvecall_3(2,k,l,m)
!             dy2=dy2*dy2
!             dz2=rxyz(3,iat)-transvecall_3(3,k,l,m)
!             dz2=dz2*dz2
!             dd2=dx2+dy2+dz2
!                if(dd2.gt.cut2) then
!                 k1=int(bb+aa*real(k,8))
!                 l1=int(bb+aa*real(l,8))
!                 m1=int(bb+aa*real(m,8))
!                 d2all_3(k1,l1,m1)=cutplus
!                endif   
!             enddo
!          enddo
!       enddo 
      
       do k=1,7
        do l=1,7
         do m=1,7
            if (d2all_3(k,l,m).ge.cut) then
            goto 1011
            else                        
            rxyzp(:,:)=rxyzexp_3(:,:,k,l,m)
      !      write(*,*) 'calling reliatp' ,iat,k,l,m
            call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
            endif
            1011 continue
         enddo
        enddo
       enddo
       infocode=0
       deallocate(d2all_3,rxyzexp_3,transvecall_3)
       nec=.false.
       call check_rep(latvec,1.d0/3.d0*cut,nec)
       if(nec) then
       return
       else
          !Calculate fourth layer of periodic images if neccesary 
          allocate(d2all_4(9,9,9),rxyzexp_4(3,nat,9,9,9),transvecall_4(3,9,9,9))
          call expand_4(rxyz,rxyzexp_4,transvecall_4,latvec,nat)
          d2all_4=cutplus
          !Calculate interaction with fourth nearest cells.
          !Check planes
          !Bottom
          call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_4(:,1,1,2),dist)
          d2all_4(:,:,1)=dist
          !Top
          call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_4(:,1,1,9),dist)
          d2all_4(:,:,9)=dist
   
          !Back
          call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_4(:,2,1,1),dist)
          d2all_4(1,:,:)=dist
          !Front
          call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_4(:,9,1,1),dist)
          d2all_4(9,:,:)=dist
   
          !Left
          call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_4(:,1,2,1),dist)
          d2all_4(:,1,:)=dist
          !Right
          call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_4(:,1,9,1),dist)
          d2all_4(:,9,:)=dist


          do k=1,9
           do l=1,9
            do m=1,9
               if (d2all_4(k,l,m).ge.cut) then
               goto 1012
               else
               rxyzp(:,:)=rxyzexp_4(:,:,k,l,m)
         !      write(*,*) 'calling reliatp' ,iat,k,l,m
               call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
               endif
               1012 continue
            enddo
           enddo
          enddo
          infocode=0
          deallocate(d2all_4,rxyzexp_4,transvecall_4)
          nec=.false.
          call check_rep(latvec,1.d0/4.d0*cut,nec)
          if (nec) then
          return
          else
          write(*,*)  "Cell periodicit exceeded the fourth order"
          stop
          endif
      endif 

    endif


 endif

 end subroutine


 subroutine expand(rxyz,rxyzout,transvecall,latvec,nat)
 !This subroutine will expand the unit cell into 26 periodic cells and store them in rxyzout
 implicit none
 real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
 integer, intent(in) :: nat
 real*8, intent(out) :: rxyzout(3,nat,3,3,3) !26 periodic images plus the main cell
 integer             :: iat,iplane,icorner,iedge,m,k,l
 real*8,intent(inout):: transvecall(3,3,3,3)!,(transvecp(3,6),transvecc(3,8),transvece(3,12)

 do m=-1,1
    do k=-1,1
       do l=-1,1
       transvecall(:,l+2,k+2,m+2)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 do m=1,3
    do k=1,3
       do l=1,3
       do iat=1,nat
       rxyzout(:,iat,m,k,l)=rxyz(:,iat)+transvecall(:,m,k,l)
       enddo
       enddo
    enddo
 enddo
 end subroutine

 
 subroutine expand_2(rxyz,rxyzout_2,transvecall_2,latvec,nat)
 !This subroutine will expand the periodic images to the second order, thus creating the periodic cells of 5*5*5 images
 !The 3*3*3 centered cells are nulled bc thay are present in the expansion of the first order (first neighboring cells)
 implicit none
 real(8), intent(in)   :: latvec(3,3),rxyz(3,nat)
 real(8), intent(inout):: rxyzout_2(3,nat,5,5,5)
 integer, intent(in)   :: nat
 integer               :: iat,iplane,icorner,iedge,m,k,l
 real(8),intent(inout) :: transvecall_2(3,5,5,5)

 do m=-2,2
    do k=-2,2
       do l=-2,2
       transvecall_2(:,l+3,k+3,m+3)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3) 
       enddo
    enddo
 enddo

 rxyzout_2=0.d0
 do m=1,5
    do k=1,5
       do l=1,5
       do iat=1,nat
!       if(m==5 .or. k==5 .or. l==5 .or. m==1 .or. k==1 .or. l==1) then
       rxyzout_2(:,iat,m,k,l)=rxyz(:,iat)+transvecall_2(:,m,k,l)
!       endif
       enddo
       enddo
    enddo
 enddo

 end subroutine

 subroutine expand_3(rxyz,rxyzout_3,transvecall_3,latvec,nat)
 !This subroutine will expand the periodic images to the second order, thus creating the periodic cells of 5*5*5 images
 !The 3*3*3 centered cells are nulled bc thay are present in the expansion of the first order (first neighboring cells)
 implicit none
 real(8), intent(in)   :: latvec(3,3),rxyz(3,nat)
 real(8), intent(inout):: rxyzout_3(3,nat,7,7,7)
 integer, intent(in)   :: nat
 integer               :: iat,iplane,icorner,iedge,m,k,l
 real(8),intent(inout) :: transvecall_3(3,7,7,7)

 do m=-3,3
    do k=-3,3
       do l=-3,3
       transvecall_3(:,l+4,k+4,m+4)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 rxyzout_3=0.d0
 do m=1,7
    do k=1,7
       do l=1,7
       do iat=1,nat
!       if(m==7 .or. k==7 .or. l==7 .or. m==1 .or. k==1 .or. l==1) then
       rxyzout_3(:,iat,m,k,l)=rxyz(:,iat)+transvecall_3(:,m,k,l)
!       endif
       enddo
       enddo
    enddo
 enddo

 end subroutine


 subroutine expand_4(rxyz,rxyzout_4,transvecall_4,latvec,nat)
 !This subroutine will expand the periodic images to the second order, thus creating the periodic cells of 5*5*5 images
 !The 3*3*3 centered cells are nulled bc thay are present in the expansion of the first order (first neighboring cells)
 implicit none
 real(8), intent(in)   :: latvec(3,3),rxyz(3,nat)
 real(8), intent(inout):: rxyzout_4(3,nat,9,9,9)
 integer, intent(in)   :: nat
 integer               :: iat,iplane,icorner,iedge,m,k,l
 real(8),intent(inout) :: transvecall_4(3,9,9,9)

 do m=-4,4
    do k=-4,4
       do l=-4,4
       transvecall_4(:,l+5,k+5,m+5)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 rxyzout_4=0.d0
 do m=1,9
    do k=1,9
       do l=1,9
       do iat=1,nat
!       if(m==7 .or. k==7 .or. l==7 .or. m==1 .or. k==1 .or. l==1) then
       rxyzout_4(:,iat,m,k,l)=rxyz(:,iat)+transvecall_4(:,m,k,l)
!       endif
       enddo
       enddo
    enddo
 enddo

 end subroutine




 subroutine reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
 !This subroutine will calculate neighbors of atom iat within the given periodic image cell and the relative coordinates
 !The relative informations are stored in relinteriat
 !The atomic index jat is stored in lstbiat
 !rxyz holds the positions of iat, rxyzp is one of the periodic images to check
 !d2plane,d2corner,d2egde hold the distance of atom iat to the walls of the cell
 real*8, intent(inout)  :: relinteriat(5,nmax)
 real*8, intent(in)     :: rxyz(3,nat),rxyzp(3,nat),cut
 integer, intent(in)    :: nat,iat,nmax
 integer, intent(inout) :: ninter,lstbiat(nmax)
 real*8                 :: cut2,xrel,yrel,zrel,rr2,tt,tti
 integer                :: jat

 cut2=cut*cut
 do jat=1,nat
 xrel= rxyz(1,iat)-rxyzp(1,jat)
 yrel= rxyz(2,iat)-rxyzp(2,jat)
 zrel= rxyz(3,iat)-rxyzp(3,jat)
 rr2=xrel**2 + yrel**2 + zrel**2
! if (rr2.le.cut2 .AND. iat.ne.jat) then
 if (rr2.lt.cut2 .AND. rr2.ne.0.d0) then
! write(*,*) xrel,yrel,zrel
  ninter=ninter+1
  if (ninter.gt.nmax) stop 'nmax exceeded'
  lstbiat(ninter)=jat
  tt=dsqrt(rr2)
  tti=1.d0/tt
  relinteriat(1,ninter)=xrel*tti
  relinteriat(2,ninter)=yrel*tti
  relinteriat(3,ninter)=zrel*tti
  relinteriat(4,ninter)=tt
  relinteriat(5,ninter)=tti
  
!  write(*,*) 'Inside reliatp,jat found',jat,tt 
 endif
 enddo
  
 
 return
 end subroutine

 subroutine makepairlist(lsta,lstb,rel,nat,rxyz,latvec,cut,nmax)
 !This subroutine will construct the pairlist lsta and lstb of the atomic positions rxyz in the cell latvec
 !The cutoff is given with cut, all relative coordinates are stored in rel
 implicit none
 integer, intent(in) :: nat
 integer, intent(inout) :: lsta(2,nat),lstb(nat*nmax),nmax
 real*8, intent(inout)  :: rel(5,nat*nmax)
 real*8, intent(in)     :: cut, rxyz(3,nat), latvec(3,3)
 real*8                 :: rxyzexp(3,nat,3,3,3),relinteriat(5,nmax),transvecall(3,3,3,3)
 integer                :: iat,lstbiat(nmax),infocode,ninter,k,l,m
 
 call expand(rxyz,rxyzexp,transvecall,latvec,nat) 
 do iat=1,nat
 ninter=0
    lstbiat=0
    relinteriat=0.d0
    call findinter(iat,nat,lstbiat,cut,rxyz,rxyzexp,transvecall,latvec,ninter,nmax,relinteriat,infocode)
    lsta(1,iat)=(iat-1)*nmax+1
    lsta(2,iat)=(iat-1)*nmax+ninter
    lstb(lsta(1,iat):lsta(1,iat)+nmax-1)=lstbiat(:)
    rel(:,lsta(1,iat):lsta(1,iat)+nmax-1)=relinteriat(:,:)
 enddo
 
! do iat=1,nat
!    k=lsta(1,iat)
!    l=lsta(2,iat)
!    do m=k,l
!    write(*,*) iat,lstb(m),m
!    enddo
! enddo
end subroutine


 subroutine checkcellsize(latvec,cut)
 !This subroutine will check if the cell size is sufficentli large to have periodic boundary conditions
 !with only one level of cell reproduction
 implicit none
 real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist
 integer:: i
 call nveclatvec(latvec,nvec)
 point0=(/0.d0,0.d0,0.d0/)
 do i=1,3
 call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist) 
! write(*,*) i,cut, dist
! if (dist.le.cut) stop "Cell too sheared for periodicity"
 if (dist.le.cut) then
 write(*,*) "Cell too sheared for periodicity"
 write(*,*) latvec
 stop
 endif
 enddo 
 end subroutine



 subroutine check_rep(latvec,cut,nec)
 !This subroutine will check if the cell size is sufficentli large to have periodic boundary conditions
 !with for the given cut. If true, it fits. If false, it does not fit--> new layer neccesary 
 implicit none
 real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist,eps
 integer:: i
 logical:: nec
! eps=1.d-6
 nec=.true.
 call nveclatvec(latvec,nvec)
 point0=(/0.d0,0.d0,0.d0/)
 do i=1,3
 call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist)
! write(*,*) "cut",i,cut, dist
 if (dist.le.cut) nec=.false.
 enddo
 end subroutine


end module
subroutine init_edip(parini)
use mod_parini, only: typ_parini
use global
implicit none
type(typ_parini), intent(in):: parini
integer:: iat
do iat=1,parini%nat
   if(int(parini%znucl(parini%typat_global(iat))).ne.14) then
     stop "EDIP only allowed with Silicon"
   endif
enddo
end subroutine
