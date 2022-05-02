subroutine random_incell_p5(LATSGP,CRYSSYS,NGUESS,RED_POS)
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS, I
real(8):: RED_POS(3,NGUESS)
  select case(LATSGP)
  case(141)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=1.0d0/8.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.25d0
     RED_POS(3,4)=5.0d0/8.0d0

     do i=5,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,i)=0.125d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)!call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
!  case(141)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.0d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=3.0d0/8.0d0
!     RED_POS(1,4)=0.0d0
!     RED_POS(2,4)=0.75d0
!     RED_POS(3,4)=0.125d0
!
!     do i=5,1*int(NGUESS/6)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)+0.25d0
!        RED_POS(3,i)=7.0/8.0d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
!        RED_POS(1,i)=0.0d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/6)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 7",i,LATSGP
!     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(142)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.25d0
     RED_POS(3,2)=0.125d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=1.0d0/8.0d0

     do i=4,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=1.0d0/8.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
!  case(142)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.125d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=3.0d0/8.0d0
!
!     do i=4,1*int(NGUESS/4)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)+0.25d0
!        RED_POS(3,i)=1.0d0/8.0d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!
!     do i=3*int(NGUESS/4)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(143)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1.0d0/3.0d0
        RED_POS(2,i)=2.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=2.0d0/3.0d0
        RED_POS(2,i)=1.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (144)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case(145)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case(146) !hexagonal axis  
     do i=1,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(147)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=1.0d0/3.0d0
        RED_POS(2,i)=2.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(148) !hexagonal axis
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(149)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1.0d0/3.0d0
     RED_POS(2,3)=2.0d0/3.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1.0d0/3.0d0
     RED_POS(2,4)=2.0d0/3.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=2.0d0/3.0d0
     RED_POS(2,5)=1.0d0/3.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=2.0d0/3.0d0
     RED_POS(2,6)=1.0d0/3.0d0
     RED_POS(3,6)=0.5d0!0.0d0
     do i=7,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1.0d0/3.0d0
        RED_POS(2,i)=2.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=2.0d0/3.0d0
        RED_POS(2,i)=1.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case(150)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=1.0d0/3.0d0
        RED_POS(2,i)=2.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(151)
     do i=1,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=1.0d0/3.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=5.0d0/6.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(152)
     do i=1,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=1.0d0/3.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=5.0d0/6.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS  !3*int(NGUESS/3)   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(153)
     do i=1,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=2.0d0/3.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=1.0d0/6.0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(154)
     do i=1,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=1.0d0/6.0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=2.0d0/3.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(155) !hexagonal axis
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(156)
     do i=1,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=1.0d0/3.0d0
        RED_POS(2,i)=2.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=2.0d0/3.0d0
        RED_POS(2,i)=1.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(157)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1.0d0/3.0d0
        RED_POS(2,i)=2.0d0/3.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (158)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=2/3
        RED_POS(2,i)=1/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (159)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (160) !hexagonal axis
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS!/3   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (161) !hexagonal axis
     do i=1,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case (162)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case (163)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=2/3
     RED_POS(2,4)=1/3
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     do i=6,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (164)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case (165)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     do i=4,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (166) !hexagonal axis
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     end do
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case (167) !hexagonal axis
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     do i=4,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     end do
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (168)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)   !BUG?
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (169)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case (170)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case (171)

     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (172)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (173)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (174)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=2/3
     RED_POS(2,5)=1/3
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=2/3
     RED_POS(2,6)=1/3
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=2/3
        RED_POS(2,i)=1/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case (175)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case (176)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=2/3
     RED_POS(2,4)=1/3
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     do i=6,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (177)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     end do
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     end do
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 8",i,LATSGP
     end do
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
  case (178)
     do i=1,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     end do
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     end do
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (179)
     do i=1,int(NGUESS/3)!1
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0!*int(NGUESS/3))=0.0d0
        !write(*,*) "case 2",i,LATSGP
     end do
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.75d0
        !write(*,*) "case 3",i,LATSGP
     end do
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case (180)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     end do
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     end do
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
  case (181)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     end do
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     end do
     do i=6*int(NGUESS/7)+1,NGUESS   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
  case (182)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.75d0
     do i=5,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=2*(RED_POS(1,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/5)+1,5*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case (183)
     do i=1,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     end do
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case (184)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (185)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (186)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (187)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=2/3
     RED_POS(2,5)=1/3
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=2/3
     RED_POS(2,6)=1/3
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/9)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/9)+1,2*int(NGUESS/9)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/9)+1,3*int(NGUESS/9)
        RED_POS(1,i)=2/3
        RED_POS(2,i)=1/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/9)+1,4*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/9)+1,5*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     end do
     do i=5*int(NGUESS/9)+1,6*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/9)+1,7*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/9)+1,8*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     end do
     do i=8*int(NGUESS/9)+1,NGUESS    !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
  case (188)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=2/3
     RED_POS(2,5)=1/3
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=2/3
     RED_POS(2,6)=1/3
     RED_POS(3,6)=0.25d0
     do i=7,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=2/3
        RED_POS(2,i)=1/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-1*(RED_POS(1,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     end do
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case (189)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=1/3
     RED_POS(2,3)=2/3
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=1/3
     RED_POS(2,4)=2/3
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     end do
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     end do
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
        RED_POS(1,i)=1/3
        RED_POS(2,i)=2/3
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
  end select
end subroutine random_incell_p5
