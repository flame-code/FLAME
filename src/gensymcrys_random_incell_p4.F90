subroutine random_incell_p4(LATSGP,CRYSSYS,NGUESS,RED_POS)
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS, I
real(8):: RED_POS(3,NGUESS)
  select case(LATSGP)
  case(124)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.0d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.25d0
     do i=7,1*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
  case(125)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0
     RED_POS(2,6)=0.25d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)+0.5d0
        call random_number(RED_POS(3,i)) 
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
!  case(125)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.0d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.75d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=0.5d0
!     RED_POS(1,4)=0.75d0
!     RED_POS(2,4)=0.25d0
!     RED_POS(3,4)=0.0d0
!     RED_POS(1,5)=0.25d0
!     RED_POS(2,5)=0.25d0
!     RED_POS(3,5)=0.5d0
!     RED_POS(1,6)=0.25d0
!     RED_POS(2,6)=0.25d0
!     RED_POS(3,6)=0.0d0
!     do i=7,1*int(NGUESS/8)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        call random_number(RED_POS(3,i)) 
!        !write(*,*) "case 8",i,LATSGP
!     enddo
!     do i=7*int(NGUESS/8)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 9",i,LATSGP
!     enddo
  case(126)!choice 1
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
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0
     do i=6,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
!  case(126)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.25d0
!     RED_POS(2,2)=0.75d0
!     RED_POS(3,2)=0.0d0
!     RED_POS(1,3)=0.25d0
!     RED_POS(2,3)=0.75d0
!     RED_POS(3,3)=0.75d0
!     RED_POS(1,4)=0.25d0
!     RED_POS(2,4)=0.25d0
!     RED_POS(3,4)=0.75d0
!     RED_POS(1,5)=0.25d0
!     RED_POS(2,5)=0.25d0
!     RED_POS(3,5)=0.25d0
!     do i=6,1*int(NGUESS/6)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.75d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.75d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/6)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 7",i,LATSGP
!     enddo
  case(127)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
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
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
!     do i=5,1*int(NGUESS/7)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)+0.5d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)+0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/7)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
  case(128)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0
     do i=5,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.0d0 !STRANGE
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(129)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.25d0
     RED_POS(2,4)=0.25d0
     RED_POS(2,4)=0.5d0
     do i=5,1*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i));RED_POS(1,i)=-RED_POS(1,i)
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
!  case(129)!choice 2
!     RED_POS(1,1)=0.75d0
!     RED_POS(2,1)=0.25d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.75d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.0d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.0d0
!     RED_POS(2,4)=0.0d0
!     RED_POS(2,4)=0.5d0
!     do i=5,1*int(NGUESS/7)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
!        RED_POS(1,i)=0.25d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/7)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
  case(130)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.25d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.0d0
     do i=4,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
!  case(130)
!     RED_POS(1,1)=0.75d0
!     RED_POS(2,1)=0.25d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.75d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.25d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.0d0
!     RED_POS(3,3)=0.0d0
!     do i=4,1*int(NGUESS/4)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/4)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
  case(131)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.25d0
     do i=7,1*int(NGUESS/12)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/12)+1,2*int(NGUESS/12)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/12)+1,3*int(NGUESS/12)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
!     do i=3*int(NGUESS/13)+1,4*int(NGUESS/13)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo

     do i=3*int(NGUESS/12)+1,4*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/12)+1,5*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/12)+1,6*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/12)+1,7*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/12)+1,8*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/12)+1,9*int(NGUESS/12)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/12)+1,10*int(NGUESS/12)!
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
     do i=10*int(NGUESS/12)+1,11*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case12",i,LATSGP
     enddo
     do i=11*int(NGUESS/12)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case13",i,LATSGP
     enddo

!     do i=7,1*int(NGUESS/13)!
!        RED_POS(1,i)=0.5d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/13)+1,2*int(NGUESS/13)!
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/13)+1,3*int(NGUESS/13)!
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/13)+1,4*int(NGUESS/13)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!
!     do i=4*int(NGUESS/13)+1,5*int(NGUESS/13)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/13)+1,6*int(NGUESS/13)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/13)+1,7*int(NGUESS/13)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 8",i,LATSGP
!     enddo
!     do i=7*int(NGUESS/13)+1,8*int(NGUESS/13)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 9",i,LATSGP
!     enddo
!     do i=8*int(NGUESS/13)+1,9*int(NGUESS/13)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case10",i,LATSGP
!     enddo
!     do i=9*int(NGUESS/13)+1,10*int(NGUESS/13)!
!        RED_POS(1,i)=0.0d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case11",i,LATSGP
!     enddo
!     do i=10*int(NGUESS/13)+1,11*int(NGUESS/13)!
!        RED_POS(1,i)=0.5d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case12",i,LATSGP
!     enddo
!     do i=11*int(NGUESS/13)+1,12*int(NGUESS/13)!
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case13",i,LATSGP
!     enddo
!     do i=12*int(NGUESS/13)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case14",i,LATSGP
!     enddo
  case(132)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.0d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.0d0

     do i=7,1*int(NGUESS/10)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/10)+1,2*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/10)+1,3*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo

     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo

     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo

     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 9",i,LATSGP
     enddo

     do i=8*int(NGUESS/10)+1,9*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/10)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
  case(133)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0

     do i=6,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.75d0
        !write(*,*) "case 5",i,LATSGP
     enddo

     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i));RED_POS(1,i)=-RED_POS(1,i)
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo

     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
!  case(133)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.75d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.75d0
!     RED_POS(1,3)=0.25d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=0.25d0
!     RED_POS(1,4)=0.75d0
!     RED_POS(2,4)=0.25d0
!     RED_POS(3,4)=0.0d0
!     RED_POS(1,5)=0.25d0
!     RED_POS(2,5)=0.255d0 !IS IT CORRECT?
!     RED_POS(3,5)=0.0d0
!
!     do i=6,1*int(NGUESS/6)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!
!     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!
!     do i=5*int(NGUESS/6)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 7",i,LATSGP
!     enddo
  case(134)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.75d0
     RED_POS(2,6)=0.75d0
     RED_POS(3,6)=0.75d0

     do i=7,1*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo

     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i));RED_POS(1,i)=-RED_POS(1,i)
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
        call random_number(RED_POS(1,i));RED_POS(1,i)=-RED_POS(1,i)
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.75d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
!  case(134)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.0d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.25d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.25d0
!     RED_POS(2,4)=0.25d0
!     RED_POS(3,4)=0.25d0
!     RED_POS(1,5)=0.75d0
!     RED_POS(2,5)=0.25d0
!     RED_POS(3,5)=0.25d0
!     RED_POS(1,6)=0.25d0
!     RED_POS(2,6)=0.75d0
!     RED_POS(3,6)=0.25d0
!
!     do i=7,1*int(NGUESS/8)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.75d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!
!     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
!     do i=7*int(NGUESS/8)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 9",i,LATSGP
!     enddo
  case(135)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0

     do i=5,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.25d0!0.5d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(136)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0

     do i=5,1*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
  case(137)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.25d0

     do i=4,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
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
!  case(137)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.75d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.25d0
!     RED_POS(1,3)=0.75d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=0.75d0
!
!     do i=4,1*int(NGUESS/5)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!
!     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
!        RED_POS(1,i)=0.25d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/5)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
  case(138)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.25d0
     RED_POS(2,4)=0.25d0
     RED_POS(3,4)=0.75d0

     do i=5,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.75d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i));RED_POS(1,i)=-RED_POS(1,i)
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS !6*int(NGUESS/6)   !¨BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
!  case(138)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.75d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.25d0
!     RED_POS(1,3)=0.75d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.0d0
!     RED_POS(2,4)=0.0d0
!     RED_POS(3,4)=0.5d0
!
!     do i=5,1*int(NGUESS/6)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!
!     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/6)+1,NGUESS !6*int(NGUESS/6)   !¨BUG?
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 7",i,LATSGP
!     enddo
  case(139)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0

     do i=6,1*int(NGUESS/10)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/10)+1,2*int(NGUESS/10)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/10)+1,3*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)
        call random_number(RED_POS(1,i));RED_POS(1,i)=-RED_POS(1,i)
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/10)+1,9*int(NGUESS/10)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/10)+1,NGUESS   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
!     do i=6,1*int(NGUESS/8)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)+0.5d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 8",i,LATSGP
!     enddo
!     do i=7*int(NGUESS/8)+1,NGUESS   !BUG?
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 9",i,LATSGP
!     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
!HERE
  case(140)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0
     do i=6,1*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS  !int(NGUESS/7)   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo

!     do i=6,1*int(NGUESS/7)!
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)!
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)+0.5d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/7)+1,NGUESS  !int(NGUESS/7)   !BUG?
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  end select
end subroutine random_incell_p4
