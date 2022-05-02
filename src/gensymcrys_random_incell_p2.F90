subroutine random_incell_p2(LATSGP,CRYSSYS,NGUESS,RED_POS)
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS, I
real(8):: RED_POS(3,NGUESS)
  select case(LATSGP)
  case(50) !choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0
     RED_POS(2,6)=0.25d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/7)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)!
        call random_number(RED_POS(1,i))
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
!  case(50) !choice 2
!     RED_POS(1,1)=0.25d0
!     RED_POS(2,1)=0.25d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.75d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.0d0
!     RED_POS(1,3)=0.75d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,3)=0.5d0
!     RED_POS(1,4)=0.25d0
!     RED_POS(2,4)=0.25d0
!     RED_POS(3,4)=0.5d0
!     RED_POS(1,5)=0.0d0
!     RED_POS(2,5)=0.0d0
!     RED_POS(3,5)=0.0d0
!     RED_POS(1,6)=0.0d0
!     RED_POS(2,6)=0.0d0
!     RED_POS(3,6)=0.5d0
!     do i=7,1*int(NGUESS/7)!
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.75d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)!
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)!
!        RED_POS(1,i)=0.25d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)!
!        RED_POS(1,i)=0.25d0
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
  case(51)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/8)!
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)!
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i)) 
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
!     do i=6*int(NGUESS/9)+1,7*int(NGUESS/9)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)!
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
!     do i=5,1*int(NGUESS/9)!
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/9)+1,2*int(NGUESS/9)!
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/9)+1,3*int(NGUESS/9)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i)) 
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/9)+1,4*int(NGUESS/9)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/9)+1,5*int(NGUESS/9)!
!        RED_POS(1,i)=0.0d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/9)+1,6*int(NGUESS/9)!
!        RED_POS(1,i)=0.0d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/9)+1,7*int(NGUESS/9)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
!     do i=7*int(NGUESS/9)+1,8*int(NGUESS/9)!
!        RED_POS(1,i)=0.25d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 9",i,LATSGP
!     enddo
!     do i=8*int(NGUESS/9)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case10",i,LATSGP
!     enddo
  case(52)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/3)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(53)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
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
  case(54)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     do i=3,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(55)
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
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(56)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/3)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.75d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(57)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     do i=3,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(58)
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
     do i=5,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(59)!choice 1
     RED_POS(1,1)=0.25d0
     RED_POS(2,1)=0.25d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.25d0
     RED_POS(2,2)=0.25d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/5)
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
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
!  case(59)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.0d0
!     RED_POS(3,2)=0.5d0
!     do i=3,1*int(NGUESS/5)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.75d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
!        RED_POS(1,i)=0.25d0
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/5)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
  case(60)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     do i=3,int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(61)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case(62)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/2)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(63)
     RED_POS(1,1)=0.0d0!6   !BUG?
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.0d0
     do i=4,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(64)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.0d0
     do i=4,1*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(65)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0
     RED_POS(2,6)=0.25d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/12)
        RED_POS(1,i)=0.0d0!0.5d0
        RED_POS(2,i)=0.0d0!0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/12)+1,2*int(NGUESS/12)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/12)+1,3*int(NGUESS/12)!
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo

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
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/12)+1,7*int(NGUESS/12)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/12)+1,8*int(NGUESS/12)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/12)+1,9*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/12)+1,10*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case11",i,LATSGP
     enddo
     do i=10*int(NGUESS/12)+1,11*int(NGUESS/12)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case12",i,LATSGP
     enddo
     do i=11*int(NGUESS/12)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case13",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(66)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.25d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0
     RED_POS(2,6)=0.75d0
     RED_POS(3,6)=0.0d0
     do i=7,1*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)!
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo

     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo

     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)!
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
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(67)
     RED_POS(1,1)=0.25d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.25d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0
     RED_POS(2,6)=0.25d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/9)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/9)+1,2*int(NGUESS/9)!
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/9)+1,3*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     enddo

     do i=3*int(NGUESS/9)+1,4*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/9)+1,5*int(NGUESS/9)!
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/9)+1,6*int(NGUESS/9)!
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/9)+1,7*int(NGUESS/9)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/9)+1,8*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/9)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(68)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.25d0
     RED_POS(3,4)=0.25d0
     do i=5,1*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
!  case(68)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.25d0
!     RED_POS(3,1)=0.25d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.75d0
!     RED_POS(1,3)=0.25d0
!     RED_POS(2,3)=0.75d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.0d0
!     RED_POS(2,4)=0.0d0
!     RED_POS(3,4)=0.0d0
!     do i=5,1*int(NGUESS/5)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
!        RED_POS(1,i)=0.0d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/5)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case (69)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.25d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0
     RED_POS(2,6)=0.25d0
     RED_POS(3,6)=0.25d0
     do i=7,1*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/10)+1,2*int(NGUESS/10)!
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/10)+1,3*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo

     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)!
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/10)+1,9*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/10)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
     !write(*,*) "BRAV", 4 ! unit2prim(RED_POS,NGUESS,2,RED_POS)
  case(70)!choice 1
     RED_POS(1,1)=0.125d0
     RED_POS(2,1)=0.125d0
     RED_POS(3,1)=0.125d0
     RED_POS(1,2)=0.625d0
     RED_POS(2,2)=0.625d0
     RED_POS(3,2)=0.625d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
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
     !write(*,*) "BRAV", 4 ! unit2prim(RED_POS,NGUESS,2,RED_POS)
!  case(70)!choice 2
!     RED_POS(1,1)=0.125d0
!     RED_POS(2,1)=0.125d0
!     RED_POS(3,1)=0.125d0
!     RED_POS(1,2)=0.125d0
!     RED_POS(2,2)=0.125d0
!     RED_POS(3,2)=0.875d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.0d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.5d0
!     RED_POS(2,4)=0.5d0
!     RED_POS(3,4)=0.5d0
!     do i=5,1*int(NGUESS/4)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.125d0
!        RED_POS(3,i)=0.125d0
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
!        RED_POS(1,i)=0.125d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.125d0
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
!        RED_POS(1,i)=0.125d0
!        RED_POS(2,i)=0.125d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/4)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     !write(*,*) "BRAV", 4 ! unit2prim(RED_POS,NGUESS,2,RED_POS)
  case(71)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0
     do i=6,1*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/10)+1,2*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/10)+1,3*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/10)+1,9*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/10)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(72)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0
     do i=6,1*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(73)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.25d0
     RED_POS(2,2)=0.25d0
     RED_POS(3,2)=0.25d0
     do i=3,1*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(74)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.25d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.25d0
     RED_POS(2,4)=0.25d0
     RED_POS(3,4)=0.75d0
     do i=5,1*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,6*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(75)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)    !BUG?
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(76)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
  case(77)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(78)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
  case(79)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
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
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(80)
     do i=1,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(81)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(82)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.75d0
     do i=5,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS!3*int(NGUESS/3)    !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(83)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.0d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case(84)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0!0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.25d0
     do i=7,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
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
  end select
end subroutine random_incell_p2
