subroutine random_incell_p3(LATSGP,CRYSSYS,NGUESS,RED_POS)
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS, I
real(8):: RED_POS(3,NGUESS)
  select case(LATSGP)
  case(85)!choice 1
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
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
!  case(85)!choice 2
!     RED_POS(1,1)=0.25d0
!     RED_POS(2,1)=0.75d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.25d0
!     RED_POS(2,2)=0.75d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.0d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.0d0
!     RED_POS(2,4)=0.0d0
!     RED_POS(3,4)=0.5d0
!     do i=5,1*int(NGUESS/3)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.75d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/3)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
  case(86)!choice 1
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
     do i=5,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
!  case(86)!choice 2
!     RED_POS(1,1)=0.25d0
!     RED_POS(2,1)=0.25d0
!     RED_POS(3,1)=0.25d0
!     RED_POS(1,2)=0.25d0
!     RED_POS(2,2)=0.25d0
!     RED_POS(3,2)=0.75d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.0d0
!     RED_POS(3,3)=0.0d0
!     RED_POS(1,4)=0.0d0
!     RED_POS(2,4)=0.0d0
!     RED_POS(3,4)=0.5d0
!     do i=5,1*int(NGUESS/3)
!        RED_POS(1,i)=0.75d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
!        RED_POS(1,i)=0.25d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/3)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
  case(87)
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
     do i=6,1*int(NGUESS/4)
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
        call random_number(RED_POS(1,i))!RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0 !call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(88)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.25d0
     RED_POS(3,3)=0.125d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.25d0
     RED_POS(3,4)=0.625d0
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
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
!  case(88)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.0d0
!     RED_POS(2,2)=0.0d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.0d0
!     RED_POS(2,3)=0.25d0
!     RED_POS(3,4)=5.0d0/8.0d0
!     RED_POS(1,5)=0.0d0
!     RED_POS(2,5)=0.25d0
!     RED_POS(3,5)=0.125d0
!     do i=6,1*int(NGUESS/2)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.25d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/2)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(89)
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
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/10)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/10)+1,2*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/10)+1,3*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/10)+1,9*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/10)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
!     do i=7,1*int(NGUESS/8)!
!        RED_POS(1,i)=0.5d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)!
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)!
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.5d0
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 8",i,LATSGP
!     enddo
!     do i=7*int(NGUESS/8)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 9",i,LATSGP
!     enddo
  case(90)
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
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(91)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=3.0d0/8.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(92)
     do i=1,1*int(NGUESS/2)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(93)
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
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.0d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.25d0
     do i=7,1*int(NGUESS/10)
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
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.75d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/10)+1,9*int(NGUESS/10)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/10)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
  case(94)
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
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
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
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(95)
     do i=1,1*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=5.0d0/8.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
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
  case(96)
     do i=1,1*int(NGUESS/2)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(97)
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
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(98)
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
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,2)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,2)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,2)=0.125d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,5*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(99)
     do i=1,1*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i))
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
  case(100)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0 !RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(101)
     do i=1,1*int(NGUESS/5)
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
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(102)
     do i=1,1*int(NGUESS/4)
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
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case(103)
     do i=1,1*int(NGUESS/4)
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
        RED_POS(1,i)=0.5d0
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
  case(104)
     do i=1,1*int(NGUESS/3)!5   !BUG?  
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)!6,10  !BUG?
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS!11,40   !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
  case(105)
     do i=1,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
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
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case(106)
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
  case(107)
     do i=1,1*int(NGUESS/5)
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
        RED_POS(2,i)=RED_POS(1,i)
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
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(108)
     do i=1,1*int(NGUESS/4)
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
        RED_POS(2,i)=RED_POS(1,i)+0.5d0 !RED_POS(1,i)
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
  case(109)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
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
  case(110)

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
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(111)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/9)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/9)+1,2*int(NGUESS/9)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/9)+1,3*int(NGUESS/9)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/9)+1,4*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/9)+1,5*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/9)+1,6*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/9)+1,7*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/9)+1,8*int(NGUESS/9)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/9)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
  case(112)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.25d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.25d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.0d0

     do i=7,1*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/8)+1,2*int(NGUESS/8)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/8)+1,3*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/8)+1,4*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)!
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
  case(113)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/4)
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
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/4)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
  case (114)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/3)
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
  case(115)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/8)
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
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo

     do i=4*int(NGUESS/8)+1,5*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/8)+1,6*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/8)+1,7*int(NGUESS/8)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/8)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
  case(116)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.25d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/6)
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
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo

     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.75d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
  case(117)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.5d0!0.0d0
     RED_POS(3,1)=0.0d0!0.25d0
     RED_POS(1,2)=0.0d0!0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.5d0!0.25d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0!0.5d0
     RED_POS(2,4)=0.0d0!0.5d0
     RED_POS(3,4)=0.5d0!0.0d0
     do i=5,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0!0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0!RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     enddo

     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0!RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(118)
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
        RED_POS(2,i)=RED_POS(1,i)+0.5d0!RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo

     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=-RED_POS(1,i)+0.5d0!-RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
  case(119)
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
     do i=5,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)!+0.5d0
        RED_POS(3,i)=0.0d0!call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)+0.5d0
        RED_POS(3,i)=0.25d0
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
!     do i=5,1*int(NGUESS/5)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
!        RED_POS(1,i)=0.0d0
!        RED_POS(2,i)=0.5d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 3",i,LATSGP
!     enddo
!     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=RED_POS(1,i)!+0.5d0
!        RED_POS(3,i)=0.0d0!call random_number(RED_POS(3,i))
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.0d0
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/5)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 6",i,LATSGP
!     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(120)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.25d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.25d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
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
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/5)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(121)
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
     do i=5,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/6)+1,5*int(NGUESS/6)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/6)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(122)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     do i=3,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,i)=0.125d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(123)
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
     RED_POS(3,5)=0.5d0
     RED_POS(1,6)=0.0d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.0d0
     do i=7,1*int(NGUESS/15)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/15)+1,2*int(NGUESS/15)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/15)+1,3*int(NGUESS/15)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/15)+1,4*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/15)+1,5*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/15)+1,6*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/15)+1,7*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/15)+1,8*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/15)+1,9*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/15)+1,10*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case11",i,LATSGP
     enddo
     do i=10*int(NGUESS/15)+1,11*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case12",i,LATSGP
     enddo
     do i=11*int(NGUESS/15)+1,12*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case13",i,LATSGP
     enddo
     do i=12*int(NGUESS/15)+1,13*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case14",i,LATSGP
     enddo
     do i=13*int(NGUESS/15)+1,14*int(NGUESS/15)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case15",i,LATSGP
     enddo
     do i=14*int(NGUESS/15)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case16",i,LATSGP
     enddo
  end select
end subroutine random_incell_p3
