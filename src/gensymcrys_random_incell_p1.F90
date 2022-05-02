subroutine random_incell_p1(LATSGP,CRYSSYS,NGUESS,RED_POS)
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS, I
real(8):: RED_POS(3,NGUESS)
  select case(LATSGP)
  case (1)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case (2)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.5d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     RED_POS(1,7)=0.0d0
     RED_POS(2,7)=0.5d0
     RED_POS(3,7)=0.5d0
     RED_POS(1,8)=0.5d0
     RED_POS(2,8)=0.5d0
     RED_POS(3,8)=0.5d0
        !write(*,*) "case 1",i,LATSGP
     do i=9,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case(3) !Unique axis b
     do i=1,int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        RED_POS(1,i)=0.5d0
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
  case(4)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
  case(5) !Unique axis b
     do i=1,int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
!     do i=1,int(NGUESS/3)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=int(NGUESS/3)+1,2*int(NGUESS/3)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 3",i,LATSGP
!     enddo
     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     !write(*,*) "ATTENTION BRAV!!!",4    !unit2prim(RED_POS,NGUESS,4,RED_POS)
  case(6) !Unique axis b
     do i=1,int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/3)+1,2*int(NGUESS/3)
        call random_number(RED_POS(1,i))
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
  case(7)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
  case(8) !Unique axis b
     do i=1,int(NGUESS/2)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     !write(*,*) "ATTENTION BRAV!!!",4    !unit2prim(RED_POS,NGUESS,4,RED_POS)
  case(9)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
     !write(*,*) "ATTENTION BRAV!!!",4    !unit2prim(RED_POS,NGUESS,4,RED_POS)
  case (10) !Unique axis b
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0!0.5d0
     RED_POS(3,3)=0.5d0!0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.0d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.5d0
     RED_POS(1,7)=0.5d0
     RED_POS(2,7)=0.0d0
     RED_POS(3,7)=0.5d0
     RED_POS(1,8)=0.5d0
     RED_POS(2,8)=0.5d0
     RED_POS(3,8)=0.5d0
     do i=9,int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
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
  case (11) !Unique axis b
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/2)
        call random_number(RED_POS(1,i))      
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case (12) !Unique axis b
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0!-0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0!-0.5d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.5d0!0.0d0
     RED_POS(1,5)=0.25d0!-0.0d0
     RED_POS(2,5)=0.25d0!0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.25d0!-0.0d0
     RED_POS(2,6)=0.25d0!0.5d0
     RED_POS(3,6)=0.5d0
     do i=7,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
!     do i=7,1*int(NGUESS/4)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.0d0
!        !write(*,*) "case 2",i,LATSGP
!     enddo
!     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=-RED_POS(1,i)
!        RED_POS(3,i)=0.5d0
!        !write(*,*) "case 3",i,LATSGP
!     enddo
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
     !write(*,*) "ATTENTION BRAV!!!",4    !unit2prim(RED_POS,NGUESS,4,RED_POS)
  case (13) !Unique axis b
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.5d0
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
  case (14) !Unique axis b
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
        !write(*,*) "case 1",i,LATSGP
     do i=5,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case (15) !Unique axis b
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0!-0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.25d0!0.0d0
     RED_POS(2,3)=0.25d0!0.5d0
     RED_POS(3,3)=0.0d0!0.0d0
     RED_POS(1,4)=0.25d0!0.0d0
     RED_POS(2,4)=0.25d0!0.5d0
     RED_POS(3,4)=0.5d0!0.5d0
     do i=5,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0 !call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))!RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     !write(*,*) "ATTENTION BRAV!!!",4    !unit2prim(RED_POS,NGUESS,4,RED_POS)
  case (16)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.5d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.5d0
     RED_POS(1,7)=0.0d0
     RED_POS(2,7)=0.5d0
     RED_POS(3,7)=0.5d0
     RED_POS(1,8)=0.5d0
     RED_POS(2,8)=0.5d0
     RED_POS(3,8)=0.5d0
     do i=9,1*int(NGUESS/13)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/13)+1,2*int(NGUESS/13)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/13)+1,3*int(NGUESS/13)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/13)+1,4*int(NGUESS/13)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/13)+1,5*int(NGUESS/13)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/13)+1,6*int(NGUESS/13)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/13)+1,7*int(NGUESS/13)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/13)+1,8*int(NGUESS/13)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/13)+1,9*int(NGUESS/13)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/13)+1,10*int(NGUESS/13)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case11",i,LATSGP
     enddo
     do i=10*int(NGUESS/13)+1,11*int(NGUESS/13)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case12",i,LATSGP
     enddo
     do i=11*int(NGUESS/13)+1,12*int(NGUESS/13)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case13",i,LATSGP
     enddo
     do i=12*int(NGUESS/13)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case14",i,LATSGP
     enddo

  case (17)
     do i=1,1*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        RED_POS(1,i)=0.5d0
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

  case(18)
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
     do i=2*int(NGUESS/3)+1,NGUESS       !BUG?
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo

  case(19)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
  case(20)
     do i=1,1*int(NGUESS/3)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0 !RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.0d0!call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i)) !RED_POS(2,i)=-RED_POS(1,i)
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/3)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(21)

     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.0d0!0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.5d0
     RED_POS(2,3)=0.0d0!-0.5d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     do i=5,1*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0!-RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0!-RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo

     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        RED_POS(1,i)=0.0d0!call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))!RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        RED_POS(1,i)=0.0d0 !call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))!RED_POS(2,i)=RED_POS(1,i)
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 5",i,LATSGP
     enddo

     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0 !RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        RED_POS(1,i)=0.0d0 !0.5d0
        RED_POS(2,i)=0.5d0 !RED_POS(1,i)
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case (22)
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
     do i=5,1*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        RED_POS(1,i)=0.25d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.25d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/7)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     !write(*,*) "BRAV", 4 ! unit2prim(RED_POS,NGUESS,2,RED_POS)
  case(23)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     do i=5,1*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)!+1   !BUG?
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
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(24)
     do i=1,1*int(NGUESS/4)
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
  case (25)
     do i=1,1*int(NGUESS/9)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/9)+1,2*int(NGUESS/9)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/9)+1,3*int(NGUESS/9)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/9)+1,4*int(NGUESS/9)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/9)+1,5*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/9)+1,6*int(NGUESS/9)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/9)+1,7*int(NGUESS/9)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/9)+1,8*int(NGUESS/9)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/9)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case10",i,LATSGP
     enddo
  case(26)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.5d0
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
  case(27)
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
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/5)+1,4*int(NGUESS/5)
        RED_POS(1,i)=0.5d0
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
  case(28)
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
        RED_POS(1,i)=0.25d0
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
  case(29)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
  case(30)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.5d0
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
  case(31)
     do i=1,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
  case(32)
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
  case(33)
     do i=1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 1",i,LATSGP
     enddo
  case(34)
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
  case(35)
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
        RED_POS(1,i)=0.25d0
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
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(36)
     do i=1,1*int(NGUESS/2)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/2)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(37)
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
        RED_POS(1,i)=0.25d0
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
     !write(*,*) "BRAV", 3 ! unit2prim(RED_POS,NGUESS,5,RED_POS)
  case(38)
     do i=1,1*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/6)+1,2*int(NGUESS/6)
        RED_POS(1,i)=0.5d0!RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0!RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/6)+1,3*int(NGUESS/6)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/6)+1,4*int(NGUESS/6)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
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
     !write(*,*) "ATTENTION BRAV!!!",3    !unit2prim(RED_POS,NGUESS,3,RED_POS)
  case(39)
     do i=1,1*int(NGUESS/4)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/4)+1,2*int(NGUESS/4)
        RED_POS(2,i)=0.5d0!RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0!RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/4)+1,3*int(NGUESS/4)
        call random_number(RED_POS(1,i))
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
     !write(*,*) "ATTENTION BRAV!!!",3    !unit2prim(RED_POS,NGUESS,3,RED_POS)
  case(40)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.25d0
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
     !write(*,*) "ATTENTION BRAV!!!",3    !unit2prim(RED_POS,NGUESS,3,RED_POS)
  case(41)
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
     !write(*,*) "ATTENTION BRAV!!!",3    !unit2prim(RED_POS,NGUESS,3,RED_POS)
  case(42)
     do i=1,1*int(NGUESS/5)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/5)+1,2*int(NGUESS/5)
        RED_POS(1,i)=0.25d0
        RED_POS(2,i)=0.25d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/5)+1,3*int(NGUESS/5)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
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
     !write(*,*) "BRAV", 4 ! unit2prim(RED_POS,NGUESS,2,RED_POS)
  case(43)
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
     !write(*,*) "BRAV", 4 ! unit2prim(RED_POS,NGUESS,2,RED_POS)
  case(44)
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
        RED_POS(2,i)=0.0d0!0.25d0
        call random_number(RED_POS(3,i))
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
     !write(*,*) "BRAV", 5 ! unit2prim(RED_POS,NGUESS,1,RED_POS)
  case(45)
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
  case(46)
     do i=1,1*int(NGUESS/3)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/3)+1,2*int(NGUESS/3)
        RED_POS(1,i)=0.25d0
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
  case(47)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.5d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.5d0
     RED_POS(3,5)=0.0d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.5d0
     RED_POS(3,6)=0.0d0
     RED_POS(1,7)=0.0d0
     RED_POS(2,7)=0.5d0
     RED_POS(3,7)=0.5d0
     RED_POS(1,8)=0.5d0
     RED_POS(2,8)=0.5d0
     RED_POS(3,8)=0.5d0
     do i=9,1*int(NGUESS/19)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/19)+1,2*int(NGUESS/19)
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/19)+1,3*int(NGUESS/19)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/19)+1,4*int(NGUESS/19)
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/19)+1,5*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/19)+1,6*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/19)+1,7*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/19)+1,8*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 9",i,LATSGP
     enddo
     do i=8*int(NGUESS/19)+1,9*int(NGUESS/19)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case10",i,LATSGP
     enddo
     do i=9*int(NGUESS/19)+1,10*int(NGUESS/19)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case11",i,LATSGP
     enddo
     do i=10*int(NGUESS/19)+1,11*int(NGUESS/19)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case12",i,LATSGP
     enddo
     do i=11*int(NGUESS/19)+1,12*int(NGUESS/19)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case13",i,LATSGP
     enddo
     do i=12*int(NGUESS/19)+1,13*int(NGUESS/19)
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case14",i,LATSGP
     enddo
     do i=13*int(NGUESS/19)+1,14*int(NGUESS/19)
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case15",i,LATSGP
     enddo
     do i=14*int(NGUESS/19)+1,15*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case16",i,LATSGP
     enddo
     do i=15*int(NGUESS/19)+1,16*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case17",i,LATSGP
     enddo
     do i=16*int(NGUESS/19)+1,17*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case18",i,LATSGP
     enddo
     do i=17*int(NGUESS/19)+1,18*int(NGUESS/19)
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.5d0
        !write(*,*) "case19",i,LATSGP
     enddo
     do i=18*int(NGUESS/19)+1,NGUESS
        call random_number(RED_POS(1,i))
        call random_number(RED_POS(2,i))
        call random_number(RED_POS(3,i))
        !write(*,*) "case20",i,LATSGP
     enddo
  case (48)!choice 1
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.0d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.0d0
     RED_POS(3,3)=0.5d0
     RED_POS(1,4)=0.0d0
     RED_POS(2,4)=0.5d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.25d0
     RED_POS(2,5)=0.25d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.75d0
     RED_POS(2,6)=0.75d0
     RED_POS(3,6)=0.75d0
     do i=7,1*int(NGUESS/7)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/7)+1,2*int(NGUESS/7)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.5d0
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/7)+1,3*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)!
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.0d0
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
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
!  case (48)!choice 2
!     RED_POS(1,1)=0.0d0
!     RED_POS(2,1)=0.0d0
!     RED_POS(3,1)=0.0d0
!     RED_POS(1,2)=0.5d0
!     RED_POS(2,2)=0.5d0
!     RED_POS(3,2)=0.5d0
!     RED_POS(1,3)=0.25d0
!     RED_POS(2,3)=0.75d0
!     RED_POS(3,3)=0.25d0
!     RED_POS(1,4)=0.25d0
!     RED_POS(2,4)=0.25d0
!     RED_POS(3,4)=0.75d0
!     RED_POS(1,5)=0.75d0
!     RED_POS(2,5)=0.25d0
!     RED_POS(3,5)=0.25d0
!     RED_POS(1,6)=0.25d0
!     RED_POS(2,6)=0.25d0
!     RED_POS(3,6)=0.25d0
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
!        RED_POS(3,i)=0.75d0
!        !write(*,*) "case 4",i,LATSGP
!     enddo
!     do i=3*int(NGUESS/7)+1,4*int(NGUESS/7)!
!        call random_number(RED_POS(1,i))
!        RED_POS(2,i)=0.25d0
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 5",i,LATSGP
!     enddo
!     do i=4*int(NGUESS/7)+1,5*int(NGUESS/7)!
!        RED_POS(1,i)=0.25d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.25d0
!        !write(*,*) "case 6",i,LATSGP
!     enddo
!     do i=5*int(NGUESS/7)+1,6*int(NGUESS/7)
!        RED_POS(1,i)=0.75d0!0.25d0
!        call random_number(RED_POS(2,i))
!        RED_POS(3,i)=0.25d0!0.75d0
!        !write(*,*) "case 7",i,LATSGP
!     enddo
!     do i=6*int(NGUESS/7)+1,NGUESS
!        call random_number(RED_POS(1,i))
!        call random_number(RED_POS(2,i))
!        call random_number(RED_POS(3,i))
!        !write(*,*) "case 8",i,LATSGP
!     enddo
  case(49)
     RED_POS(1,1)=0.0d0
     RED_POS(2,1)=0.0d0
     RED_POS(3,1)=0.0d0
     RED_POS(1,2)=0.5d0
     RED_POS(2,2)=0.5d0
     RED_POS(3,2)=0.0d0
     RED_POS(1,3)=0.0d0
     RED_POS(2,3)=0.5d0
     RED_POS(3,3)=0.0d0
     RED_POS(1,4)=0.5d0
     RED_POS(2,4)=0.0d0
     RED_POS(3,4)=0.0d0
     RED_POS(1,5)=0.0d0
     RED_POS(2,5)=0.0d0
     RED_POS(3,5)=0.25d0
     RED_POS(1,6)=0.5d0
     RED_POS(2,6)=0.0d0
     RED_POS(3,6)=0.25d0
     RED_POS(1,7)=0.0d0
     RED_POS(2,7)=0.5d0
     RED_POS(3,7)=0.25d0
     RED_POS(1,8)=0.5d0
     RED_POS(2,8)=0.5d0
     RED_POS(3,8)=0.25d0
     do i=9,1*int(NGUESS/10)
        RED_POS(1,i)=0.0d0!0.25d0
        RED_POS(2,i)=0.0d0!0.75d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 2",i,LATSGP
     enddo
     do i=1*int(NGUESS/10)+1,2*int(NGUESS/10)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 3",i,LATSGP
     enddo
     do i=2*int(NGUESS/10)+1,3*int(NGUESS/10)!
        RED_POS(1,i)=0.5d0
        RED_POS(2,i)=0.0d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 4",i,LATSGP
     enddo
     do i=3*int(NGUESS/10)+1,4*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        RED_POS(2,i)=0.5d0
        call random_number(RED_POS(3,i))
        !write(*,*) "case 5",i,LATSGP
     enddo
     do i=4*int(NGUESS/10)+1,5*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.0d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 6",i,LATSGP
     enddo
     do i=5*int(NGUESS/10)+1,6*int(NGUESS/10)!
        call random_number(RED_POS(1,i))
        RED_POS(2,i)=0.5d0
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 7",i,LATSGP
     enddo
     do i=6*int(NGUESS/10)+1,7*int(NGUESS/10)!
        RED_POS(1,i)=0.5d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
        !write(*,*) "case 8",i,LATSGP
     enddo
     do i=7*int(NGUESS/10)+1,8*int(NGUESS/10)!
        RED_POS(1,i)=0.0d0
        call random_number(RED_POS(2,i))
        RED_POS(3,i)=0.25d0
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
  end select
end subroutine random_incell_p1
