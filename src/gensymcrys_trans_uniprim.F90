!These subroutines contain the matrices to perform the transition from the conventional unit cell 
!to the primitive unit cell, or the other way around. The matrices are from appendix A.5 from
!the manual to CRYSTAL2003, originally from table 5.1 of the international tables for crystallography A, 1987

subroutine conv2prim(choice,transmat)
implicit none
real(8):: transmat(3,3)
integer:: choice

select case(choice)
case(1) !A --> P
transmat(1,:)=(/ 1.d0, 0.d0, 0.d0/)
transmat(2,:)=(/ 0.d0, 1.d0, 1.d0/)
transmat(3,:)=(/ 0.d0,-1.d0, 1.d0/)
case(2) !B --> P
transmat(1,:)=(/ 1.d0, 0.d0,-1.d0/)
transmat(2,:)=(/ 0.d0, 1.d0, 0.d0/)
transmat(3,:)=(/ 1.d0, 0.d0, 1.d0/)
case(3) !C --> P
transmat(1,:)=(/ 1.d0, 1.d0, 0.d0/)
transmat(2,:)=(/-1.d0, 1.d0, 0.d0/)
transmat(3,:)=(/ 0.d0, 0.d0, 1.d0/)
case(4) !F --> P
transmat(1,:)=(/-1.d0, 1.d0, 1.d0/)
transmat(2,:)=(/ 1.d0,-1.d0, 1.d0/)
transmat(3,:)=(/ 1.d0, 1.d0,-1.d0/)
case(5) !I --> P 
transmat(1,:)=(/ 0.d0, 1.d0, 1.d0/)
transmat(2,:)=(/ 1.d0, 0.d0, 1.d0/)
transmat(3,:)=(/ 1.d0, 1.d0, 0.d0/)
case(6) !H --> R
transmat(1,:)=(/ 1.d0, 0.d0, 1.d0/)
transmat(2,:)=(/-1.d0, 1.d0, 1.d0/)
transmat(3,:)=(/ 0.d0,-1.d0, 1.d0/)
case(7) !P --> P
transmat(1,:)=(/ 1.d0, 0.d0, 0.d0/)
transmat(2,:)=(/ 0.d0, 1.d0, 0.d0/)
transmat(3,:)=(/ 0.d0, 0.d0, 1.d0/)
case default

stop "Wrong choice in conv2prim"
end select
end subroutine

subroutine prim2conv(choice,transmat)
implicit none
real(8):: transmat(3,3)
integer:: choice

select case(choice)
case(1) !P --> A
transmat(1,:)=(/ 1.0d0, 0.0d0, 0.0d0/)
transmat(2,:)=(/ 0.0d0, 0.5d0,-0.5d0/)
transmat(3,:)=(/ 0.0d0, 0.5d0, 0.5d0/)
case(2) !P --> B
transmat(1,:)=(/ 0.5d0, 0.0d0, 0.5d0/)
transmat(2,:)=(/ 0.0d0, 1.0d0, 0.0d0/)
transmat(3,:)=(/-0.5d0, 0.0d0, 0.5d0/)
case(3) !P --> C
transmat(1,:)=(/ 0.5d0,-0.5d0, 0.0d0/)
transmat(2,:)=(/ 0.5d0, 0.5d0, 0.0d0/)
transmat(3,:)=(/ 0.0d0, 0.0d0, 1.0d0/)
case(4) !P --> F
transmat(1,:)=(/ 0.0d0, 0.5d0, 0.5d0/)
transmat(2,:)=(/ 0.5d0, 0.0d0, 0.5d0/)
transmat(3,:)=(/ 0.5d0, 0.5d0, 0.0d0/)
case(5) !P --> I 
transmat(1,:)=(/-0.5d0, 0.5d0, 0.5d0/)
transmat(2,:)=(/ 0.5d0,-0.5d0, 0.5d0/)
transmat(3,:)=(/ 0.5d0, 0.5d0,-0.5d0/)
case(6) !R --> H
transmat(1,:)=(/ 2.d0/3.d0,-1.d0/3.d0,-1.d0/3.d0/)
transmat(2,:)=(/ 1.d0/3.d0, 1.d0/3.d0,-2.d0/3.d0/)
transmat(3,:)=(/ 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0/)
case(7) !P --> P
transmat(1,:)=(/ 1.d0, 0.d0, 0.d0/)
transmat(2,:)=(/ 0.d0, 1.d0, 0.d0/)
transmat(3,:)=(/ 0.d0, 0.d0, 1.d0/)
case default
stop "Wrong choice in prim2conv"
end select
end subroutine
