subroutine correct_latvec(latvec,pos_red,nat,correctalg,iout)
use yaml_output
implicit none
integer:: correctalg,nat,iproc,iout
real(8):: latvec(3,3),pos_red(3,nat),latvec0(3,3),diff(9)
latvec0=latvec
select case(correctalg)
 case(1)
  !The oganov method
  iproc=1
  call correct_latvec_oganov(latvec,pos_red,nat,iproc)
 case(2)
  !Niggli reduction
  call fixcell_niggli(nat,latvec,pos_red)
 case default
  stop "Currently not implemented method for cell correction!"
end select 
diff(1:3)=latvec(:,1)-latvec0(:,1)
diff(4:6)=latvec(:,2)-latvec0(:,2)
diff(7:9)=latvec(:,3)-latvec0(:,3)
if(dot_product(diff,diff).lt.1.d-14) then
    call yaml_mapping_open('Cell not changed',flow=.true.)
    call yaml_map('latvec',latvec,fmt='(es20.12)')
    call yaml_map('latvec0',latvec0,fmt='(es20.12)')
    call yaml_mapping_close()
 !write(*,*) "Cell not changed"
 !write(*,*) latvec
 !write(*,*) latvec0
 iout=0
else
 iout=1
    call yaml_mapping_open('Cell changed',flow=.true.)
    call yaml_map('latvec',latvec,fmt='(es20.12)')
    call yaml_map('latvec0',latvec0,fmt='(es20.12)')
    call yaml_mapping_close()
 !write(*,*) "Cell changed"
 !write(*,*) latvec
 !write(*,*) latvec0
endif


end subroutine
