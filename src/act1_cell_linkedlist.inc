do iz=1,linked_lists%mz+linked_lists%mlimnb3
do iy=1-linked_lists%mlimnb2,linked_lists%my+linked_lists%mlimnb2
do ix=1-((isign(1,iy-1)-1)/2)*linked_lists%mx,linked_lists%mx+linked_lists%mlimnb1
ip=linked_lists%prime(ix,iy,iz)
il=linked_lists%last(ix,iy,iz)
do jz=iz,min(linked_lists%mz+linked_lists%mlimnb3,iz+linked_lists%mlimnb)
do jy=iy+linked_lists%limnby(1,jz-iz),iy+linked_lists%limnby(2,jz-iz)
