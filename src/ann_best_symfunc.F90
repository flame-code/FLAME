!*****************************************************************************************
subroutine ann_best_symfunc(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann):: ann, ann_t
    type(typ_atoms_arr):: atoms_train
    type(typ_symfunc_arr):: symfunc_train
    integer:: n_tot, i_rc, i_eta, i_zeta, i_lambda
    real(8), allocatable:: his(:,:)
    real(8):: disparity
    logical:: file_exists
    stop 'ERROR: routine not ready since ann -> ann_arr, take care of lines with HERE'
    !---------------------------------------------
    !call read_input_ann(parini,iproc,ann)
    inquire(file="list_posinp_train.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_train.yaml',atoms_train)
    else
        call read_data_old(parini,'list_posinp_train',atoms_train)
    endif
    if(iproc==0) then
        write(*,'(a34,i8)') 'number of training data points:   ',atoms_train%nconf
    endif

    allocate(his(1000,10000),source=0.d0)
    ann_t%ng1=0
    ann_t%ng2=0
    ann_t%ng3=0
    ann_t%ng4=0
    ann_t%ng5=0
    ann_t%ng6=0
    n_tot=0
    ann_t%nn(0)=1
    !-------------------------------------------------------
    !varying parameters in symmetry function G2
    ann_t%g2rs(1)=0.d0
    !ann_t%nn(0)=ann_t%ng1+ann_t%ng2+ann_t%ng3+ann_t%ng4+ann_t%ng5+ann_t%ng6
    ann_t%ng2=1
    do i_rc=1,5
!HERE        ann_t%g2rc(1)=8.d0+(i_rc-1)*0.5d0
        do i_eta=1,101
            ann_t%g2eta(1)=0.001d0+(i_eta-1)*2.d-2
!HERE            call set_gbounds(ann_t,atoms_train,'bounds_train',symfunc_train)
            call gbounds_distro(ann_t,atoms_train,'train')
            call cal_symfunc_diversity(n_tot,his,ann_t,disparity)
            !write(61,'(2i5,es14.5)') i_rc,i_eta,disparity
            if(disparity>0.06d0 .or. n_tot==0) then
                n_tot=n_tot+1
                his(1:1000,n_tot)=ann_t%his(1:1000,1)
                write(61,'(2i5,es14.5)') i_rc,i_eta,disparity
!HERE                write(71,'(3f8.4)') ann_t%g2rc(1),ann_t%g2eta(1),ann_t%g2rs(1)
            endif
        enddo
    enddo
    ann_t%ng2=0
    !-------------------------------------------------------
    !varying parameters in symmetry function G5
    ann_t%ng5=1
    do i_rc=1,3
!HERE        ann_t%g5rc(1)=8.d0+(i_rc-1)*1.d0
        do i_eta=1,21
            ann_t%g5eta(1)=0.001d0+(i_eta-1)*2.d-2
            do i_zeta=1,16
                ann_t%g5zeta(1)=real(i_zeta,8)
                do i_lambda=1,2
                    ann_t%g5lambda=real((-1)**i_lambda,8)
!HERE                    call set_gbounds(ann_t,atoms_train,'bounds_train',symfunc_train)
                    call gbounds_distro(ann_t,atoms_train,'train')
                    call cal_symfunc_diversity(n_tot,his,ann_t,disparity)
                    if(disparity>0.06d0 .or. n_tot==0) then
                        n_tot=n_tot+1
                        his(1:1000,n_tot)=ann_t%his(1:1000,1)
                        write(61,'(4i5,es14.5)') i_rc,i_eta,i_zeta,i_lambda,disparity
!HERE                        write(71,'(4f8.4)') ann_t%g5rc(1),ann_t%g5eta(1),ann_t%g5zeta(1),ann_t%g5lambda(1)
                    endif
                enddo
            enddo
        enddo
    enddo
    ann_t%ng5=0
    !-------------------------------------------------------
    deallocate(his)
end subroutine ann_best_symfunc
!*****************************************************************************************
subroutine cal_symfunc_diversity(n_tot,his,ann,disparity)
    use mod_ann, only: typ_ann
    implicit none
    integer, intent(in):: n_tot
    real(8), intent(in):: his(1000,n_tot)
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: disparity
    !local variables 
    integer:: ibin, i
    real(8):: tt
    if(n_tot==0) then
        disparity=0.d0
        return
    endif
    disparity=1.d20
    do i=1,n_tot
        tt=0.d0
        do ibin=1,1000
            tt=tt+(his(ibin,i)-ann%his(ibin,1))**2
        enddo
        tt=sqrt(tt)
        disparity=min(disparity,tt)
        !write(*,'(f5.2)',advance='no') disparity
    enddo
end subroutine cal_symfunc_diversity
!*****************************************************************************************
subroutine gbounds_distro(ann,atoms_arr,strmess)
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann
    use mod_processors, only: iproc
    implicit none
    type(typ_ann), intent(inout):: ann
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    !local variables 
    !type(typ_ann):: ann_t
    integer:: ig, iconf, iat, i0, i, j, ibin
    real(8):: hbin(100), realn, tt
    do i0=1,ann%nn(0)
        write(51,'(i5,2es14.5)') i0,ann%gbounds(1,i0),ann%gbounds(2,i0)
        hbin(i0)=(ann%gbounds(2,i0)-ann%gbounds(1,i0))/1000.d0
    enddo
    !The following values for bounds will keep the symmetry functions intact
    ann%his(1:1000,1:ann%nn(0))=0.d0
    realn=0.d0
    do iconf=1,atoms_arr%nconf
        do iat=1,atoms_arr%atoms(iconf)%nat
            realn=realn+1
            !This line commented when input.ann became type depedent, so it needs to be fixed
            !call symmetry_functions(parini,ann,iat,atoms_arr%atoms(iconf),.false.)
            do i=1,ann%nn(0)
                ibin=int((ann%y(i,0)-ann%gbounds(1,i))/hbin(i))+1
                ibin=min(ibin,1000)
                !if(ibin<1 .or. ibin>1000) then 
                !    write(*,'(a,i7)') 'ERROR: ibin out of bounds',ibin
                !    stop
                !endif
                ann%his(ibin,i)=ann%his(ibin,i)+1.d0
            enddo
        enddo
    enddo
    do i=1,ann%nn(0)
        do ibin=1,1000
            ann%his(ibin,i)=ann%his(ibin,i)/realn
            !write(100+i,'(i6,es14.5)') ibin,ann%his(ibin,i)
        enddo
    enddo
end subroutine gbounds_distro
!*****************************************************************************************
