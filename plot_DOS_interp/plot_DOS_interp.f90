program plot_DOS_interp
    !******************************************
    !cc     Written by Lin-Wang Wang, March 30, 2001.  
    !cc     Copyright 2001 The Regents of the University of California
    !cc     The United States government retains a royalty free license in this work
    !******************************************

    !****************************************
    !****  It stores the wavefunction in G space, only in half
    !****  of the E_cut sphere. 
    !******************************************
    use pseudo_potential_upf
    use omp_lib

    implicit double precision (a-h,o-z)

    !parameter (matom_1=50000,mtype_1=20,nE=4000)
    parameter (matom_1=50000,mtype_1=20)
    !******************************************
    real*8 AL(3,3),AL_t(3,3)
    !***********************************************
    !*************************************************
    real*8 xatom(3,matom_1),w(matom_1),weight(8,mtype_1), weight_all(8,mtype_1)
    real*8 weight_ll(9,8,mtype_1)
    integer iatom(matom_1),imov_at(3,matom_1),iiatom(mtype_1),&
        &    ityatom(matom_1)
    integer numref(matom_1),iref_start_8(matom_1)
    integer is_ref(mtype_1),ip_ref(mtype_1),id_ref(mtype_1),if_ref(mtype_1)
    character*3 psp_head
    real*8 Z_type(mtype_1),Ecut_type(mtype_1),amass_type(mtype_1)
    real*8 error_ii(8)
    integer ind_ii(0:1,0:1,0:1)

    !ccccccccc not really used  
    integer smatr(3,3,48),nrot
    integer iCGmth0(100),iCGmth1(100),iscfmth0(100),iscfmth1(100)
    real*8 FermidE0(100),FermidE1(100),totNel
    integer itypeFermi0(100),itypeFermi1(100)
    integer kpt_dens(2),ispin_dens(2),iw_dens(2)
    integer niter0,nline0,niter1,nline1
    integer icoul
    real*8 xcoul(3)
    integer ipsp_type(mtype_1),nref_type(mtype_1)
    character*200 vwr_atom(mtype_1),fforce_out,fdens_out
    character*20 fwg_in(2),fwg_out(2),frho_in(2),frho_out(2),&
        & fvr_in(2),fvr_out(2),f_tmp,fxatom_out,fvext_in
    character*20 f_xatom,sym_file,kpt_file
    character(len=200) :: filename
    !cccccccccccccccccccccccccccccc
    real*8, allocatable, dimension(:):: weighkpt_2,akx_2,aky_2,akz_2 
    real*8, allocatable, dimension(:,:,:):: E_st
    real*8, allocatable, dimension(:,:,:,:) ::  weight_ii, weight_ii_all
    real*8, allocatable, dimension(:,:,:,:,:) ::  weight_ll_ii
    real*8, allocatable, dimension(:,:,:) ::  weight_8, weight_8_all
    real*8, allocatable, dimension(:,:,:,:) ::  weight_ll_8

    complex*16, allocatable, dimension(:,:) :: beta_psi_tmp
    complex*16, allocatable, dimension(:,:,:,:) :: beta_psi_8
    complex*16, allocatable, dimension(:,:,:) :: S_m
    complex*16, allocatable, dimension(:,:,:,:,:) :: hh
    complex*16, allocatable, dimension(:) :: work
    real*8, allocatable, dimension(:):: rwork
    real*8, allocatable, dimension(:,:,:,:):: EE
    integer j123(3,8),ispin_ind(8),ikpt_ind(8),ikpt_ext_ind(8)
    integer hse_nq1,hse_nq2,hse_nq3
    complex*16 cc

    integer :: nE
    real*8,allocatable,dimension(:,:) :: DOS
    real*8,allocatable,dimension(:,:,:,:) :: DOS_ref, DOS_ref_all
    real*8,allocatable,dimension(:,:,:,:,:) :: DOS_ref_ll

    !real*8 DOS(nE,2),DOS_ref(nE,8,mtype_1,2), DOS_ref_all(nE,8,mtype_1,2)
    !real*8 DOS_ref_ll(nE,9,8,mtype_1,2)
    integer  lll(8,mtype_1),nchi(mtype_1)
    character*9 fname
    character*20 :: label

    !cccccccccccccccccccccccc
    character*200 message
    character*10 JOB, temp_job

    REAL       RDUM(20)
    INTEGER    IDUM(20)
    CHARACTER(60)  CHARAC
    INTEGER    N1_NUM(20)
    COMPLEX    CDUM(20)
    LOGICAL    LDUM(20),FOUND,CONT,YINTST,LOPEN
    LOGICAL    L
    INTEGER    N,IERR
    character(len=200) temp_char
    CHARACTER(LEN=200) :: right, key_words, right_efec
    LOGICAL :: right_logical
    INTEGER :: flag_dens

    LOGICAL :: readit

    CHARACTER(len=2) :: atom_name
    integer :: ll
    
    integer :: istep
    !
    real*8, allocatable, dimension(:,:,:):: occ
    real*8 :: stime, etime


    istep=0
    !**************************************************
    !**************************************************
    !       write(6,*) "input: 0 (total DOS) or 1 (atom part DOS, need weight in atom.config)"
    !       read(5,*) ipart_DOS
    !       write(6,*) "input: 0 (no-interp) or 1(interp,needs OUT.overlap_uk)"
    !       read(5,*) interp

    write(6,*) "input param from DOS.input"
    open(10,file="DOS.input",status="old",action="read",iostat=ierr)
    if(ierr.ne.0) then
        write(6,*) "The DOS.input file not found, stop"
        stop
    end if
    rewind(10)
    !       write(6,*) "input: 0 (total DOS) or 1 (atom part DOS, need weight in atom.config)"
    read(10,*) ipart_DOS
    !       write(6,*) "input: 0 (no-interp) or 1(interp,needs OUT.overlap_uk)"
    read(10,*) interp
    !read(10,*) E_b   ! energy smearing in eV
    read(10,'(A200)') right
    !
    i_index=index(right,'#')
    if(i_index>0) then
        right=right(1:i_index-1)
    endif
    !
    read(right,*,iostat=ierr) E_b, nE
    if(ierr.ne.0) then
        read(right,*,iostat=ierr) E_b
        nE=4000
    endif
    !
    read(10,*) nm1,nm2,nm3
    !    nm1,nm2,nm3 are the interpolation grid for each k-point in the
    !    MP_N123 grid (within each box, interpolation using 8-k points in
    !    the MP123 grid). Larger nm1,2,3 (e.g., 10,10,10) will generate
    !    smoother DOS, but it will be more time consuming. 8x8x8 might be
    !    fine. It should go with E_b. Larger E_b (e.g., 0.1) requires less
    !    nm1,2,3. Recommand: E_b=0.05, nm123=6,6,6 
    read(10,*,iostat=ierr) iocc
    if(ierr.ne.0) then
        iocc=0
    endif
    close(10)
    
    allocate(DOS(nE,2))
    allocate(DOS_ref(nE,8,mtype_1,2))
    allocate(DOS_ref_all(nE,8,mtype_1,2))
    allocate(DOS_ref_ll(nE,9,8,mtype_1,2))

    open(9,file='REPORT',status='old',action='read',iostat=ierr) 
    if(ierr.ne.0) then
        write(6,*) "file ***REPORT*** does not exist, stop"
        stop
    end if
    rewind(9)

    CALL read_key_words ( 9, 'IN.ATOM',LEN('IN.ATOM'), right, readit )
    if(readit) then
        READ ( right, "(A20)" ) f_xatom
    else
        f_xatom = "atom.config"
    end if

    call readf_xatom_new(w)

    ntype = 0
    DO i = 1, 100
        !temp_char = "IN.PSP"//CHAR(i+48)
        write(temp_char,*) i
        temp_char="IN.PSP"//ADJUSTL(trim(temp_char))
        temp_char = ADJUSTL(temp_char)
        CALL read_key_words ( 9, TRIM(temp_char),len_trim(temp_char) &
            &                           , right, readit )
        IF( readit ) THEN
            ntype = ntype + 1
        ELSE
            EXIT
        ENDIF
    ENDDO


    allocate(upfpsp(ntype))
    !cccccccccccccccccccccccccccccccccccccccccccccccc
    do I = 1,ntype
        !temp_char = "IN.PSP"//char(48+I)
        write(temp_char,*) i
        temp_char="IN.PSP"//ADJUSTL(trim(temp_char))
        temp_char = ADJUSTL(temp_char)
        CALL read_key_words ( 9, TRIM(temp_char),len_trim(temp_char) &
            &                          , right, readit )
        READ ( right, "(A200)" ) vwr_atom(i)
        psp_head = "NUL"

!        call trim_string(vwr_atom(i), psp_head)
!        if(psp_head.eq."usp") then          ! NEED TO BE FIXED LATER
!            ipsp_type(i)=2
!        elseif (psp_head .eq."vwr") then
!            ipsp_type(i)=1
!        else
            !  iflag=1, for DOS, inside read_upf_v2, it will use the wave
            !  function as the projector
            call read_upf_v2(vwr_atom(i), upfpsp(i),1)
            if (upfpsp(i)%tvanp) then
                ipsp_type(i) = 2
            else
                ipsp_type(i) = 1
            end if
            !               message =&
            !     &         "IN.PSP# must be vwr.xxx(norm),or uspp.xxx(ultra-soft)"
            !               call error_stop(message)
!        endif
            !if(inode_tot.eq.1) then
            !    write(22,*) TRIM(temp_char)//" = ", vwr_atom(I)
            !endif
    end do

    ipsp_all=1
    do ia=1,ntype
        if(ipsp_type(ia).eq.2) ipsp_all=2
    end do

    !cccccccccccccccccccccccccccccccccc
    do ia=1,ntype
        if(index(vwr_atom(ia),'vwr')>0) then
            call readvwr_head()
        elseif (index(vwr_atom(ia),'usp')>0) then
            call readusp_head(vwr_atom(ia),iiatom(ia),nref_type(ia),&
                &                          Z_type(ia),Ecut_type(ia),amass_type(ia))
        else
            nref_type(ia) = 0
            iiatom(ia) = upfpsp(ia)%num
            is_ref(ia)=0
            ip_ref(ia)=0
            id_ref(ia)=0
            if_ref(ia)=0
            !               do j = 1, upfpsp(ia)%nbeta

            if(upfpsp(ia)%nwfc.gt.0) then
                !
                nchi(ia)=upfpsp(ia)%nwfc
                nref_type(ia)=0
                do j=1,nchi(ia)
                    lll(j,ia)=upfpsp(ia)%lchi(j)
                    nref_type(ia)=nref_type(ia)+lll(j,ia)*2+1
                enddo


                !do j = 1, upfpsp(ia)%nwfc
                !    nref_type(ia) = nref_type(ia) + upfpsp(ia)%lchi(j)*2+1
                !    if(upfpsp(ia)%lchi(j).eq.0) is_ref(ia)=is_ref(ia)+1
                !    if(upfpsp(ia)%lchi(j).eq.1) ip_ref(ia)=ip_ref(ia)+1
                !    if(upfpsp(ia)%lchi(j).eq.2) id_ref(ia)=id_ref(ia)+1
                !    if(upfpsp(ia)%lchi(j).eq.3) if_ref(ia)=if_ref(ia)+1
                !end do
            else
                !nref_type(ia)=16
                !is_ref(ia)=1
                !ip_ref(ia)=1
                !id_ref(ia)=1
                !if_ref(ia)=1
                nchi(ia)=4
                nref_type(ia)=16
                lll(1,ia)=0 
                lll(2,ia)=1 
                lll(3,ia)=2 
                lll(4,ia)=3 
            endif


            z_type(ia) = upfpsp(ia)%zp
            Ecut_type(ia) =  upfpsp(ia)%ecutwfc
            amass_type(ia) = upfpsp(ia)%mass
        end if
    end do


    nref_tot_8=0
    mref_8=0
    totNel=0.d0
    do ia=1,natom
        iref_start_8(ia)=nref_tot_8
        iitype=0
        do itype=1, ntype
            if(iatom(ia).eq.iiatom(itype)) iitype=itype
        enddo
        if(iitype.eq.0) then
            write(6,*) "itype not found, stop", iatom(ia),ia
            stop
        endif
        ityatom(ia)=iitype
        numref(ia)=nref_type(iitype)
        totNel = totNel + z_type(iitype)
        nref_tot_8=nref_tot_8+nref_type(iitype)
        if(nref_type(iitype).gt.mref_8) mref_8=nref_type(iitype)
    end do

    close(9)
    !*************************************************
    nref_tot_tmp=nref_tot_8
    natom_tmp=natom

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    open(23,file="OUT.EIGEN",form="unformatted")
    rewind(23)
    read(23) islda,nkpt,mx,nref_tot_8,natom,nnodes


    !cccc Note, usually, eigen_all.store is from previous step calculation, its nref_tot_8 is smaller than in the bpsiifilexxxx

    nref_tot_8=nref_tot_tmp          ! nref_tot_8 in eigen_all, bpsiiofile100xx are supposed to be different
    !cccccc nref_tot_8 will be used below, we need to use the right nref_tot_8
    allocate(weighkpt_2(nkpt))
    allocate(akx_2(nkpt))
    allocate(aky_2(nkpt))
    allocate(akz_2(nkpt))
    allocate(E_st(mx,nkpt,islda))
    do iislda=1,islda
        do kpt=1,nkpt
            read(23) iislda_tmp,kpt_tmp,weighkpt_2(kpt),&
                &       akx_2(kpt),aky_2(kpt),akz_2(kpt)
            read(23) (E_st(i,kpt,iislda),i=1,mx)
        enddo
    enddo
    close(23)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    Emin=100000
    Emax=-100000
    do iislda=1,islda
        do kpt=1,nkpt
            do m=1,mx
                if(E_st(m,kpt,iislda).lt.Emin) Emin=E_st(m,kpt,iislda)
                if(E_st(m,kpt,iislda).gt.Emax) Emax=E_st(m,kpt,iislda)
            enddo
        enddo
    end do
    Emin=Emin-2.d0
    Emax=Emax+4.d0

    !       write(6,*) "input E_b (eV), Gaussian broadening"
    !       read(5,*) E_b
    !       E_b=0.025

    Emin=Emin-4*E_b
    Emax=Emax+4*E_b


    dE=(Emax-Emin)/nE
    nb=2*E_b/dE
    if(nb<1) then
        nb=1
    endif

    write(6,*) "Emin,Emax,nb=",Emin,Emax,nb

    !do ia=1,ntype
    !    write(6,*)"ia,spd_DOS",iiatom(ia),is_ref(ia),ip_ref(ia),id_ref(ia)
    !enddo


    sum0=0.d0
    do i=-nb,nb
        fact=exp(-(i*2.d0/nb)**2)
        sum0=sum0+fact
    enddo
    fact_norm=1.d0/(sum0*dE)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    mx_n=mx/nnodes+1
    allocate(beta_psi_tmp(nref_tot_8,mx_n))
    allocate(beta_psi_8(nref_tot_8,mx,nkpt,islda))

    write(6,*) "nref_tot_8=",nref_tot_8

    do iislda=1,islda
        do kpt=1,nkpt

            fname="bpsiiofil"
            kpt1=mod(kpt,10)
            kpt2=mod((kpt-kpt1)/10,10)
            kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
            kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

            open(10,file=&
                & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//&
                &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),&
                &       form="unformatted")
            !        read(10) nref_tot_8_tmp,mx_n_tmp,nnodes_tmp   ! to be added in later
            !       if(nref_tot_8_tmp.ne.nref_tot_8) then
            !       write(6,*)
            !       "nref_tot_8_tmp.ne.nref_tot_8,stop",nref_tot_8_tmp,nref_tot_8
            !        stop
            !       endif

            do inode=1,nnodes
                read(10) beta_psi_tmp
                do im=1,mx_n
                    m=im+(inode-1)*mx_n
                    if(m.le.mx) then
                        beta_psi_8(:,m,kpt,iislda)=beta_psi_tmp(:,im)


                    endif
                enddo
            enddo
            close(10)
        enddo
    enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if(iocc.ne.0) then
        allocate(occ(mx,nkpt,islda))
        open(23,file="IN.OCC_ADIA",form="unformatted",iostat=ierr)
        if(ierr .ne.0) then
            write(*,*) "Error reading IN.OCC_ADIA"
            stop
        else
            rewind(23)
            read(23) islda_tmp,nkpt_tmp,mx_tmp
            if(islda_tmp .ne. islda .or. nkpt_tmp .ne. nkpt .or. mx_tmp .ne. mx) then
                write(*,*) "Wrong file IN.OCC_ADIA"
                stop
            endif
            do iislda=1,islda
                do kpt=1,nkpt
                    read(23) (occ(i,kpt,iislda),i=1,mx)
                enddo
            enddo
        endif
        close(23)
    endif

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if(interp.eq.1) then

        !       nm1=4
        !       nm2=4
        !       nm3=4
        open(25,file="OUT.overlap_uk",form="unformatted",action="READ")
        rewind(25)
        read(25) hse_nq1,hse_nq2,hse_nq3,num_ext_kpt,nkpt_t,islda_t,mx_t
        if(nkpt_t.ne.nkpt.or.islda_t.ne.islda.or.mx_t.ne.mx) then
            write(6,*) "nkpt,islda,mx changed, stop", &
                    &     nkpt_t,nkpt,islda_t,islda,mx_t,mx
            stop
        endif
        allocate(S_m(mx,mx,8))
        allocate(hh(mx,mx,0:nm1,0:nm2,0:nm3))
        Lwork=mx*(2*mx-1)
        allocate(work(Lwork))
        allocate(rwork(3*mx-2))
        allocate(EE(mx,0:nm1,0:nm2,0:nm3))

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

        DOS=0.d0
        DOS_ref=0.d0
        DOS_ref_all=0.d0
        DOS_ref_ll=0.d0
        write(6,*)  "ikpt,err: ikpt,uk_ov_err,err(1:8,nearest_neigh)"
        allocate(weight_ii(8,ntype,mx,8))
        allocate(weight_ll_ii(9,8,ntype,mx,8))
        allocate(weight_ii_all(8,ntype,mx,8))
        allocate(weight_8(8,ntype,mx))
        allocate(weight_ll_8(9,8,ntype,mx))
        allocate(weight_8_all(8,ntype,mx))


        !call omp_set_num_threads(4)
        call omp_set_dynamic(.true.)
        do iislda=1,islda
            do iikpt=1,num_ext_kpt

                !stime=omp_get_wtime()
                do ii=1,8
                    read(25) k1,k2,k3,j123(1,ii),j123(2,ii),j123(3,ii),&
                    & ispin_ind(ii),ikpt_ind(ii),ikpt_ext_ind(ii)
                    !ccccc ikpt_ind(ii): the original k-pint index
                    !ccccc ikpt_ext_ind(ii): the extended k-point index
                    ! S_m(m1,m2,ii): is dconjg(u_k(m1,kpt=1))*u_k(m2,kpt=ii)
                    ind_ii(j123(1,ii),j123(2,ii),j123(3,ii))=ii
                    read(25) (S_m(:,:,ii))
                enddo
                !etime=omp_get_wtime()
                !write(*,*) "time1 = ", etime-stime


                !ccccccccccccccccccccccccccccccccccccccccccc
                !stime=omp_get_wtime()
                weight_ii=0.d0
                weight_ii_all=0.d0
                weight_ll_ii=0.d0
                do m=1,mx
                    do ii=1,8
                        kpt=ikpt_ind(ii)
                        iref=0
                        do ia=1,natom
                            itype=ityatom(ia) 
                            do ichi=1,nchi(itype)
                                do ll=1,lll(ichi,itype)*2+1   
                                    iref=iref+1
                                    if(iocc==0) then
                                        weight_ii(ichi,itype,m,ii)=weight_ii(ichi,itype,m,ii)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)
                                        weight_ll_ii(ll,ichi,itype,m,ii)=weight_ll_ii(ll,ichi,itype,m,ii)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)
                                    else
                                        weight_ii(ichi,itype,m,ii)=weight_ii(ichi,itype,m,ii)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)*occ(m,kpt,iislda)/weighkpt_2(kpt)/(2.d0/islda)
                                        weight_ll_ii(ll,ichi,itype,m,ii)=weight_ll_ii(ll,ichi,itype,m,ii)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)*occ(m,kpt,iislda)/weighkpt_2(kpt)/(2.d0/islda)
                                    endif
                                    weight_ii_all(ichi,itype,m,ii)=weight_ii_all(ichi,itype,m,ii)+abs(beta_psi_8(iref,m,kpt,iislda))**2
                                enddo
                            enddo
                        enddo    ! ia
                    enddo    ! ii
                enddo    ! m
                !etime=omp_get_wtime()
                !write(*,*) "time2 = ", etime-stime
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !stime=omp_get_wtime()
                !$OMP parallel &
                !$omp private(ii,m1,m2,cc,m3,sum0) &
                !$omp shared(mx,S_m,error_ii)
                !$omp do 
                do ii=1,8
                    error_ii(ii)=0.d0
                    do m1=1,mx
                        do m2=1,m1-1
                            cc=cmplx(0.d0,0.d0)
                            do m3=1,mx
                                cc=cc+S_m(m3,m1,ii)*conjg(S_m(m3,m2,ii))
                            enddo
                            do m3=1,mx
                                S_m(m3,m1,ii)=S_m(m3,m1,ii)-cc*S_m(m3,m2,ii)
                            enddo
                        enddo ! m2
                        sum0=0.d0
                        do m3=1,mx
                            sum0=sum0+abs(S_m(m3,m1,ii))**2
                        enddo
                        error_ii(ii)=error_ii(ii)+abs(sum0-1.d0)

                        sum0=1/dsqrt(sum0)
                        do m3=1,mx
                            S_m(m3,m1,ii)=S_m(m3,m1,ii)*sum0
                        enddo
                    enddo  ! m1 
                enddo  ! ii
                !$omp end do 
                !$omp end parallel
                !etime=omp_get_wtime()
                !write(*,*) "time3 = ", etime-stime
                error=sum(error_ii(1:8))/mx/8
                error_ii=error_ii/mx
                write(6,778) iikpt,error,(error_ii(ii),ii=1,8)
                778    format("ikpt,err ",i5,2x,f8.5,2x,8(f7.4,1x))


                       !stime=omp_get_wtime()
                hh=0.d0
                !$OMP parallel &
                !$omp private(i1,i2,i3,x1,x2,x3,ii,w1,kpt,m1,m2,m3) &
                !$omp shared(nm1,nm2,nm3,j123,ikpt_ind,mx,hh,S_m,E_st,iislda)
                !$omp do collapse(3)
                do i1=0,nm1-1
                    do i2=0,nm2-1
                        do i3=0,nm3-1
                            x1=i1*1.d0/nm1
                            x2=i2*1.d0/nm2
                            x3=i3*1.d0/nm3
                            do ii=1,8
                                w1=abs(1-x1-j123(1,ii))*abs(1-x2-j123(2,ii))*abs(1-x3-j123(3,ii))
                                kpt=ikpt_ind(ii)
                                do m1=1,mx
                                    do m2=1,mx
                                        do m3=1,mx
                                            hh(m1,m2,i1,i2,i3)=hh(m1,m2,i1,i2,i3)+w1*E_st(m3,kpt,iislda)*conjg(S_m(m1,m3,ii))*S_m(m2,m3,ii)
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
                !$omp end do 
                !$omp end parallel
                !etime=omp_get_wtime()
                !write(*,*) "time4 = ", etime-stime

                !stime=omp_get_wtime()

                !$OMP parallel &
                !$omp private(i1,i2,i3,work,rwork,info) &
                !$omp shared(nm1,nm2,nm3,mx,hh,EE,lwork)
                !$omp do collapse(3)
                do i1=0,nm1-1
                    do i2=0,nm2-1
                        do i3=0,nm3-1

                            !  S_m(m1,m3) is the m3th state at ii, in the basis of ii=0 (m1 is for
                            !  ii=0, m3 is for ii)

                            !       call zheev('N','U',mx,hh,mx,EE,work,lwork,rwork,info)   ! eigen energy only
                            call zheev('V','U',mx,hh(1,1,i1,i2,i3),mx,EE(1,i1,i2,i3),work,lwork,rwork,info)   ! eigen vectors
                        enddo
                    enddo
                enddo
                !$omp end do 
                !$omp end parallel

                !etime=omp_get_wtime()
                !write(*,*) "time5 = ", etime-stime

                !stime=omp_get_wtime()
                !$OMP parallel &
                !$omp private(i1,i2,i3,weight_8,weight_ll_8,weight_8_all,ii,w1,m1,m2,m3,cc,itype,ichi,Ep,iEp,i,fact) &
                !$omp shared(nm1,nm2,nm3,mx,hh,EE,lwork,j123,ntype,nchi,weight_ii,weight_ii_all,Emin,dE,nb) &
                !$omp shared(fact_norm,num_exp_kpt,DOS_ref,DOS_ref_all,DOS_ref_ll,upfpsp,iislda)
                !$omp do collapse(3)
                do i1=0,nm1-1
                    do i2=0,nm2-1
                        do i3=0,nm3-1

                            weight_8=0.d0
                            weight_ll_8=0.d0
                            weight_8_all=0.d0
                            do ii=1,8
                                w1=abs(1-x1-j123(1,ii))*abs(1-x2-j123(2,ii))*abs(1-x3-j123(3,ii))
                                do m1=1,mx

                                    do m2=1,mx    ! m2 is the uk basis at ii
                                        cc=cmplx(0.d0,0.d0)

                                        do m3=1,mx
                                            cc=cc+hh(m3,m1,i1,i2,i3)*conjg(S_m(m3,m2,ii))    !  S_m is a unitary transformation
                                        enddo
                                        do itype=1,ntype 
                                            do ichi=1,nchi(itype)
                                                weight_8(ichi,itype,m1)=weight_8(ichi,itype,m1)+weight_ii(ichi,itype,m2,ii)*abs(cc)**2*w1
                                                weight_8_all(ichi,itype,m1)=weight_8_all(ichi,itype,m1)+weight_ii_all(ichi,itype,m2,ii)*abs(cc)**2*w1
                                                weight_ll_8(:,ichi,itype,m1)=weight_ll_8(:,ichi,itype,m1)+weight_ll_ii(:,ichi,itype,m2,ii)*abs(cc)**2*w1
                                            enddo
                                        enddo ! itype
                                    enddo ! m2
                                enddo ! m1
                            enddo ! ii

                            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                            do m3=1,mx
                                Ep=EE(m3,i1,i2,i3)
                                iEp=(Ep-Emin)/dE+1
                                do i=-nb,nb
                                    fact=exp(-(i*2.d0/nb)**2)*fact_norm/(nm1*nm2*nm3*num_ext_kpt)
                                    ip=iEp+i
                                    if(ip<1) cycle
                                    if(ip>nE) cycle
                                    DOS(ip,iislda)=DOS(ip,iislda)+fact

                                    do itype=1,ntype
                                        do ichi=1,nchi(itype)
                                            DOS_ref(ip,ichi,itype,iislda)=DOS_ref(ip,ichi,itype,iislda)+fact*weight_8(ichi,itype,m3)
                                            DOS_ref_all(ip,ichi,itype,iislda)=DOS_ref_all(ip,ichi,itype,iislda)+fact*weight_8_all(ichi,itype,m3)
                                            DOS_ref_ll(ip,:,ichi,itype,iislda)=DOS_ref_ll(ip,:,ichi,itype,iislda)+fact*weight_ll_8(:,ichi,itype,m3)
                                        enddo

                                    enddo
                                enddo
                                !cccccccccccccccccccccccccccccccccccccccccccc
                            end do

                        end do
                    end do
                end do
                !$omp end do 
                !$omp end parallel
                !etime=omp_get_wtime()
                !write(*,*) "time6 = ", etime-stime

            enddo
        end do


    endif   ! interp.eq.1
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if(interp.eq.0) then


        DOS=0.d0
        DOS_ref=0.d0
        DOS_ref_all=0.d0
        DOS_ref_ll=0.d0

        do iislda=1,islda
            do kpt=1,nkpt
                do m=1,mx
                    weight=0.d0
                    weight_ll=0.d0
                    weight_all=0.d0
                    !cccccccccccccccccccccccccccccccc
                    iref=0
                    do ia=1,natom
                        itype=ityatom(ia) 
                        do ichi=1,nchi(itype)
                            do ll=1,lll(ichi,itype)*2+1   
                                iref=iref+1
                                if(iocc==0) then
                                    weight(ichi,itype)=weight(ichi,itype)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)
                                    weight_ll(ll,ichi,itype)=weight_ll(ll,ichi,itype)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)
                                else
                                    weight(ichi,itype)=weight(ichi,itype)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)*occ(m,kpt,iislda)/weighkpt_2(kpt)/(2.d0/islda)
                                    weight_ll(ll,ichi,itype)=weight_ll(ll,ichi,itype)+abs(beta_psi_8(iref,m,kpt,iislda))**2*w(ia)*occ(m,kpt,iislda)/weighkpt_2(kpt)/(2.d0/islda)
                                endif
                                weight_all(ichi,itype)=weight_all(ichi,itype)+abs(beta_psi_8(iref,m,kpt,iislda))**2
                            enddo
                        enddo
                    enddo    ! ia

                    !ccccccccccccccccccccccccccccccccccccccc
                    Ep=E_st(m,kpt,iislda)
                    iEp=(Ep-Emin)/dE+1

                    do i=-nb,nb
                        fact=exp(-(i*2.d0/nb)**2)*fact_norm
                        ip=iEp+i
                        DOS(ip,iislda)=DOS(ip,iislda)+fact*weighkpt_2(kpt)

                        do itype=1,ntype
                            do ichi=1,nchi(itype)
                            DOS_ref(ip,ichi,itype,iislda)=DOS_ref(ip,ichi,itype,iislda)+fact*weighkpt_2(kpt)*weight(ichi,itype)
                            DOS_ref_all(ip,ichi,itype,iislda)=DOS_ref_all(ip,ichi,itype,iislda)+fact*weighkpt_2(kpt)*weight_all(ichi,itype)
                            DOS_ref_ll(ip,:,ichi,itype,iislda)=DOS_ref_ll(ip,:,ichi,itype,iislda)+fact*weighkpt_2(kpt)*weight_ll(:,ichi,itype)
                            enddo
                        enddo

                    enddo   ! i
                    !cccccccccccccccccccccccccccccccccccccccccccc

                end do
            end do
        end do

    endif   ! interp.eq.0
    !cccccccccccccccccccccccccccccccccccccccccccc
    !ccccc s,p,d, need to rescale by the (psi*dV(r))^2 factor
    !ccccccccccccc   rescale
!    if(ipart_DOS.eq.0) then
        do iislda=1,islda
            do i=1,nE
                sum0=0.d0
                do itype=1,ntype
                    do ichi=1,nchi(itype)
                       sum0=sum0+DOS_ref_all(i,ichi,itype,iislda)
                    enddo
                enddo
                if(sum0.gt.1.D-20) then
                    fact=DOS(i,iislda)/sum0
                else
                    fact=1.d0
                endif

               do itype=1,ntype
                   do ichi=1,nchi(itype)
                       DOS_ref(i,ichi,itype,iislda)=DOS_ref(i,ichi,itype,iislda)*fact
                       DOS_ref_ll(i,:,ichi,itype,iislda)=DOS_ref_ll(i,:,ichi,itype,iislda)*fact
                   enddo
               enddo
            enddo
        enddo

        !!!! DOS recalculated here
        
        DOS = 0.0d0
        do iislda=1,islda
            do i=1,nE
                do itype=1,ntype
                    do ichi=1,nchi(itype)
                    DOS(i,iislda) = DOS(i,iislda) + DOS_ref(i,ichi,itype,iislda)
                    enddo
                enddo
            enddo
        enddo

!    endif
    !cccccccccccccccccccccccccccccccccccccc

    if(islda.eq.1) then
        open(10,file="DOS.totalspin")
        rewind(10)

        write(10,'(A15,$)') " # Energy"
        write(10,'(1x,A15,$)') "Total"
        do itype = 1, ntype
            do ichi=1,nchi(itype)
                atom_name= upfpsp(itype)%psd
                atom_name = adjustr(atom_name)
                if(upfpsp(itype)%els(ichi).ne."") then
                    label="-"//adjustl(trim(upfpsp(itype)%els(ichi)))
                else
                    if(lll(ichi,itype) .eq. 0) label="-s"
                    if(lll(ichi,itype) .eq. 1) label="-p"
                    if(lll(ichi,itype) .eq. 2) label="-d"
                    if(lll(ichi,itype) .eq. 3) label="-f"
                endif
                write(10,'(1x,A15,$)') trim(trim(atom_name)//trim(label))
            enddo
        end do
        write(10,*) 

         do i = 1, nE
             write (10,200) Emin+(i-1)*dE, DOS(i,1),((DOS_ref(i,ichi,itype,1),ichi=1,nchi(itype)),itype=1,ntype)
         end do
        close(10)
    elseif(islda.eq.2) then
        open(10,file="DOS.totalspin")
        rewind(10)

        write(10,'(A15,$)') " # Energy"
        write(10,'(1x,A15,$)') "Total"
        do itype = 1, ntype
        do ichi = 1, nchi(itype)
            atom_name= upfpsp(itype)%psd
            atom_name = adjustr(atom_name)
            if(upfpsp(itype)%els(ichi).ne."") then
                label="-"//adjustl(upfpsp(itype)%els(ichi))
            else
                if(lll(ichi,itype) .eq. 0) label="-s"
                if(lll(ichi,itype) .eq. 1) label="-p"
                if(lll(ichi,itype) .eq. 2) label="-d"
                if(lll(ichi,itype) .eq. 3) label="-f"
            endif
            write(10,'(1x,A15,$)') trim(trim(atom_name)//trim(label))
        enddo
        end do
        write(10,*) 

         do i = 1, nE
             write (10,200) Emin+(i-1)*dE, DOS(i,1)+DOS(i,2),((DOS_ref(i,ichi,itype,1)+DOS_ref(i,ichi,itype,2),ichi=1,nchi(itype)),itype=1,ntype)
         end do
        close(10)
        open(10,file="DOS.spinup")
        rewind(10)

        write(10,'(A15,$)') " # Energy"
        write(10,'(1x,A15,$)') "Total"
        do itype = 1, ntype
            do ichi = 1, nchi(itype)
                atom_name= upfpsp(itype)%psd
                atom_name = adjustr(atom_name)
                if(upfpsp(itype)%els(ichi).ne."") then
                    label="-"//adjustl(upfpsp(itype)%els(ichi))
                else
                    if(lll(ichi,itype) .eq. 0) label="-s"
                    if(lll(ichi,itype) .eq. 1) label="-p"
                    if(lll(ichi,itype) .eq. 2) label="-d"
                    if(lll(ichi,itype) .eq. 3) label="-f"
                endif
                write(10,'(1x,A15,$)') trim(trim(atom_name)//trim(label))
            enddo
        end do
        write(10,*) 

         do i = 1, nE
             write (10,200) Emin+(i-1)*dE, DOS(i,1),((DOS_ref(i,ichi,itype,1),ichi=1,nchi(itype)),itype=1,ntype)
         end do
        close(10)
        open(10,file="DOS.spindown")
        rewind(10)

        write(10,'(A15,$)') " # Energy"
        write(10,'(1x,A15,$)') "Total"
        do itype = 1, ntype
            do ichi = 1, nchi(itype)
                atom_name= upfpsp(itype)%psd
                atom_name = adjustr(atom_name)
                if(upfpsp(itype)%els(ichi).ne."") then
                    label="-"//adjustl(upfpsp(itype)%els(ichi))
                else
                    if(lll(ichi,itype) .eq. 0) label="-s"
                    if(lll(ichi,itype) .eq. 1) label="-p"
                    if(lll(ichi,itype) .eq. 2) label="-d"
                    if(lll(ichi,itype) .eq. 3) label="-f"
                endif
                write(10,'(1x,A15,$)') trim(trim(atom_name)//trim(label))
            enddo
        end do
        write(10,*)

         do i = 1, nE
             write (10,200) Emin+(i-1)*dE, -DOS(i,2),((-DOS_ref(i,ichi,itype,2),ichi=1,nchi(itype)),itype=1,ntype)
         end do
        close(10)
    endif

200    format(E15.5,1x,E15.5,1x,*(E15.5,1x))
    !200    format(E12.5,2x,E11.4,2x,15(E11.4,1x))
    !200    format(E12.5,2x,F11.7,2x,15(F11.7,1x))

   !plot projected DOS
   if(islda==2) then
       k=1
       filename='DOS.totalspin_projected'
       write(6,*) "DOS data in file: ", ADJUSTL(trim(filename))
       open (unit = 44+k, file = trim(filename))
       write(44+k,'(A15,$)') " # Energy"
       write(44+k,'(1x,A15,$)') "Total"

       do itype = 1, ntype
           do ichi=1,nchi(itype)
               atom_name= upfpsp(itype)%psd
               if(upfpsp(itype)%els(ichi).ne."") then
                   label="-"//adjustl(upfpsp(itype)%els(ichi))
               else
                   if(lll(ichi,itype) .eq. 0) label="-s"
                   if(lll(ichi,itype) .eq. 1) label="-p"
                   if(lll(ichi,itype) .eq. 2) label="-d"
                   if(lll(ichi,itype) .eq. 3) label="-f"
               endif
               if(lll(ichi,itype) .eq. 0)  write(44+k,'(1x,A15,$)') trim(trim(atom_name)//trim(label))
               if(lll(ichi,itype) .eq. 1) then
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"x")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"y")
               endif
               if(lll(ichi,itype) .eq. 2) then
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z2")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xz")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"yz")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"(x^2-y^2)")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xy")
               endif
               if(lll(ichi,itype) .eq. 3) then
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z^3")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xz^2")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"yz^2")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z(x^2-y^2)")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xyz")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"x(x^2-3y^2)")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"y(3x^2-y^2)")
               endif
           enddo
       enddo
       write(44+k,*) 

       do i = 1, nE
           write(44+k,'(E15.5,1x,E15.5,$)')  Emin+(i-1)*dE,DOS(i,1)+DOS(i,2)
           do itype=1,ntype
               do ichi=1,nchi(itype)
                   ll=lll(ichi,itype)
                   write(44+k,'(1x,E15.5,$)') DOS_ref_ll(i,1:(2*ll+1),ichi,itype,1)+DOS_ref_ll(i,1:(2*ll+1),ichi,itype,2)
               enddo
           enddo
           write(44+k,*)
       end do
       close(10)
   endif
   do k = 1, islda
       if(islda==1) then
           filename='DOS.totalspin_projected'
       else
           if(k==1) filename='DOS.spinup_projected'
           if(k==2) filename='DOS.spindown_projected'
       endif
       write(6,*) "DOS data in file: ", ADJUSTL(trim(filename))
       open (unit = 44+k, file = trim(filename))
       write(44+k,'(A15,$)') " # Energy"
       write(44+k,'(1x,A15,$)') "Total"

       do itype = 1, ntype
           do ichi=1,nchi(itype)
               atom_name= upfpsp(itype)%psd
               if(upfpsp(itype)%els(ichi).ne."") then
                   label="-"//adjustl(upfpsp(itype)%els(ichi))
               else
                   if(lll(ichi,itype) .eq. 0) label="-s"
                   if(lll(ichi,itype) .eq. 1) label="-p"
                   if(lll(ichi,itype) .eq. 2) label="-d"
                   if(lll(ichi,itype) .eq. 3) label="-f"
               endif
               if(lll(ichi,itype) .eq. 0)  write(44+k,'(1x,A15,$)') trim(trim(atom_name)//trim(label))
               if(lll(ichi,itype) .eq. 1) then
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"x")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"y")
               endif
               if(lll(ichi,itype) .eq. 2) then
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z2")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xz")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"yz")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"(x^2-y^2)")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xy")
               endif
               if(lll(ichi,itype) .eq. 3) then
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z^3")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xz^2")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"yz^2")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"z(x^2-y^2)")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"xyz")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"x(x^2-3y^2)")
                   write(44+k,'(1x,A15,$)') trim(trim(trim(atom_name)//trim(label))//"y(3x^2-y^2)")
               endif
           enddo
       enddo
       write(44+k,*) 

       do i = 1, nE
           if(k.eq.1) then
               write(44+k,'(E15.5,1x,E15.5,$)')  Emin+(i-1)*dE,DOS(i,k)
           else
               write(44+k,'(E15.5,1x,E15.5,$)')  Emin+(i-1)*dE,-DOS(i,k)
           endif
           do itype=1,ntype
               do ichi=1,nchi(itype)
                   ll=lll(ichi,itype)
                   if(k.eq.1) then
                       write(44+k,'(1x,E15.5,$)') DOS_ref_ll(i,1:(2*ll+1),ichi,itype,k)
                   else
                       write(44+k,'(1x,E15.5,$)') -DOS_ref_ll(i,1:(2*ll+1),ichi,itype,k)
                   endif
               enddo
           enddo
           write(44+k,*)
       end do
       close(10)

   enddo


    !*************************************************
contains

    !*******************************************
    subroutine readvwr_head()

        implicit double precision (a-h,o-z)


        open(10,file=vwr_atom(ia),status='old',action='read',iostat=ierr)
        if(ierr.ne.0) then
            write(6,*) "vwr_file ***",filename,"*** does not exist, stop"
            stop
        endif
        read(10,*) nrr_t,ic_t,iiatom(ia),zatom_t,iloc_t,occ_s_t,&
            &  occ_p_t,occ_d_t
        read(10,*) is_ref_t,ip_ref_t,id_ref_t,&
            &  is_TB_t,ip_TB_t,id_TB_t
        close(10)


        if(iloc_t.eq.1) is_ref_t=0
        if(iloc_t.eq.2) ip_ref_t=0
        if(iloc_t.eq.3) id_ref_t=0
        !ccccccccccccccccccccccccccccccccccccccccccccccc
        is_ref_t=1
        ip_ref_t=1
        id_ref_t=1
        !cccccccccccccccccccccccccccccccccccccccccccccccccc

        nref_type(ia)=is_ref_t+ip_ref_t*3+id_ref_t*5
        is_ref(ia)=is_ref_t
        ip_ref(ia)=ip_ref_t
        id_ref(ia)=id_ref_t

        return
    end subroutine readvwr_head


    !ccccccccccc

    subroutine readf_xatom_new(w)
        implicit double precision (a-h,o-z)
        real*8, allocatable, dimension(:,:) :: xatom_tmp
        integer, allocatable, dimension(:,:) :: imov_at_tmp
        integer, allocatable, dimension(:)   ::  iatom_tmp
        real*8, allocatable, dimension(:) :: w_tmp
        real*8 :: w(:)


        if(ipart_DOS.eq.1) then
            write(6,*)  "For atom based partial DOS, there must be an"//&
                &  " extra last column (atomic weight) in ", f_xatom 
        endif

        open(10,file=f_xatom,status='old',action='read',iostat=ierr)

        if(ierr.ne.0) then
            if(inode.eq.1) write(6,*) "IN.ATOM",f_xatom,"not exist, stop"
            stop
        endif

        rewind(10)
        read(10,*) natom
        if(natom.gt.matom_1) then
            if(inode.eq.1) then
                write(6,*) "natom.gt.matom_1, increase matom_1 in data.f, stop"
                write(6,*) f_xatom,natom,matom_1
            endif
            stop
        endif

        allocate(xatom_tmp(3,natom))
        allocate(imov_at_tmp(3,natom))
        allocate(iatom_tmp(natom))
        allocate(w_tmp(natom))

        !      CALL read_key_words ( 10, 'Lattice',LEN('Lattice'), right,IERR)
        CALL scan_key_words ( 10, 'Lattice',LEN('Lattice'), right,readit)
        !      if(ierr < 0) then
        !         rewind(10)
        !         read(10, *)
        !      endif

        read(10,*) (AL(i,1),i=1,3)
        read(10,*) (AL(i,2),i=1,3)
        read(10,*) (AL(i,3),i=1,3)

        !AL = AL/A_AU_1 !zhilin

        !      CALL read_key_words(10,'Position',LEN('Position'),right, IERR)
        CALL scan_key_words(10,'Position',LEN('Position'),right, readit)
        !      if(ierr < 0) then
        !         write(*,*) "keyword 'position' is needed at",f_xatom
        !         stop
        !      endif
        if(ipart_DOS.eq.0) then
            do i=1,natom
                read(10,*) iatom_tmp(i),&
                    &     xatom_tmp(1,i),xatom_tmp(2,i),xatom_tmp(3,i),&
                    &     imov_at_tmp(1,i),imov_at_tmp(2,i),imov_at_tmp(3,i)
                w_tmp(i)=1.d0
            enddo
        else
            do i=1,natom
                read(10,*) iatom_tmp(i),&
                    &     xatom_tmp(1,i),xatom_tmp(2,i),xatom_tmp(3,i),&
                    &     imov_at_tmp(1,i),imov_at_tmp(2,i),imov_at_tmp(3,i),&
                    &     w_tmp(i)
            enddo
        endif
        close(10)


        !cccccccccccccccccccccccccccccccccccccccccccccccccccc
        !cccc Now, re-arrange xatom, so the same atoms are consequentive together. 
        !cccc This is useful to speed up the getwmask.f
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ii=0
        100    continue
        ncount=0
        itype_tmp=-2
        do i=1,natom
            if(itype_tmp.eq.-2.and.iatom_tmp(i).ne.-1)&
                &  itype_tmp=iatom_tmp(i)

            if(iatom_tmp(i).eq.itype_tmp) then
                ii=ii+1
                ncount=ncount+1
                iatom(ii)=iatom_tmp(i)
                iatom_tmp(i)=-1
                xatom(:,ii)=xatom_tmp(:,i)
                imov_at(:,ii)=imov_at_tmp(:,i)
                w(ii)=w_tmp(i)
            endif
        enddo
        if(ncount.gt.0) goto 100
        if(ii.ne.natom) then
            write(6,*) "something wrong to rearrange xatom, stop"
            stop
        endif


        deallocate(xatom_tmp)
        deallocate(imov_at_tmp)
        deallocate(iatom_tmp)
        deallocate(w_tmp)

        return
    end subroutine readf_xatom_new
    !**************************************************

    subroutine error_stop(message)
        implicit double precision(a-h,o-z)
        character*200 message
        write(6,*) "error in etot.input,stop"
        write(6,*) message
        stop

        return
    end subroutine error_stop



    end



