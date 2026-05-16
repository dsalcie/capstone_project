!
!   Copyright (c) 2013-2022 Stefan Typel, Marco Mancini, Micaela Oertel
!
!   This file is part of CompOSE
!
!   CompOSE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   CompOSE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with CompOSE; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!


module m_get_tables
 implicit none

 private

 public :: init_eos_table_term, read_eos_4_tables


 private :: read_eos_tables_tnyb, read_eos_table_thermo,&
   &        read_eos_table_compo, read_eos_table_micro


contains
 !***********************************************************************
 subroutine read_eos_4_tables(iwr,nbl,iyq,ii_thermo,ii_tynb,unit,iinit)
  ! MMancini for the CompOSE core team, version 2.03, 2016/11/08
  use omp_lib
  integer,intent(in) :: iwr,nbl,ii_thermo,ii_tynb,unit,iinit
  integer,intent(out) :: iyq
  integer :: nthreads
#ifdef DEBUG
  double precision :: timei,timef
#endif


  call read_eos_tables_tnyb(iwr,ii_tynb,unit,iinit,iyq)

#ifdef DEBUG
  timei = omp_get_wtime()
#endif

  if(unit/=0)then
    nthreads = 1
  else
    nthreads = 3
  endif

  ! open omp parallel region with only nthreads threads

  !$OMP PARALLEL num_threads(nthreads)
  !$OMP SINGLE
  !$OMP TASK
  call read_eos_table_thermo(iwr,nbl,ii_thermo,unit,iyq)
  !$OMP END TASK
  !$OMP TASK
  call read_eos_table_compo(iwr,nbl,ii_thermo,unit)
  !$OMP END TASK
  !$OMP TASK
  call read_eos_table_micro(iwr,nbl,ii_thermo,unit)
  !$OMP END TASK
  !$OMP END SINGLE
  !$OMP END PARALLEL

#ifdef DEBUG
  timef = omp_get_wtime()
  print '(2x,a,f12.5,a,i3,a)','TIME for read tables : ',&
    &    timef-timei,  '(s) with ', nthreads, ' threads'
#endif


 end subroutine read_eos_4_tables
 !***********************************************************************
 subroutine init_eos_table_term(iwr,iterm,iinit)
  ! Stefan Typel for the CompOSE core team, version 1.04, 2017/11/16
  !use compose_internal, only : error_msg
  integer,intent(in) :: iwr,iterm
  integer,intent(out) ::  iinit
  integer :: nbl,iunit2,ierror,inew,iyq

  ! maximum number of subtables
  nbl = 10

  iunit2 = 31

  open(unit=iunit2,file='eos.init',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    !   write(*,*)
    !   write(*,*) ' The file eos.init exists already.'
    !   write(*,*) ' Do you want to generate a new file eos.init?'
    !   write(*,*) ' Please select:'
    !   write(*,*) ' 1: Yes'
    !   write(*,*) ' else: No'
    !   read(5,*) inew
    inew = 0
  else
    inew = 1
  end if
  close(unit=iunit2)
  if (iterm == 999) inew = 1

  if (inew == 1) then
    open(unit=iunit2,file='eos.init',&
      status='unknown',action='write',iostat=ierror)

    if (iterm /= 999) then
      write(*,*)
      write(*,*) ' The file eos.init does not exist.'
    end if
    write(*,*)
    write(*,*) ' generating a new file eos.init'

    call read_eos_4_tables(iwr,nbl,iyq,ii_thermo=1,ii_tynb=1,unit=iunit2,iinit=0)

    write(*,*)
    write(*,*) ' new file eos.init generated'
    write(*,*)

    close(unit=iunit2)
  end if

  iinit = inew

 end subroutine init_eos_table_term


 !***********************************************************************
 subroutine read_eos_tables_tnyb(iwr,ii,unit_in,iinit,iyq)
  use compose_internal, only : error_msg, para_min, para_max, dim_p, &
    &                          dim_a, tab_para, inbyq, imap, irpl, &
    &                          min_idx, dim_err, idx_arg, dim_idx, &
    &                          jmap, val_rpl

  ! Stefan Typel for the CompOSE core team, version 2.05, 2017/11/16
  integer,intent(in) :: iwr,ii,unit_in,iinit
  integer,intent(out) :: iyq
  integer :: it,in,iy,ib,ip,ipp,ippp,iunit,ierror,ierr,ipar
  integer :: alloc_status,itmp(4),itmpp(4),icnt(4),ivar(4)
  double precision :: tmp
  integer, parameter :: unit0 = 500

  ! counter for total number of read table entries
  icnt(1:4) = 0
  ivar(1:4) = 0

  ! error counter
  ierr = 0

  ! maximum dimension of parameter vector
  dim_a = 0

  ! temperature
  ip = 1
  iunit = unit0 + ip
  open(unit=iunit,file='eos.t',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    if (iwr == 1) then
      write(*,*)
      write(*,*) ' reading minimum and maximum index ',&
        'from parameter table for temperature T'
    end if
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    if (itmpp(ip) < itmp(ip)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 3
    end if
    if (itmpp(ip) > itmp(ip)) ivar(ip) = 1
  else
    ierr = ierr+1
    if (ierr < dim_err) error_msg(ierr) = 1
  end if
  close(unit=iunit)

  ! baryon number density
  ip = 2
  iunit = unit0 + ip
  open(unit=iunit,file='eos.nb',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    if (iwr == 1) then
      write(*,*)
      write(*,*) ' reading minimum and maximum index ',&
        'from parameter table for baryon number density n_b'
    end if
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    if (itmpp(ip) < itmp(ip)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 8
    end if
    if (itmpp(ip) > itmp(ip)) ivar(ip) = 1
  else
    ierr = ierr+1
    error_msg(ierr) = 6
  end if
  close(unit=iunit)

  ! hadronic charge fraction
  iyq = 0
  ip = 3
  iunit = unit0 + ip
  open(unit=iunit,file='eos.yq',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    if (iwr == 1) then
      write(*,*)
      write(*,*) ' reading minimum and maximum index ',&
        'from parameter table for hadronic charge fraction Y_q'
    end if
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    if (itmpp(ip) < itmp(ip)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 13
    end if
    if (itmpp(ip) > itmp(ip)) then
      ivar(ip) = 1
      iyq = 1
    end if
  else
    ierr = ierr+1
    if (ierr < dim_err) error_msg(ierr) = 11
  end if
  close(unit=iunit)

  ! magnetic field strength
  ip = 4
  iunit = unit0 + ip
  open(unit=iunit,file='eos.b',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    if (iwr == 1) then
      write(*,*)
      write(*,*) ' reading minimum and maximum index ',&
        'from parameter table for magnetic field strength B'
      write(*,*)
    end if
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    if (itmpp(ip) < itmp(ip)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 18
    end if
    if (itmpp(ip) > itmp(ip)) ivar(ip) = 1
  else
    itmp(ip) = 0
    itmpp(ip) = 0
    write(*,*)
    write(*,*) ' no file eos.b'
  end if
  close(unit=iunit)

  if (ivar(4) == 1) then
    if ((ivar(1)+ivar(2)+ivar(3)) == 3) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 15
    else
      if (ivar(1) == 0) then
        irpl = 1
        write(*,*) ' no dependence on temperature T'
      endif
      if (ivar(2) == 0) then
        irpl = 2
        write(*,*) ' no dependence on baryon density n_b'
      endif
      if (ivar(3) == 0) then
        irpl = 3
        write(*,*) ' no dependence on hadronic charge fraction Y_q'
      endif
      write(*,*) ' dependence on magnetic field strength B'
    end if
  else
    irpl = 0
  end if

  if ((irpl < 0).or.(irpl > 3)) then
    ierr = ierr+1
    if (ierr < dim_err) error_msg(ierr) = 5
  end if

  ipar = 0
  do ip=1,4,1
    ! dimension of parameter vectors
    dim_p(ip) = itmpp(ip)-itmp(ip)+1
    if (dim_p(ip) > dim_a) then
      dim_a = dim_p(ip)
    end if
    if (dim_p(ip) == 1) ipar = ipar+1
  end do
  if (ipar == 0) then
    ! no dependence on any parameter
    ierr = ierr+1
    if (ierr < dim_err) error_msg(ierr) = 19
  end if

  call write_errors(ierr)

  ! standard index ordering: T, n_b, Y_q
  do ip=1,3,1
    imap(ip) = ip
    jmap(ip) = ip
  end do
  jmap(4) = 0
  inbyq = 1

  if (irpl == 1) then
    ! index ordering: B, n_b, Y_q
    imap(1) = 4
    imap(2) = 2
    imap(3) = 3
    jmap(1) = 0
    jmap(2) = 2
    jmap(3) = 3
    jmap(4) = 1
    inbyq = 1
  endif
  if (irpl == 2) then
    ! index ordering: T, B, Y_q
    imap(1) = 1
    imap(2) = 4
    imap(3) = 3
    jmap(1) = 1
    jmap(2) = 0
    jmap(3) = 3
    jmap(4) = 2
    inbyq = 0
  endif
  if (irpl == 3) then
    ! index ordering: T, n_b, B
    imap(1) = 1
    imap(2) = 2
    imap(3) = 4
    jmap(1) = 1
    jmap(2) = 2
    jmap(3) = 0
    jmap(4) = 3
    inbyq = 0
  endif
  ! mapping of indices
  do ip=1,3,1
    dim_idx(ip) = dim_p(imap(ip))
  end do
  if (inbyq == 1) then
    if ((dim_idx(2) < 2).or.(dim_idx(3) < 2)) inbyq = 0
  end if

  ! allocation of arrays
  if (iinit == 0) then
    allocate(tab_para(dim_a,1:8),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 28
    end if

    allocate(idx_arg(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),0:4),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 113
    end if
  end if

  call write_errors(ierr)
  ierr = 0

  ! initialization
  if (dim_a > 0) tab_para(1:dim_a,1:8) = 0.d00
  if ((dim_idx(1) > 0).and.(dim_idx(2) > 0).and.(dim_idx(3) > 0)) then
    idx_arg(:,:,:,0:3) = 0
    idx_arg(:,:,:,4) = -1
  end if

  ! temperature
  ip = 1
  iunit = unit0 + ip
  open(unit=iunit,file='eos.t',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    do it=1,dim_p(ip),1
      read(iunit,*) tmp
      tab_para(it,ip) = tmp
      icnt(ip) = icnt(ip)+1
    end do
    do it=2,dim_p(ip),1
      if ((ierr == 0).and.(tab_para(it,ip) < tab_para(it-1,ip))) then
        ierr = ierr+1
        if (ierr < dim_err) error_msg(ierr) = 2
      end if
    end do
    if (ierr == 0) then
      min_idx(ip) = itmp(ip)
    end if
  end if
  close(unit=iunit)
  if (iwr == 1) then
    write(*,*)
    write(*,*) icnt(ip),&
      'entries of parameter table for temperatures read'
  end if
  ipp = ip+4
  do it=1,(dim_p(ip)-1),1
    tab_para(it,ipp) = tab_para(it+1,ip)-tab_para(it,ip)
  end do

  ! baryon number density
  ip = 2
  iunit = unit0 + ip
  open(unit=iunit,file='eos.nb',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    do in=1,dim_p(ip),1
      read(iunit,*) tmp
      tab_para(in,ip) = tmp
      icnt(ip) = icnt(ip)+1
    end do
    do in=2,dim_p(ip),1
      if ((ierr == 0).and.(tab_para(in,ip) < tab_para(in-1,ip))) then
         write(*,*) 'in,ip',in,ip,tab_para(in,ip),tab_para(in-1,ip)
         ierr = ierr+1
        if (ierr < dim_err) error_msg(ierr) = 7
      end if
    end do
    if (ierr == 0) then
      min_idx(ip) = itmp(ip)
    end if
  end if
  close(unit=iunit)
  if (iwr == 1) then
    write(*,*) icnt(ip),&
      'entries of parameter table for baryon number densities read'
  end if
  ipp = ip+4
  do in=1,(dim_p(ip)-1),1
    tab_para(in,ipp) = tab_para(in+1,ip)-tab_para(in,ip)
  end do

  ! hadronic charge fraction
  ip = 3
  iunit = unit0 + ip
  open(unit=iunit,file='eos.yq',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    do iy=1,dim_p(ip),1
      read(iunit,*) tmp
      tab_para(iy,ip) = tmp
      icnt(ip) = icnt(ip)+1
    end do
    do iy=2,dim_p(ip),1
      if ((ierr == 0).and.(tab_para(iy,ip) < tab_para(iy-1,ip))) then
        ierr = ierr+1
        if (ierr < dim_err) error_msg(ierr) = 12
      end if
    end do
    if (ierr == 0) then
      min_idx(ip) = itmp(ip)
    end if
  end if
  close(unit=iunit)
  if (iwr == 1) then
    write(*,*) icnt(ip),&
      'entries of parameter table for charge fraction read'
  end if
  ipp = ip+4
  do iy=1,(dim_p(ip)-1),1
    tab_para(iy,ipp) = tab_para(iy+1,ip)-tab_para(iy,ip)
  end do

  ! magnetic field strength
  ip = 4
  iunit = unit0 + ip
  open(unit=iunit,file='eos.b',&
    status='old',action='read',iostat=ierror)
  if (ierror == 0) then
    read(iunit,*) itmp(ip)
    read(iunit,*) itmpp(ip)
    do ib=1,dim_p(ip),1
      read(iunit,*) tmp
      tab_para(ib,ip) = tmp
      icnt(ip) = icnt(ip)+1
    end do
    do ib=2,dim_p(ip),1
      if ((ierr == 0).and.(tab_para(ib,ip) < tab_para(ib-1,ip))) then
        ierr = ierr+1
        if (ierr < dim_err) error_msg(ierr) = 17
      end if
    end do
    if (ierr == 0) then
      min_idx(ip) = itmp(ip)
      !      idx_min(ip) = itmp(ip)
      !      idx_max(ip) = itmpp(ip)
    end if
    if (iwr == 1) then
      write(*,*) icnt(ip),&
        'entries of parameter table for magnetic field strength read'
    end if
    ipp = ip+4
    do ib=1,(dim_p(ip)-1),1
      tab_para(ib,ipp) = tab_para(ib+1,ip)-tab_para(ib,ip)
    end do
  else
    icnt(ip) = 1
    min_idx(ip) = 1
  end if
  close(unit=iunit)

  ! consistency of dimensions
  do ip=1,4,1
    if (icnt(ip) /= dim_p(ip)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 119+ip
    end if
  end do

  call write_errors(ierr)

  if (iwr == 1) then
    write(*,*)
    write(*,*) ' maximum dimension of parameter files =',dim_a
  end if

  ! minimum and maximum values
  do ip=1,3,1
    para_min(ip) = tab_para(1,ip)
    para_max(ip) = tab_para(dim_p(ip),ip)
  end do
  ip = 4
  if (irpl > 0) then
    para_min(ip) = tab_para(1,ip)
    para_max(ip) = tab_para(dim_p(ip),ip)
  else
    para_min(ip) = 0.d00
    para_max(ip) = 0.d00
  end if

  if (iwr == 1) then
    if (ivar(1) == 1) then
      write(*,*)
      write(*,*) ' minimum temperature:',para_min(1),'MeV'
      write(*,*) ' maximum temperature:',para_max(1),'MeV'
    end if
    if (ivar(2) == 1) then
      write(*,*)
      write(*,*) ' minimum baryon number density:',para_min(2),'fm^-3'
      write(*,*) ' maximum baryon number density:',para_max(2),'fm^-3'
    end if
    if (ivar(3) == 1) then
      write(*,*)
      write(*,*) ' minimum hadronic charge fraction:',para_min(3)
      write(*,*) ' maximum hadronic charge fraction:',para_max(3)
    end if
    if (ivar(4) == 1) then
      write(*,*)
      write(*,*) ' minimum magnetic field strength:',para_min(4),'G'
      write(*,*) ' maximum magnetic field strength:',para_max(4),'G'
    end if
  end if

  ! mapping of parameter values, save value of constant parameter
  if (irpl == 0) then
    val_rpl = tab_para(1,4)
  else
    val_rpl = tab_para(1,irpl)
    do ip=1,3,1
      ipp = ip+4
      ippp = imap(ip)+4
      do ib=1,dim_p(4),1
        tab_para(ib,ip) = tab_para(ib,imap(ip))
        tab_para(ib,ipp) = tab_para(ib,ippp)
      end do
    end do
  end if

  if (ii == 1) then
    write(unit_in,*) (ivar(in),in=1,4,1),irpl
    do ip=1,4,1
      write(unit_in,*) para_min(ip),para_max(ip)
    end do
  end if

  return
 end subroutine read_eos_tables_tnyb

 !***********************************************************************
 subroutine read_eos_table_thermo(iwr,nbl,ii,unit_in,iyq)
  ! Stefan Typel for the CompOSE core team, version 2.04, 2017/11/16
  use compose_internal, only : error_msg,idx_ex, tab_thermo,idx_thermo&
    &,v_thermo,eos_thermo_add, idx_add, idx_arg, dim_err, dim_reg,&
    & incl_l,  nadd_max, m_p, m_n, nall_max, dim_idx, min_idx, imap
  implicit none
  integer,intent(in) :: iwr,nbl,ii,unit_in,iyq
  integer :: ibl,ierror,icnt,it,in,iy,it2,in2,iy2,alloc_status,&
    iv,iunit,iflag,itmp,nqty,iw,ierr,ibeta
  double precision :: tmp(dim_reg)
  double precision, allocatable :: qty(:)
  integer, parameter :: unit0 = 300

  ! indices and stored quantities in tab_thermo
  ! 1: p/n_b
  ! 2: s/n_b
  ! 3: mu_b/m_n-1
  ! 4: mu_q/m_n
  ! 5: mu_e/m_n
  ! 6: f/(n_b*m_n)-1
  ! 7: e/(n_b_m_n)-1

  ! counter for total number of read table entries
  icnt = 0


  ! error counter
  ierr = 0

  idx_ex(1) = 0
  if (ii /= 2) then
    open(unit=unit0,file='eos.thermo', status='old',action='read',iostat=ierror)
    if (ierror == 0) then
      itmp = 0
    else
      itmp = nbl-1
    end if

    do ibl=0,itmp,1
      ! flag for stop reading file
      iflag = 0
      iunit = unit0+ibl

      choose_file_0: select case (ibl)
      case (0)
        if (itmp > 0) then
          open(unit=iunit,file='eos.thermo0',&
            status='old',action='read',iostat=ierror)
        end if
      case (1)
        open(unit=iunit,file='eos.thermo1',&
          status='old',action='read',iostat=ierror)
      case (2)
        open(unit=iunit,file='eos.thermo2',&
          status='old',action='read',iostat=ierror)
      case (3)
        open(unit=iunit,file='eos.thermo3',&
          status='old',action='read',iostat=ierror)
      case (4)
        open(unit=iunit,file='eos.thermo4',&
          status='old',action='read',iostat=ierror)
      case (5)
        open(unit=iunit,file='eos.thermo5',&
          status='old',action='read',iostat=ierror)
      case (6)
        open(unit=iunit,file='eos.thermo6',&
          status='old',action='read',iostat=ierror)
      case (7)
        open(unit=iunit,file='eos.thermo7',&
          status='old',action='read',iostat=ierror)
      case (8)
        open(unit=iunit,file='eos.thermo8',&
          status='old',action='read',iostat=ierror)
      case (9)
        open(unit=iunit,file='eos.thermo9',&
          status='old',action='read',iostat=ierror)
      end select choose_file_0

      nadd_max = 0
      if (ierror == 0) then
        idx_ex(1) = 1
        if (iwr == 1) then
          write(*,*)
          write(*,*) ' reading eos table(s) with thermodynamic properties'
        end if
        read(iunit,*,iostat=ierror) m_n,m_p,incl_l
        if ((ierror /= 0).or.(m_n < 0.d00).or.(m_p < 0.d00)) then
          iflag = 1
          ierr = ierr+1
          if (ierr < dim_err) error_msg(ierr) = 27
        end if
        do while (iflag == 0)
          read(iunit,*,iostat=ierror) it,in,iy,(tmp(iv),iv=1,dim_reg,1),&
            nqty
          if (ierror /= 0) then
            iflag = 1
          else
            if (nqty < 0) then
              iflag = 1
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 21
            end if
            if (nqty > nadd_max) nadd_max = nqty
            icnt = icnt+1
          end if
        end do
      end if

      close(unit=iunit)
    end do
    if (iwr == 1) then
      write(*,*)
      write(*,*) ' maximum number of additional quantities in eos.thermo =',nadd_max
    end if

    ! number of stored quantities in eos.thermo
    nall_max = dim_reg+nadd_max
  end if

  if (ii /= 1) then
    allocate(tab_thermo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:nall_max),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 21
    end if
    allocate(idx_thermo(1:nall_max),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 34
    end if
    allocate(v_thermo(1:nall_max,0:4),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 35
    end if
    if (nadd_max < 2) then
      allocate(eos_thermo_add(1:2),stat=alloc_status)
    else
      allocate(eos_thermo_add(1:nadd_max),stat=alloc_status)
    end if
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 36
    end if
    allocate(idx_add(1:nadd_max),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 37
    end if
    allocate(qty(1:nadd_max),stat=alloc_status)
    if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 38
    end if

    ! initialization
    if ((dim_idx(1) > 0).and.(dim_idx(2) > 0).and.(dim_idx(3) > 0).and.(nall_max > 0)) then
      tab_thermo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:nall_max) = 0.d00
    end if

    call write_errors(ierr)

    icnt = 0
    open(unit=20,file='eos.thermo',&
      status='old',action='read',iostat=ierror)
    if (ierror == 0) then
      itmp = 0
    else
      itmp = nbl-1
    end if

    do ibl=0,itmp,1
      ! flag for stop reading file
      iflag = 0
      iunit = 20+ibl

      choose_file: select case (ibl)
      case (0)
        if (itmp > 0) then
          open(unit=iunit,file='eos.thermo0',&
            status='old',action='read',iostat=ierror)
        end if
      case (1)
        open(unit=iunit,file='eos.thermo1',&
          status='old',action='read',iostat=ierror)
      case (2)
        open(unit=iunit,file='eos.thermo2',&
          status='old',action='read',iostat=ierror)
      case (3)
        open(unit=iunit,file='eos.thermo3',&
          status='old',action='read',iostat=ierror)
      case (4)
        open(unit=iunit,file='eos.thermo4',&
          status='old',action='read',iostat=ierror)
      case (5)
        open(unit=iunit,file='eos.thermo5',&
          status='old',action='read',iostat=ierror)
      case (6)
        open(unit=iunit,file='eos.thermo6',&
          status='old',action='read',iostat=ierror)
      case (7)
        open(unit=iunit,file='eos.thermo7',&
          status='old',action='read',iostat=ierror)
      case (8)
        open(unit=iunit,file='eos.thermo8',&
          status='old',action='read',iostat=ierror)
      case (9)
        open(unit=iunit,file='eos.thermo9',&
          status='old',action='read',iostat=ierror)
      end select choose_file

      if (ierror == 0) then
        idx_ex(1) = 1
        read(iunit,*,iostat=ierror) m_n,m_p,incl_l
        do while (iflag == 0)
          read(iunit,*,iostat=ierror) it,in,iy,(tmp(iv),iv=1,dim_reg,1),&
            nqty,(qty(iw),iw=1,nqty,1)
          if (ierror /= 0) then
            iflag = 1
          else
            it2 = it-min_idx(imap(1))+1
            in2 = in-min_idx(imap(2))+1
            iy2 = iy-min_idx(imap(3))+1
            idx_arg(it2,in2,iy2,1) = 1
            do iv=1,dim_reg,1
              tab_thermo(it2,in2,iy2,iv) = tmp(iv)
            end do
            do iw=1,nqty,1
              tab_thermo(it2,in2,iy2,dim_reg+iw) = qty(iw)
            end do
            icnt = icnt+1
          end if
        end do
      end if

      close(unit=iunit)
    end do

    if (allocated(qty)) deallocate(qty)
  end if

  call write_errors(ierr)

  if (incl_l /= 1) then
    ibeta = 0
  else
    ibeta = iyq
  end if

  if (ii == 1) then
    write(unit_in,*) incl_l,nadd_max,nall_max,ibeta
  end if

  if (iwr == 1) then
    write(*,*)
    write(*,*) icnt,'entries of thermodynamic table read'
  end if



 end subroutine read_eos_table_thermo
 !***********************************************************************
 subroutine read_eos_table_compo(iwr,nbl,ii,unit_in)
  ! Stefan Typel for the CompOSE core team, version 2.04, 2016/11/08
  use compose_internal, only : eos_q,  error_msg,idxp_compo, idxq_compo, &
    &   tabp_compo, tabq_compo, idx_q, idx_p, min_idx, idx_ex, dim_idx,&
    &   idx_compo_p,eos_compo_p,idx_compo_q,eos_compo_q, idx_arg,&
    & dim_err, dim_ip, dim_iq, np_max, nq_max, imap
  implicit none
  integer,intent(in) :: iwr,nbl,ii,unit_in
  integer :: ibl,ierror,icnt,it,in,iy,it2,in2,iy2,&
    iphase,iunit,iflag,itmp,ierr,ic,ic_min,ic_max,&
    iv,iw,np,nq,alloc_status,idum0(dim_ip),idum1(dim_iq),&
    idxp(dim_ip),idxq(dim_iq),ip,iq,jflag
  integer, allocatable :: iqtyp(:),iqtyq(:)
  double precision :: dum0,dum1,dum2,dum3
  double precision, allocatable :: qtyp(:),qtyq(:,:)
  integer, parameter :: unit0 = 200

  ! counter for total number of read table entries
  icnt = 0

  ! error counter
  ierr = 0

  !if (ii /= 2) then
  !   np_max = 0
  !   nq_max = 0
  !end if

  if (dim_ip > 0) idxp(1:dim_ip) = -1

  if (dim_iq > 0) idxq(1:dim_iq) = -1

  ! reading cycles
  idx_ex(2) = 0
  ic_min = 1
  ic_max = 2
  
  if (ii == 1) then
    ic_max = 1
  end if
  if (ii == 2) then
    ic_min = 2
  end if
  do ic=ic_min,ic_max,1
    if (ic == 1) then
      np_max = 0
      nq_max = 0
    end if
    open(unit=unit0,file='eos.compo', status='old',action='read',iostat=ierror)
    if (ierror == 0) then
      itmp = 0
    else
      itmp = nbl-1
    end if

    do ibl=0,itmp,1
      ! flag for stop reading file
      iflag = 0
      iunit = unit0 + ibl

      choose_file: select case (ibl)
      case (0)
        if (itmp > 0) then
          open(unit=iunit,file='eos.compo0',&
            status='old',action='read',iostat=ierror)
        end if
      case (1)
        open(unit=iunit,file='eos.compo1',&
          status='old',action='read',iostat=ierror)
      case (2)
        open(unit=iunit,file='eos.compo2',&
          status='old',action='read',iostat=ierror)
      case (3)
        open(unit=iunit,file='eos.compo3',&
          status='old',action='read',iostat=ierror)
      case (4)
        open(unit=iunit,file='eos.compo4',&
          status='old',action='read',iostat=ierror)
      case (5)
        open(unit=iunit,file='eos.compo5',&
          status='old',action='read',iostat=ierror)
      case (6)
        open(unit=iunit,file='eos.compo6',&
          status='old',action='read',iostat=ierror)
      case (7)
        open(unit=iunit,file='eos.compo7',&
          status='old',action='read',iostat=ierror)
      case (8)
        open(unit=iunit,file='eos.compo8',&
          status='old',action='read',iostat=ierror)
      case (9)
        open(unit=iunit,file='eos.compo9',&
          status='old',action='read',iostat=ierror)
      end select choose_file

      if (ierror == 0) then
        idx_ex(2) = 1
        if (ic == 1) then
          if (iwr == 1) then
            write(*,*)
            write(*,*) ' reading eos table(s) with composition'
          end if
          do while (iflag == 0)
            read(iunit,*,iostat=ierror) it,in,iy,iphase,&
              np,(idum0(iv),dum0,iv=1,np,1),&
              nq,(idum1(iw),dum1,dum2,dum3,iw=1,nq,1)
            if (ierror /= 0) then
              iflag = 1
            else
              icnt = icnt+1
              do ip=1,np,1
                jflag = 0
                do iq=1,np_max,1
                  if (idum0(ip) == idxp(iq)) then
                    jflag = 1
                  end if
                end do
                if (jflag == 0) then
                  np_max = np_max+1
                  if (np_max > dim_ip) then
                    ierr = ierr + 1
                    if (ierr < dim_err) error_msg(ierr) = 175
                    call write_errors(ierr)
                  end if
                  idxp(np_max) = idum0(ip)
                end if
              end do
              do ip=1,nq,1
                jflag = 0
                do iq=1,nq_max,1
                  if (idum1(ip) == idxq(iq)) then
                    jflag = 1
                  end if
                end do
                if (jflag == 0) then
                  nq_max = nq_max+1
                  if (nq_max > dim_iq) then
                    ierr = ierr + 1
                    if (ierr < dim_err) error_msg(ierr) = 176
                    call write_errors(ierr)
                  end if
                  idxq(nq_max) = idum1(ip)
                end if
              end do
            end if
          end do

          if ((idx_ex(2) == 1).and.(iwr == 1)) then
            write(*,*)
            write(*,*) ' maximum number of pairs in eos.compo      =',np_max
            write(*,*) ' maximum number of quadruples in eos.compo =',nq_max
          end if
        else
          icnt = 0
          if (ii /= 1) then
            call allocate_int(iqtyp, np_max, 0, alloc_status, ierr, dim_err, error_msg, 102)
            call allocate_int(iqtyq, nq_max, 0, alloc_status, ierr, dim_err, error_msg, 102)

            allocate(qtyp(1:np_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 103
            end if
            allocate(qtyq(1:nq_max,1:3),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 104
            end if
            allocate(idxp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:np_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 105
            end if
            idxp_compo = 0
            allocate(idxq_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:nq_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 106
            end if
            idxq_compo = 0
            allocate(tabp_compo(1:dim_idx(1),1:dim_idx(2),0:dim_idx(3),1:np_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 107
            end if
            tabp_compo = 0
            allocate(tabq_compo(0:dim_idx(1),1:dim_idx(2),0:dim_idx(3),1:3*nq_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 108
            end if
            tabq_compo = 0

            allocate(idx_p(1:np_max),stat=alloc_status)
            call allocate_int(idx_q, nq_max, 0, alloc_status, ierr, dim_err, error_msg, 132)
            call allocate_int(idx_compo_p, np_max, 0, alloc_status, ierr, dim_err, error_msg, 134)
            call allocate_int(idx_compo_q, nq_max, 0, alloc_status, ierr, dim_err, error_msg, 135)


            allocate(eos_compo_p(1:np_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 137
            end if
            eos_compo_p = 0.d00
            allocate(eos_compo_q(1:nq_max,1:4),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 138
            end if
            eos_compo_q = 0.d00
            allocate(eos_q(1:4*nq_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 145
            end if
            eos_q = 0.d00

            ! initialization
            if ((dim_idx(1) > 0).and.(dim_idx(2) > 0).and.(dim_idx(3) > 0)) then
              idx_arg(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),2) = 0
              if (np_max > 0) then
                idxp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                  1:np_max) = -1
                tabp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                  1:np_max) = 0.d00
              end if
              if (nq_max > 0) then
                tabp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                  1:np_max) = 0.d00
                idxq_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                  1:nq_max) = -1
                tabq_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                  1:3*nq_max) = 0.d00
              end if
            end if
            call write_errors(ierr)
            do while (iflag == 0)
              read(iunit,*,iostat=ierror) it,in,iy,iphase,&
                np,(iqtyp(iv),qtyp(iv),iv=1,np,1),&
                nq,(iqtyq(iw),qtyq(iw,1),qtyq(iw,2),qtyq(iw,3),iw=1,nq,1)

              if (ierror /= 0) then
                iflag = 1
              else
                it2 = it-min_idx(imap(1))+1
                in2 = in-min_idx(imap(2))+1
                iy2 = iy-min_idx(imap(3))+1
                idx_arg(it2,in2,iy2,2) = 1
                idx_arg(it2,in2,iy2,4) = iphase
                do iv=1,np,1
                  idxp_compo(it2,in2,iy2,iv) = iqtyp(iv)
                  tabp_compo(it2,in2,iy2,iv) = qtyp(iv)
                end do
                do iw=1,nq,1
                  idxq_compo(it2,in2,iy2,iw) = iqtyq(iw)
                  tabq_compo(it2,in2,iy2,3*(iw-1)+1) = qtyq(iw,1)
                  tabq_compo(it2,in2,iy2,3*(iw-1)+2) = qtyq(iw,2)
                  tabq_compo(it2,in2,iy2,3*(iw-1)+3) = qtyq(iw,3)
                end do
                icnt = icnt+1
              end if
            end do
          end if
          if ((idx_ex(2) == 1).and.(iwr == 1)) then
            write(*,*)
            write(*,*) icnt,'entries of composition table read'
          end if
        end if
      else
        if ((ic == 1).and.(ibl == 0)) then
          call allocate_int(iqtyp, np_max, 0, alloc_status, ierr, dim_err, error_msg, 101)
          call allocate_int(iqtyq, nq_max, 0, alloc_status, ierr, dim_err, error_msg, 102)
          
          allocate(qtyp(1:np_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 103
          end if
          qtyp = 0.d00
          allocate(qtyq(1:nq_max,1:3),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 104
          end if
          qtyq = 0.d00
          allocate(idxp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:np_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 105
          end if
          idxp_compo = 0
          allocate(idxq_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:nq_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 106
          end if
          idxq_compo = 0
          allocate(tabp_compo(1:dim_idx(1),1:dim_idx(2),0:dim_idx(3),1:np_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 107
          end if
          tabp_compo = 0.d00
          allocate(tabq_compo(0:dim_idx(1),1:dim_idx(2),0:dim_idx(3),1:3*nq_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 108
          end if
          tabq_compo = 0.d00
          
          call allocate_int(idx_p, np_max, 0, alloc_status, ierr, dim_err, error_msg, 131)
          call allocate_int(idx_q, nq_max, 0, alloc_status, ierr, dim_err, error_msg, 132)
          call allocate_int(idx_compo_p, np_max, 0, alloc_status, ierr, dim_err, error_msg, 134)
          call allocate_int(idx_compo_q, nq_max, 0, alloc_status, ierr, dim_err, error_msg, 135)

          allocate(eos_compo_p(1:np_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 137
          end if
          eos_compo_p = 0.d00
          allocate(eos_compo_q(1:nq_max,1:4),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 138
          end if
          eos_compo_q = 0.d00
          allocate(eos_q(1:4*nq_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 145
          end if
          eos_q = 0.d00
          ! initialization

          if ((dim_idx(1) > 0).and.(dim_idx(2) > 0).and.(dim_idx(3) > 0)) then
            idx_arg(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),2) = 0
            if (np_max > 0) then
              idxp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                1:np_max) = -1
              tabp_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                1:np_max) = 0.d00
            end if
            if (nq_max > 0) then
              idxq_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                1:nq_max) = -1
              tabq_compo(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                1:3*nq_max) = 0.d00
            end if
          end if
          call write_errors(ierr)
        end if
      end if
      close(unit=iunit)
    end do
  end do

  call write_errors(ierr)

  if (ii == 1) then
    write(unit_in,*) np_max,nq_max
    write(unit_in,*) (idxp(iv),iv=1,np_max,1)
    write(unit_in,*) (idxq(iv),iv=1,nq_max,1)
  end if

  if (ii /= 1) then
    if (allocated(iqtyp)) deallocate(iqtyp)
    if (allocated(iqtyq))deallocate(iqtyq)
    if (allocated(qtyp)) deallocate(qtyp)
    if (allocated(qtyq)) deallocate(qtyq)
  end if


  return
 end subroutine read_eos_table_compo


 !***********************************************************************
 subroutine read_eos_table_micro(iwr,nbl,ii,unit_in)
  ! Stefan Typel for the CompOSE core team, version 2.03, 2016/11/08
  USE compose_internal
  implicit none
  integer,intent(in) :: iwr,nbl,ii,unit_in
  integer :: ibl,ierror,icnt,it,in,iy,it2,in2,iy2,&
    iunit,iflag,itmp,ierr,ic,ic_min,ic_max,iv,nm,alloc_status,&
    idum(dim_im),im,jflag,idxm(dim_im)
  integer, allocatable :: iqtym(:)
  double precision :: dum
  double precision, allocatable :: qtym(:)
  integer, parameter :: unit0 = 400

  !nm_max = 0

  if (dim_im > 0) idxm(1:dim_im) = -1

  ! counter for total number of read table entries
  icnt = 0

  ! error counter
  ierr = 0

  ! reading cycles
  idx_ex(3) = 0
  ic_min = 1
  ic_max = 2
  if (ii == 1) then
    ic_max = 1
  end if
  if (ii == 2) then
    ic_min = 2
  end if
  do ic=ic_min,ic_max,1
    if (ic == 1) nm_max = 0
    open(unit=unit0,file='eos.micro',&
      status='old',action='read',iostat=ierror)
    if (ierror == 0) then
      itmp = 0
    else
      itmp = nbl-1
    end if

    do ibl=0,itmp,1
      ! flag for stop reading file
      iflag = 0
      iunit = unit0 + ibl

      choose_file: select case (ibl)
      case (0)
        if (itmp > 0) then
          open(unit=iunit,file='eos.micro0',&
            status='old',action='read',iostat=ierror)
        end if
      case (1)
        open(unit=iunit,file='eos.micro1',&
          status='old',action='read',iostat=ierror)
      case (2)
        open(unit=iunit,file='eos.micro2',&
          status='old',action='read',iostat=ierror)
      case (3)
        open(unit=iunit,file='eos.micro3',&
          status='old',action='read',iostat=ierror)
      case (4)
        open(unit=iunit,file='eos.micro4',&
          status='old',action='read',iostat=ierror)
      case (5)
        open(unit=iunit,file='eos.micro5',&
          status='old',action='read',iostat=ierror)
      case (6)
        open(unit=iunit,file='eos.micro6',&
          status='old',action='read',iostat=ierror)
      case (7)
        open(unit=iunit,file='eos.micro7',&
          status='old',action='read',iostat=ierror)
      case (8)
        open(unit=iunit,file='eos.micro8',&
          status='old',action='read',iostat=ierror)
      case (9)
        open(unit=iunit,file='eos.micro9',&
          status='old',action='read',iostat=ierror)
      end select choose_file

      if (ierror == 0) then
        idx_ex(3) = 1
        if (ic == 1) then
          if (iwr == 1) then
            write(*,*)
            write(*,*) ' reading eos table(s) with microscopic information'
          end if

          do while (iflag == 0)
            read(iunit,*,iostat=ierror) it,in,iy,&
              nm,(idum(iv),dum,iv=1,nm,1)
            if (ierror /= 0) then
              iflag = 1
            else
              icnt = icnt+1
              do im=1,nm,1
                jflag = 0
                do iv=1,nm_max,1
                  if (idum(im) == idxm(iv)) then
                    jflag = 1
                  end if
                end do
                if (jflag == 0) then
                  nm_max = nm_max+1
                  if (nm_max > dim_im) then
                    ierr = ierr + 1
                    if (ierr < dim_err) error_msg(ierr) = 177
                    call write_errors(ierr)
                  end if
                  idxm(nm_max) = idum(im)
                end if
              end do
            end if
          end do
          if ((idx_ex(3) == 1).and.(iwr == 1)) then
            write(*,*)
            write(*,*) ' maximum number of pairs in eos.micro      =',nm_max
          end if
        else
          icnt = 0
          if (ii /= 1) then
            call allocate_int(iqtym, nm_max, 0, alloc_status, ierr, dim_err, error_msg, 109)
            
            allocate(qtym(1:nm_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 110
            end if
            qtym = 0.d00
            allocate(idx_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),nm_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 111
            end if
            idx_mic = 0

            allocate(tab_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),nm_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 112
            end if
            tab_mic = 0.d00

            call allocate_int(idx_m, nm_max, 0, alloc_status, ierr, dim_err, error_msg, 133)
            call allocate_int(idx_micro, nm_max, 0, alloc_status, ierr, dim_err, error_msg, 136)

            
            allocate(eos_micro(1:nm_max),stat=alloc_status)
            if (alloc_status /= 0) then
              ierr = ierr+1
              if (ierr < dim_err) error_msg(ierr) = 139
            end if
            eos_micro = 0.d00

            ! initialization
            if ((dim_idx(1) > 0).and.(dim_idx(2) > 0).and.(dim_idx(3) > 0)) then
              idx_arg(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),3) = 0
              if (nm_max > 0) then
                idx_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:nm_max) = -1
                tab_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:nm_max) = 0.d00
              end if
            end if
          end if
          !         else
          do while (iflag == 0)
            read(iunit,*,iostat=ierror) it,in,iy,&
              nm,(iqtym(iv),qtym(iv),iv=1,nm,1)
            if (ierror /= 0) then
              iflag = 1
            else
              it2 = it-min_idx(imap(1))+1
              in2 = in-min_idx(imap(2))+1
              iy2 = iy-min_idx(imap(3))+1
              idx_arg(it2,in2,iy2,3) = 1
              do iv=1,nm,1
                idx_mic(it2,in2,iy2,iv) = iqtym(iv)
                tab_mic(it2,in2,iy2,iv) = qtym(iv)
              end do
              icnt = icnt+1
            end if
          end do
          if ((idx_ex(3) == 1).and.(iwr == 1)) then
            write(*,*)
            write(*,*) icnt,'entries of microscopic table read'
          end if
        end if
      else
        if ((ic == 1).and.(ibl == 0)) then
          call allocate_int(iqtym, nm_max, 0, alloc_status, ierr, dim_err, error_msg, 109)
          
          allocate(qtym(1:nm_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 110
          end if
          qtym = 0.d00
          allocate(idx_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),nm_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 111
          end if
          idx_mic = 0
          allocate(tab_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),nm_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 112
          end if
          tab_mic = 0.d00
          call allocate_int(idx_m, nm_max, 0, alloc_status, ierr, dim_err, error_msg, 133)
          call allocate_int(idx_micro, nm_max, 0, alloc_status, ierr, dim_err, error_msg, 136)

          allocate(eos_micro(1:nm_max),stat=alloc_status)
          if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 139
          end if
          eos_micro = 0.d00

          ! initialization

          if ((dim_idx(1) > 0).and.(dim_idx(2) > 0).and.(dim_idx(3) > 0)) then
            idx_arg(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),3) = 0
            if (nm_max > 0) then
              idx_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                1:nm_max) = -1
              tab_mic(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),&
                1:nm_max) = 0.d00
            end if
          end if
        end if
      end if
    end do
    close(unit=iunit)
  end do
  call write_errors(ierr)

  if (ii == 1) then
    write(unit_in,*) nm_max
    write(unit_in,*) (idxm(iv),iv=1,nm_max,1)
  end if

  if (ii /= 1) then
    if (allocated(iqtym)) deallocate(iqtym)
    if (allocated(qtym)) deallocate(qtym)
  end if


 end subroutine read_eos_table_micro

 subroutine allocate_int(array, size_array, value_init, alloc_status, ierr, dim_err, error_msg, error_code)
  integer, allocatable, intent(inout) :: array(:)
  integer, intent(in) :: size_array, value_init, dim_err, error_code
  integer, intent(out) :: alloc_status
  integer, intent(inout) :: ierr
  integer, intent(inout) :: error_msg(:)

  allocate(array(1:size_array), stat=alloc_status)
  if (alloc_status /= 0) then
    ierr = ierr + 1
    if (ierr < dim_err) error_msg(ierr) = error_code
  else
    array = value_init
  end if
end subroutine allocate_int


end module m_get_tables
