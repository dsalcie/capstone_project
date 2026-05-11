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


! Sample routine for reading the EoS data from a Compose generated table in
! Hdf5 format

subroutine read_hdf5

  use hdf5
  use modhdf5
  use precision
  IMPLICIT NONE


  real(dp), allocatable :: nb_hdf5(:),t_hdf5(:),y_q_hdf5(:)
  real(dp), allocatable :: thermo_hdf5(:,:,:,:), thermo_hdf5_add(:,:,:,:)
  real(dp), allocatable :: yi_hdf5(:,:,:,:),aav_hdf5(:,:,:,:),zav_hdf5(:,:,:,:),yav_hdf5(:,:,:,:)
  real(dp), allocatable :: nav_hdf5(:,:,:,:)
  real(dp), allocatable :: micro_hdf5(:,:,:,:),err_hdf5(:,:,:,:)
  integer,allocatable :: index_thermo(:),index_thermo_add(:),index_err(:)
  integer, allocatable :: index_yi(:),index_av(:),index_micro(:)

  integer itab

  integer m_nb,n_temp,o_y_q,n_qty,n_add,n_p,n_q,n_m,n_err,i,j,k
  integer :: n4d(4)

  character(len=500) :: file_read

  integer(hid_t) ::  h5id, h5file_read
  character(LEN=1),parameter :: as = '"'
  logical :: it_exist,group_exist


  ! First read the number of different quantities stored


  write(*,*) 'File name for reading data?'
  read(*,*) file_read


  inquire(file=file_read,exist=it_exist)
  if(.not. it_exist) then
     print '(3a)', "Error : the file '",trim(file_read),"' does not exist"
     return
  endif

  call hdf5_init()

  call hdf5_open_file(file_read,h5file_read,'r')

  ! reading thermo parameters and array dimensions

  call hdf5_open_group(h5file_read,'Parameters', h5id)
  call hdf5_read_attr(h5id,'pointsnb',m_nb)
  call hdf5_read_attr(h5id,'tabulation_scheme',itab)
  call hdf5_read_attr(h5id,'pointst',n_temp)
  call hdf5_read_attr(h5id,'pointsyq',o_y_q)
  allocate(nb_hdf5(m_nb))
  call hdf5_read_data(h5id,'nb',m_nb,nb_hdf5)

  if(itab.eq.0) then
     allocate(t_hdf5(m_nb))
     allocate(y_q_hdf5(m_nb))
     call hdf5_read_data(h5id,'t',m_nb,t_hdf5)
     call hdf5_read_data(h5id,'yq',m_nb,y_q_hdf5)
  else
     allocate(t_hdf5(n_temp))
     allocate(y_q_hdf5(o_y_q))
     call hdf5_read_data(h5id,'t',n_temp,t_hdf5)
     call hdf5_read_data(h5id,'yq',o_y_q,y_q_hdf5)
  end if
  call hdf5_close_group(h5id)

  n4d(1) = m_nb
  n4d(2) = n_temp
  n4d(3) = o_y_q
! read thermo quantities

  call hdf5_query_group(h5file_read, 'Thermo_qty', group_exist)
  if(.not. group_exist) then
     n_qty = 0
     write(*,*) 'No thermo quantities in file'
  else
     call hdf5_open_group(h5file_read,'Thermo_qty', h5id)
     call hdf5_read_attr(h5id,'pointsqty',n_qty)
     write(*,*) 'Reading', n_qty,'thermo quantities from file'
     n4d(4) = n_qty
     allocate(thermo_hdf5(m_nb,n_temp,o_y_q,n_qty))
     allocate(index_thermo(n_qty))
     call hdf5_read_data(h5id,'thermo',n4d,thermo_hdf5)
     call hdf5_read_data(h5id,'index_thermo',n_qty,index_thermo)
     call hdf5_close_group(h5id)
     write(*,*) 'thermo',thermo_hdf5(1,1,1,1),index_thermo(1)

  end if

! read additional thermo quantities

  call hdf5_query_group(h5file_read, 'Thermo_add', group_exist)
  if(.not. group_exist) then
     n_add = 0
     write(*,*) 'No additonal thermo quantities in file'
  else
     call hdf5_open_group(h5file_read,'Thermo_add', h5id)
     call hdf5_read_attr(h5id,'pointsadd',n_add)
     write(*,*) 'Reading', n_add,'additional thermo quantities from file'
     n4d(4) = n_add
     allocate(thermo_hdf5_add(m_nb,n_temp,o_y_q,n_add))
     allocate(index_thermo_add(n_add))
     call hdf5_read_data(h5id,'thermo_add',n4d,thermo_hdf5_add)
     call hdf5_read_data(h5id,'index_thermo_add',n_add,index_thermo_add)
     call hdf5_close_group(h5id)
     write(*,*) 'error',thermo_hdf5_add(1,1,1,1),index_thermo_add(1)

  end if

  ! read pairs for compositional information


  call hdf5_query_group(h5file_read, 'Composition_pairs', group_exist)
  if(.not. group_exist) then
     n_p = 0
     write(*,*) 'No pairs in file'
  else

     call hdf5_open_group(h5file_read,'Composition_pairs', h5id)
     call hdf5_read_attr(h5id,'pointspairs',n_p)
     n4d(4) = n_p
     write(*,*) 'Reading', n_p,'pairs from file'
     allocate(yi_hdf5(m_nb,n_temp,o_y_q,n_p))
     allocate(index_yi(n_p))
     call hdf5_read_data(h5id,'yi',n4d,yi_hdf5)
     call hdf5_read_data(h5id,'index_yi',n_p,index_yi)
     call hdf5_close_group(h5id)
     write(*,*) 'P',yi_hdf5(1,1,1,1),index_yi(1)

  end if

  ! read quadrupels for compositional information


  call hdf5_query_group(h5file_read, 'Composition_quadrupels', group_exist)
  if(.not. group_exist) then
     n_q = 0
     write(*,*) 'No quadrupels in file'
  else
     call hdf5_open_group(h5file_read,'Composition_quadrupels', h5id)
     call hdf5_read_attr(h5id,'pointsav',n_q)
     n4d(4) = n_q
     write(*,*) 'Reading', n_q,'quadrupels from file'
     allocate(yav_hdf5(m_nb,n_temp,o_y_q,n_q))
     allocate(zav_hdf5(m_nb,n_temp,o_y_q,n_q))
     allocate(aav_hdf5(m_nb,n_temp,o_y_q,n_q))
     allocate(nav_hdf5(m_nb,n_temp,o_y_q,n_q))
     allocate(index_av(n_q))

     call hdf5_read_data(h5id,'yav',n4d,yav_hdf5)
     call hdf5_read_data(h5id,'aav',n4d,aav_hdf5)
     call hdf5_read_data(h5id,'zav',n4d,zav_hdf5)
     call hdf5_read_data(h5id,'nav',n4d,nav_hdf5)
     call hdf5_read_data(h5id,'index_av',n_q,index_av)
     call hdf5_close_group(h5id)
     write(*,*) 'Q',yav_hdf5(1,1,1,1),index_av(1)

  end if

  ! read microscopic information


  call hdf5_query_group(h5file_read, 'Micro_qty', group_exist)
  if(.not. group_exist) then
     n_m = 0
     write(*,*) 'No microscopic information in file'
  else
     call hdf5_open_group(h5file_read,'Micro_qty', h5id)
     call hdf5_read_attr(h5id,'pointsmicro',n_m)
     n4d(4) = n_m
     write(*,*) 'Reading', n_m,'microscopic quantities from file'
     allocate(micro_hdf5(m_nb,n_temp,o_y_q,n_m))
     allocate(index_micro(n_m))
     call hdf5_read_data(h5id,'micro',n4d,micro_hdf5)
     call hdf5_read_data(h5id,'index_micro',n_q,index_micro)
     call hdf5_close_group(h5id)
     write(*,*) 'micro',micro_hdf5(1,1,1,1),index_micro(1)
  end if


  ! read error information


  call hdf5_query_group(h5file_read, 'Error_qty', group_exist)
  if(.not. group_exist) then
     n_err = 0
     write(*,*) 'No error information in file'
  else

     call hdf5_open_group(h5file_read,'Error_qty', h5id)
     call hdf5_read_attr(h5id,'pointserr',n_err)
     n4d(4) = n_err
     write(*,*) 'Reading', n_err,'error quantities from file'
     allocate(err_hdf5(m_nb,n_temp,o_y_q,n_err))
     allocate(index_err(n_m))
     call hdf5_read_data(h5id,'error',n4d,err_hdf5)
     call hdf5_read_data(h5id,'index_err',n_err,index_err)
     call hdf5_close_group(h5id)
     write(*,*) 'error',err_hdf5(1,1,1,1),index_err(1)
  end if

  call hdf5_close_file(h5file_read)

  ! write the data in ascii format for checking purposes
  open(33,file='readtest.d',status='unknown')
  IF(itab.eq.0) then
     DO i = 1,m_nb
        write(33,*) nb_hdf5(i),t_hdf5(i),y_q_hdf5(i)
     end DO
  else
     DO i = 1,n_temp
        DO j = 1,m_nb
           DO k = 1,o_y_q
              write(33,*) t_hdf5(i),nb_hdf5(j),y_q_hdf5(k)
           end DO
        end DO
     end DO
  end IF


  IF(allocated(thermo_hdf5)) deallocate(thermo_hdf5)
  IF(allocated(thermo_hdf5_add)) deallocate(thermo_hdf5_add)
  IF(allocated(yi_hdf5)) deallocate(yi_hdf5)
  IF(allocated(yav_hdf5)) deallocate(yav_hdf5)
  IF(allocated(zav_hdf5)) deallocate(zav_hdf5)
  IF(allocated(nav_hdf5)) deallocate(nav_hdf5)
  IF(allocated(aav_hdf5)) deallocate(aav_hdf5)
  IF(allocated(micro_hdf5)) deallocate(micro_hdf5)
  IF(allocated(err_hdf5)) deallocate(err_hdf5)
  IF(allocated(index_thermo)) deallocate(index_thermo)
  IF(allocated(index_thermo_add)) deallocate(index_thermo_add)
  IF(allocated(index_err)) deallocate(index_err)
  IF(allocated(index_yi)) deallocate(index_yi)
  IF(allocated(index_av)) deallocate(index_av)
  IF(allocated(index_micro)) deallocate(index_micro)
  IF(allocated(nb_hdf5)) deallocate(nb_hdf5)
  IF(allocated(t_hdf5)) deallocate(t_hdf5)
  IF(allocated(y_q_hdf5)) deallocate(y_q_hdf5)


end subroutine read_hdf5


program read_eos

  write(*,*) 'Sample program to read data from the EoS table in HDF5 format '
  write(*,*) 'as generated by the Compose software. Please see the manual '
  write(*,*) 'for the definition of the different quantities '
  call read_hdf5
end program read_eos
