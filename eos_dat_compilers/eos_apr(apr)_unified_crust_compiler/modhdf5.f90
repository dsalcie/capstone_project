!
!   Copyright (c) 2013-2022 Fabrice Roy
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

!> @file
!! This file contains wrappers for HDF5 serial I/O functions.

!> Authors: F. Roy

module modhdf5

  use hdf5
  use iso_c_binding

 implicit none

 integer(kind=4), parameter :: H5STRLEN = 32

 private
 integer(kind=4) :: h5err !< hdf5 error code


 public :: hdf5_open_file, &
   hdf5_create_file, &
   hdf5_close_file, &
   hdf5_open_group, &
   hdf5_create_group, &
   hdf5_close_group, &
   hdf5_find_groups, &
   hdf5_query_group,       &
   hdf5_query_attr,    &
   hdf5_write_data, &
   hdf5_read_data, &
   hdf5_write_attr, &
   hdf5_read_attr, &
   hdf5_init, &
   hdf5_finalize

 public :: hid_t  !< integer kind for hdf5 id, comes from hdf5 module and must be public to be used outside this module.
 public :: H5STRLEN

 !> Generic inteface used to write attributes in a hdf5 file.
 !> @param[in] id id of the file/group where the attribute will be written
 !> @param[in] name name of the attribute
 !> @param[in] n1 first dimension of the attribute array (optional)
 !> @param[in] n2 second dimension of the attribute array (optional)
 !> @param[in] data 0-D, 1-D or 2-D attribute of type Integer(4), Real(4), Real(8) or Character
 interface hdf5_write_attr
   module procedure hdf5_write_int4_attr0D   ! (id, name, data)
   module procedure hdf5_write_int4_attr1D   ! (id, name, n1, data)
   module procedure hdf5_write_int4_attr2D   ! (id, name, n1, n2, data)
   module procedure hdf5_write_int4_attr3D   ! (id, name, n(3), data)
   module procedure hdf5_write_real4_attr0D  ! (id, name, data)
   module procedure hdf5_write_real4_attr1D  ! (id, name, n1, data)
   module procedure hdf5_write_real4_attr2D  ! (id, name, n1, n2, data)
   module procedure hdf5_write_real8_attr0D  ! (id, name, data)
   module procedure hdf5_write_real8_attr1D  ! (id, name, n1, data)
   module procedure hdf5_write_real8_attr2D  ! (id, name, n1, n2, data)
   module procedure hdf5_write_char_attr     ! (id, name, data)
   module procedure hdf5_write_char_attr1D   ! (id, name, n1, data)
 end interface hdf5_write_attr

 !> Generic inteface used to read attributes from a hdf5 file.
 !> @param[in] id id of the file/group where the attribute will be read
 !> @param[in] name name of the attribute
 !> @param[in] n1 first dimension of the attribute array (optional)
 !> @param[in] n2 second dimension of the attribute array (optional)
 !> @param[in,out] data 0-D, 1-D or 2-D attribute of type Integer(4), Real(4), Real(8) or Character
 interface hdf5_read_attr
   module procedure hdf5_read_int4_attr0D   ! (id, name, data)
   module procedure hdf5_read_int4_attr1D   ! (id, name, n1, data)
   module procedure hdf5_read_int4_attr2D   ! (id, name, n1, n2, data)
   module procedure hdf5_read_int4_attr3D   ! (id, name, n(3) data)
   module procedure hdf5_read_real4_attr0D  ! (id, name, data)
   module procedure hdf5_read_real4_attr1D  ! (id, name, n1, data)
   module procedure hdf5_read_real4_attr2D  ! (id, name, n1, n2, data)
   module procedure hdf5_read_real8_attr0D  ! (id, name, data)
   module procedure hdf5_read_real8_attr1D  ! (id, name, n1, data)
   module procedure hdf5_read_real8_attr2D  ! (id, name, n1, n2, data)
   module procedure hdf5_read_char_attr     ! (id, name, data)
   module procedure hdf5_read_char_attr1D   ! (id, name, len , n1, data)
 end interface hdf5_read_attr

 !> Generic inteface used to write data in a hdf5 file.
 !> @param[in] id id of the file/group where the dataset will be written
 !> @param[in] name name of the dataset
 !> @param[in] n1 first dimension of the data array
 !> @param[in] n2 second dimension of the data array (optional)
 !> @param[in] data 1-D or 2-D data array of type Integer(4), Integer(8), Real(4) or Real(8)
 interface hdf5_write_data
   module procedure hdf5_write_int4_1D       ! (id, name, n1, data)
   module procedure hdf5_write_int4_2D       ! (id, name, n1, n2, data)
   module procedure hdf5_write_int4_3D       ! (id, name, n, data)
   module procedure hdf5_write_int8_0D       ! (id, name, data)
   module procedure hdf5_write_int8_1D       ! (id, name, n1, data)
   module procedure hdf5_write_int8_2D       ! (id, name, n1, n2, data)
   module procedure hdf5_write_real4_1D      ! (id, name, n1, data)
   module procedure hdf5_write_real4_2D      ! (id, name, n1, n2, data)
   module procedure hdf5_write_real8_1D      ! (id, name, n1, data)
   module procedure hdf5_write_real8_2D      ! (id, name, n1, n2, data)
   module procedure hdf5_write_real8_3D      ! (id, name, n(3), data)
   module procedure hdf5_write_real8_4D      ! (id, name, n(4), data)
   module procedure hdf5_write_real8_5D      ! (id, name, n(5), data)
 end interface hdf5_write_data

 !> Generic inteface used to read data from a hdf5 file.
 !> @param[in] id id of the file/group where the dataset will be read
 !> @param[in] name name of the dataset
 !> @param[in] n1 first dimension of the data array
 !> @param[in] n2 second dimension of the data array (optional)
 !> @param[in,out] data 1-D or 2-D data array of type Integer(4), Integer(8), Real(4) or Real(8)
 interface hdf5_read_data
   module procedure hdf5_read_int4_1D        ! (id, name, n1, data)
   module procedure hdf5_read_int4_2D        ! (id, name, n1, n2, data)
   module procedure hdf5_read_int4_3D        ! (id, name, n1, n2, data)
   module procedure hdf5_read_int8_0D        ! (id, name, data)
   module procedure hdf5_read_int8_1D        ! (id, name, n1, data)
   module procedure hdf5_read_int8_2D        ! (id, name, n1, n2, data)
   module procedure hdf5_read_real4_1D       ! (id, name, n1, data)
   module procedure hdf5_read_real4_2D       ! (id, name, n1, n2, data)
   module procedure hdf5_read_real4_3D       ! (id, name, n1, n2, data)
   module procedure hdf5_read_real8_1D       ! (id, name, n1, data)
   module procedure hdf5_read_real8_2D       ! (id, name, n1, n2, data)
   module procedure hdf5_read_real8_3D       ! (id, name, n(3), data)
   module procedure hdf5_read_real8_4D       ! (id, name, n(4), data)
   module procedure hdf5_read_real8_5D       ! (id, name, n(5), data)
 end interface hdf5_read_data


contains

 !===============================================!
 !         INIT/CLOSE ROUTINES
 !===============================================!

 !=======================================================================
 !> Open the hdf5 interface
 subroutine hdf5_init()

  implicit none

  call h5open_f(h5err)

#ifdef DEBUG
  if(h5err /= 0) then
    print *,'HDF5 error while initializing the HDF5 interface'
    stop
  end if
#endif

 end subroutine hdf5_init


 !=======================================================================
 !> Close the hdf5 interface
 subroutine hdf5_finalize()

  implicit none

  call h5close_f(h5err)

#ifdef DEBUG
  if(h5err /= 0) then
    print *,'HDF5 error while closing the HDF5 interface'
    stop
  end if
#endif

 end subroutine hdf5_finalize


 !===============================================!
 !         FILE MANAGEMENT ROUTINES
 !===============================================!

 !=======================================================================
 !> Open an existing hdf5 file
 subroutine hdf5_open_file(filename, file_id, rw)

  character(len=*), intent(in) :: filename           !< name of the file to open
  integer(kind=hid_t), intent(out) :: file_id                      !< hdf5 id of the opened file
  character(len=1), intent(in), optional :: rw                     !< read/write selector

  integer :: h5err                                                 ! hdf5 error code
  character(len=1) :: rw_loc


#ifdef DEBUG
  write(*,*) 'Hdf5_open_file begins for file ', trim (filename)
#endif

  if(.not. present(rw)) then
    rw_loc = 'r'
  else
    rw_loc = rw
  end if

  ! open the file
  select case(rw_loc)
  case('r')
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err)
  case('w')
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)
  case default
    write(*,*) 'HDF5 error in opening option for file ',trim(filename)
  end select

#ifdef DEBUG
  if(h5err /= 0) then
    write(*,*) 'HDF5 error while opening file ',trim(filename)
    stop
  end if
#endif

#ifdef DEBUG
  write(*,*) 'Hdf5_open_file ends for file ', trim (filename)
#endif

 end subroutine Hdf5_open_file



 !=======================================================================
 !> Create and open a hdf5 file
 subroutine hdf5_create_file(filename, file_id)

  implicit none

  character(len=*), intent(in) :: filename           !< name of the file to create
  integer(hid_t), intent(out) :: file_id              !< hdf5 identifier of the file

  integer :: rank                                     ! Nb of dimensions for the attributes
  integer(hsize_t), dimension(1) :: adims             ! Dimensions of the attributes
  integer(hid_t) :: gr_id                        ! Root group identifier
  character(len=H5STRLEN) :: aname                          ! Attribute name
  character(len=8) :: simdate                         ! date of creation of the file
  character(len=10) :: simtime                        ! time of creation of the file
  character(len=10) :: attrdate
  character(len=12) :: attrtime

#ifdef DEBUG
  print *,'hdf5_create_file: ',trim(filename)
#endif

  ! Create a new file using default properties.
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)

#ifdef DEBUG
  if(h5err < 0) then
    print *,'HDF5 error while creating the HDF5 file ',filename
    stop
  end if
#endif

  aname = 'metadata'
  ! Open the '/' group to add some info as attributes
  call hdf5_create_group(file_id, aname, gr_id)

  rank = 1
  adims(1) = 1

  ! Get date and time
  call date_and_time(date=simdate,time=simtime)

  attrdate(1:4) = simdate(1:4)
  attrdate(5:5) = '/'
  attrdate(6:7) = simdate(5:6)
  attrdate(8:8) = '/'
  attrdate(9:10) = simdate(7:8)

  attrtime(1:2) = simtime(1:2)
  attrtime(3:3) = 'h'
  attrtime(4:5) = simtime(3:4)
  attrtime(6:6) = 'm'
  attrtime(7:8) = simtime(5:6)
  attrtime(9:9) = 's'
  attrtime(10:12) = simtime(8:10)


  ! Write date as attribute
  aname = 'date'
  call hdf5_write_attr(gr_id, aname, attrdate)

  ! Write time as attribute
  aname = 'time'
  call hdf5_write_attr(gr_id, aname, attrtime)

  ! Close the group.
  call h5gclose_f(gr_id, h5err)

 end subroutine hdf5_create_file


 !=======================================================================
 !> Close a h5 file
 subroutine hdf5_close_file(file_id)

  implicit none

  integer(hid_t), intent(in) :: file_id              !< hdf5 identifier of the file to close

#ifdef DEBUG
  print *,'hdf5_close_file'
#endif

  ! Terminate access to the file.
  call h5fclose_f(file_id, h5err)

#ifdef DEBUG
  if(h5err /= 0) then
    print *,'HDF5 error while closing the file'
    stop
  end if
#endif

 end subroutine hdf5_close_file


 !===============================================!
 !         GROUP MANAGEMENT ROUTINES
 !===============================================!

 !=======================================================================
 !> Open an existing hdf5 group
 subroutine hdf5_open_group(id, name, group_id)

  implicit none

  integer(hid_t), intent(in) :: id                        !< "parent" identifier
  character(len=*), intent(in) :: name             !< group name
  integer(hid_t), intent(out) :: group_id                 !< group identifier

  ! Open the group
  call h5gopen_f(id, name, group_id, h5err)

#ifdef DEBUG
  if(h5err < 0) then
    print *,'HDF5 error while opening the group ',trim(name)
    stop
  end if
#endif

 end subroutine hdf5_open_group


 !=======================================================================
 !> Create and open a hdf5 group
 subroutine hdf5_create_group(id, name, group_id)

  implicit none

  integer(hid_t), intent(in) :: id                        !< "parent" identifier
  character(len=*), intent(in) :: name             !< group name
  integer(hid_t), intent(out) :: group_id                 !< group identifier

  ! Open the '/' group to add some info as attributes
  call h5gcreate_f(id, name, group_id, h5err)

#ifdef DEBUG
  if(h5err < 0) then
    print *,'HDF5 error while creating the group ',trim(name)
    stop
  end if
#endif

 end subroutine hdf5_create_group

 !=======================================================================
 !> Test element
 subroutine hdf5_query_group(id, name, link_exists)
  integer(hid_t), intent(in) :: id                        !< "parent" identifier
  character(len=*), intent(in) :: name             !< group name
  logical, intent(out) :: link_exists                     !< group identifier

  integer :: hdferr                                       !< Error code:

  if(name /= '/') then
    call h5lexists_f(id, name, link_exists, hdferr)
  else
    link_exists = .true.
  endif
 end subroutine hdf5_query_group


 !=======================================================================
 !> Test attribute
 subroutine hdf5_query_attr(id, name, link_exists)
  integer(hid_t), intent(in) :: id                        !< "parent" identifier
  character(len=*), intent(in) :: name             !< group name
  logical, intent(out) :: link_exists                     !< group identifier

  integer :: hdferr                                       !< Error code:

  call h5aexists_f(id, name, link_exists, hdferr)

 end subroutine hdf5_query_attr


 !=======================================================================
 !> Find childrens of a group
 subroutine hdf5_find_groups(file_id, name, subgroups_name,subgroups_type)
  integer(hid_t), intent(in) :: file_id            !< "file" identifier
  character(len=*), intent(in) :: name             !< group name
  character(len=H5STRLEN),intent(inout),allocatable :: subgroups_name(:)
  integer,intent(inout),allocatable,optional :: subgroups_type(:)

  character(len=H5STRLEN) :: name_buffer ! Buffer to hold object's name
  integer :: type ! Type of the object
  integer :: hdferr,nmembers,i            !< Error code:
  logical :: typeout, link_exists

#ifdef DEBUG
  print *,'Enter hdf5_find_groups for group ',name
#endif

  nmembers = 0
  typeout = present(subgroups_type)

  if(allocated(subgroups_name)) deallocate(subgroups_name)
  if(typeout .and. allocated(subgroups_type)) deallocate(subgroups_type)

  call hdf5_query_group(file_id, name, link_exists)
  if(link_exists)then

    call h5gn_members_f(file_id, name, nmembers, hdferr)
    allocate(subgroups_name(nmembers))

    if(typeout) allocate(subgroups_type(nmembers))
    do i = 0, nmembers - 1
      call h5gget_obj_info_idx_f(file_id, name, i, name_buffer, type, hdferr)
      subgroups_name(i+1) = name_buffer
      if(typeout) subgroups_type(i+1) = type
    end do

  endif

#ifdef DEBUG
  write(*,'(a,a,a4,i0)') ' Number of members of group ',name,' is ',nmembers
  write(*,'(a,*(a))') ' Groups found:', subgroups_name

  print *,'Exit hdf5_find_groups ', name
#endif

 end subroutine hdf5_find_groups



 !=======================================================================
 !> Close a hdf5 group
 subroutine hdf5_close_group(group_id)

  implicit none

  integer(hid_t), intent(in) :: group_id                 !< group identifier

  ! close group
  call h5gclose_f(group_id, h5err)

 end subroutine hdf5_close_group


 !===============================================!
 !         SERIAL DATA OUTPUT ROUTINES
 !===============================================!


 !=======================================================================
 !> Write a 1-D integer(kind=4) array in a serial file
 subroutine hdf5_write_int4_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                       !< name of the dataset
  integer, intent(in) :: n1                                   !< dimension of the data
  integer(kind=4), dimension(n1), intent(in), target :: data  !< data array

  integer, parameter :: rank = 1                   ! dimension of the array = 1
  integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr1D                             ! pointer to the array


  h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)
  ptr1D = c_loc(data(1))
  dim1D(1) = n1
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_int4_1D


 !=======================================================================
 !> Write an integer(kind=8) scalar in a serial file
 subroutine hdf5_write_int8_0D(id, name, data)

  implicit none

  integer(hid_t), intent(in) :: id                 !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name            !< name of the dataset
  integer(kind=8), intent(in), target :: data      !< data scalar

  integer, parameter :: rank = 1                   ! dimension of the array = 1
  integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr1D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)
  ptr1D = c_loc(data)
  dim1D(1) = 1
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_int8_0D


 !=======================================================================
 !> Write a 1-D integer(kind=8) array in a serial file
 subroutine hdf5_write_int8_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                      !< name of the dataset
  integer, intent(in) :: n1                                  !< dimension of the data
  integer(kind=8), dimension(n1), intent(in), target :: data !< data array

  integer, parameter :: rank = 1                   ! dimension of the array = 1
  integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr1D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)
  ptr1D = c_loc(data(1))
  dim1D(1) = n1
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_int8_1D


 !=======================================================================
 !> Write a 2-D integer(kind=4) array in a serial file
 subroutine hdf5_write_int4_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                              !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                         !< name of the dataset
  integer, intent(in) :: n1                                     !< first dimension of the data
  integer, intent(in) :: n2                                     !< second dimension of the data
  integer(kind=4), dimension(n1,n2), intent(in), target :: data !< data array

  integer, parameter :: rank = 2                   ! dimension of the array = 1
  integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr2D                             ! pointer to the array


  h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)
  ptr2D = c_loc(data(1,1))
  dim2D(1) = n1
  dim2D(2) = n2
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_int4_2D


 !=======================================================================
 !> Write a 2-D integer(kind=4) array in a serial file
 subroutine hdf5_write_int4_3D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                              !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                   !< name of the dataset
  integer, intent(in) :: n(3)                                   !< first dimension of the data
  integer(kind=4), dimension(n(1),n(2),n(3)), intent(in), target :: data !< data array

  integer, parameter :: rank = 3                   ! dimension of the array = 1
  integer(hsize_t), dimension(3) :: dim3D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr3D                             ! pointer to the array


  h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)
  ptr3D = c_loc(data(1,1,1))
  dim3D(1) = n(1)
  dim3D(2) = n(2)
  dim3D(3) = n(3)
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim3D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr3D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_int4_3D


 !=======================================================================
 !> Write a 2-D integer(kind=8) array in a serial file
 subroutine hdf5_write_int8_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t),    intent(in) :: id                            !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer,           intent(in) :: n1                            !< first dimension of the data
  integer,           intent(in) :: n2                            !< second dimension of the data
  integer(kind=8), dimension(n1,n2), intent(in), target :: data  !< data array

  integer, parameter :: rank = 2                   ! dimension of the array = 1
  integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr2D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)
  ptr2D = c_loc(data(1,1))
  dim2D(1) = n1
  dim2D(2) = n2
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_int8_2D


 !=======================================================================
 !> Write a 1-D real(kind=4) array in a serial file
 subroutine hdf5_write_real4_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                        !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                   !< name of the dataset
  integer, intent(in) :: n1                               !< dimension of the data
  real(kind=4), dimension(n1), intent(in), target :: data !< data array

  integer, parameter :: rank = 1                   ! dimension of the array = 1
  integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr1D                             ! pointer to the array


  h5_kind = h5kind_to_type(4,H5_REAL_KIND)
  ptr1D = c_loc(data(1))
  dim1D(1) = n1
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real4_1D


 !=======================================================================
 !> Write a 1-D real(kind=8) array in a serial file
 subroutine hdf5_write_real8_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                         !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                    !< name of the dataset
  integer, intent(in) :: n1                                !< dimension of the data
  real(kind=8), dimension(n1), intent(in), target :: data  !< data array

  integer, parameter :: rank = 1                   ! dimension of the array = 1
  integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr1D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_REAL_KIND)
  ptr1D = c_loc(data(1))
  dim1D(1) = n1
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real8_1D


 !=======================================================================
 !> Write a 2-D real(kind=4) array in a serial file
 subroutine hdf5_write_real4_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                      !< name of the dataset
  integer, intent(in) :: n1                                  !< first dimension of the data
  integer, intent(in) :: n2                                  !< second dimension of the data
  real(kind=4), dimension(n1,n2), intent(in), target :: data !< data array

  integer, parameter :: rank = 2                   ! dimension of the array = 1
  integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr2D                             ! pointer to the array


  h5_kind = h5kind_to_type(4,H5_REAL_KIND)
  ptr2D = c_loc(data(1,1))
  dim2D(1) = n1
  dim2D(2) = n2
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real4_2D


 !=======================================================================
 !> Write a 2-D real(kind=8) array in a serial file
 subroutine hdf5_write_real8_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                       !< name of the dataset
  integer, intent(in) :: n1                                   !< first dimension of the data
  integer, intent(in) :: n2                                   !< second dimension of the data
  real(kind=8), dimension(n1,n2), intent(in), target :: data  !< data array

  integer, parameter :: rank = 2                   ! dimension of the array = 1
  integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr2D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_REAL_KIND)
  ptr2D = c_loc(data(1,1))
  dim2D(1) = n1
  dim2D(2) = n2
  ! Create the dataspace.
  call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real8_2D


 !=======================================================================
 !> Write a 3-D real(kind=8) array in a serial file
 subroutine hdf5_write_real8_3D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                 !< name of the dataset
  integer, intent(in) :: n(3)                                 !< first dimension of the data
  real(kind=8), dimension(n(1),n(2),n(3)), intent(in), target :: data  !< data array

  integer, parameter :: rank = 3                   ! dimension of the array = 1
  integer(hsize_t), dimension(3) :: dim3D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr3D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_REAL_KIND)
  ptr3D = c_loc(data(1,1,1))
  dim3D = n

  ! Create the dataspace.
  call h5screate_simple_f(rank, dim3D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr3D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real8_3D

  !=======================================================================
 !> Write a 4-D real(kind=8) array in a serial file
 subroutine hdf5_write_real8_4D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                 !< name of the dataset
  integer, intent(in) :: n(4)                                 !< first dimension of the data
  real(kind=8), dimension(n(1),n(2),n(3),n(4)), intent(in), target :: data  !< data array

  integer, parameter :: rank = 4                   ! dimension of the array = 1
  integer(hsize_t), dimension(4) :: dim4D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr4D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_REAL_KIND)
  ptr4D = c_loc(data(1,1,1,1))
  dim4D = n

  ! Create the dataspace.
  call h5screate_simple_f(rank, dim4D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr4D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real8_4D

  !=======================================================================
 !> Write a 5-D real(kind=8) array in a serial file
 subroutine hdf5_write_real8_5D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
  character(len=*), intent(in) :: name                 !< name of the dataset
  integer, intent(in) :: n(5)                                 !< first dimension of the data
  real(kind=8), dimension(n(1),n(2),n(3),n(4),n(5)), intent(in), target :: data  !< data array

  integer, parameter :: rank = 5                   ! dimension of the array = 1
  integer(hsize_t), dimension(5) :: dim5D          ! number of elements = size
  integer(hid_t) :: dset_id                        ! id of the dataset
  integer(hid_t) :: dspace_id                      ! id of the dataspace
  integer(hid_t) :: h5_kind                        ! HDF5 real type
  type(c_ptr) :: ptr5D                             ! pointer to the array


  h5_kind = h5kind_to_type(8,H5_REAL_KIND)
  ptr5D = c_loc(data(1,1,1,1,1))
  dim5D = n

  ! Create the dataspace.
  call h5screate_simple_f(rank, dim5D, dspace_id, h5err)

  ! Create the dataset with default properties.
  call h5dcreate_f(id, name, h5_kind, dspace_id, &
    dset_id, h5err)

  call h5dwrite_f(dset_id, h5_kind, ptr5D, h5err)

  ! End access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, h5err)

  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, h5err)

 end subroutine hdf5_write_real8_5D


 !===============================================!
 !         SERIAL DATA INPUT ROUTINES
 !===============================================!


 !=======================================================================
 !> Read a 1-D integer(kind=4) array from a serial file
 subroutine hdf5_read_int4_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n1                              !< dimension of the array to read
  integer(kind=4), dimension(n1), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

#ifdef DEBUG
  print *,'Enter hdf5_read_int4_1D for data ',name
#endif

  ptr = c_loc(data(1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)
  if(h5err < 0) then
    print *,'h5dopen_f failed for ',name
  end if

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)
  if(h5err < 0) then
    print *,'h5dread_f failed for ',name
    print *,'size of the data array ',size(data)
  end if

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

#ifdef DEBUG
  print *,'Exit hdf5_read_int4_1D for data ',name
#endif

 end subroutine hdf5_read_int4_1D


 !=======================================================================
 !> Read a 1-D integer(kind=8) scalar from a serial file
 subroutine hdf5_read_int8_0D(id, name, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=8), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data)
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_int8_0D


 !=======================================================================
 !> Read a 1-D integer(kind=8) array from a serial file
 subroutine hdf5_read_int8_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n1                              !< dimension of the array to read
  integer(kind=8), dimension(n1), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_int8_1D


 !=======================================================================
 !> Read a 2-D integer(kind=4) array from a serial file
 subroutine hdf5_read_int4_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                                   !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                              !< name of the dataset
  integer(kind=4), intent(in) :: n1                                  !< first dimension of the array to read
  integer(kind=4), intent(in) :: n2                                  !< second dimension of the array to read
  integer(kind=4), dimension(n1,n2), intent(inout), target :: data   !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_int4_2D


 !=======================================================================
 !> Read a 2-D integer(kind=4) array from a serial file
 subroutine hdf5_read_int4_3D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                                   !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                              !< name of the dataset
  integer(kind=4), intent(in) :: n(3)                                  !< first dimension of the array to read
  integer(kind=4), dimension(n(1),n(2),n(3)), intent(inout), target :: data   !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type


  ptr = c_loc(data(1,1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_int4_3D

 !=======================================================================
 !> Read a 2-D integer(kind=8) array from a serial file
 subroutine hdf5_read_int8_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                                  !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                             !< name of the dataset
  integer(kind=4), intent(in) :: n1                                 !< first dimension of the array to read
  integer(kind=4), intent(in) :: n2                                 !< second dimension of the array to read
  integer(kind=8), dimension(n1,n2), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_int8_2D


 !=======================================================================
 !> Read a 1-D real(kind=4) array from a serial file
 subroutine hdf5_read_real4_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                       !< name of the dataset
  integer(kind=4), intent(in) :: n1                           !< dimension of the array to read
  real(kind=4), dimension(n1), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(4,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real4_1D


 !=======================================================================
 !> Read a 1-D real(kind=8) array from a serial file
 subroutine hdf5_read_real8_1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                       !< name of the dataset
  integer(kind=4), intent(in) :: n1                           !< dimension of the array to read
  real(kind=8), dimension(n1), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real8_1D


 !=======================================================================
 !> Read a 2-D real(kind=4) array from a serial file
 subroutine hdf5_read_real4_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n1                              !< first dimension of the array to read
  integer(kind=4), intent(in) :: n2                              !< second dimension of the array to read
  real(kind=4), dimension(n1,n2), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(4,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real4_2D

 !=======================================================================
 !> Read a 3-D real(kind=4) array from a serial file
 subroutine hdf5_read_real4_3D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n(3)                              !< first dimension of the array to read
  real(kind=4), dimension(n(1),n(2),n(3)), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(4,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real4_3D


 !=======================================================================
 !> Read a 2-D real(kind=8) array from a serial file
 subroutine hdf5_read_real8_2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n1                              !< first dimension of the array to read
  integer(kind=4), intent(in) :: n2                              !< second dimension of the array to read
  real(kind=8), dimension(n1,n2), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real8_2D


 !=======================================================================
 !> Read a 3-D real(kind=8) array from a serial file
 subroutine hdf5_read_real8_3D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n(3)                              !< first dimension of the array to read
  real(kind=8), dimension(n(1),n(2),n(3)), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real8_3D


  !=======================================================================
 !> Read a 4-D real(kind=8) array from a serial file
 subroutine hdf5_read_real8_4D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n(4)                              !< first dimension of the array to read
  real(kind=8), dimension(n(1),n(2),n(3),n(4)), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1,1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real8_4D

  !=======================================================================
 !> Read a 5-D real(kind=8) array from a serial file
 subroutine hdf5_read_real8_5D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
  character(len=*), intent(in) :: name                          !< name of the dataset
  integer(kind=4), intent(in) :: n(5)                              !< first dimension of the array to read
  real(kind=8), dimension(n(1),n(2),n(3),n(4),n(5)), intent(inout), target :: data  !< array

  type(c_ptr) :: ptr
  integer(hid_t) :: dset_id                               ! id of the dataset
  integer(hid_t) :: h5_kind                               ! HDF5 integer type

  ptr = c_loc(data(1,1,1,1,1))
  ! hdf5 type corresponding to the integer type to read
  h5_kind = h5kind_to_type(8,H5_REAL_KIND)

  ! open the dataset
  call h5dopen_f(id, name, dset_id, h5err)

  ! read the dataset
  call h5dread_f(dset_id, h5_kind, ptr, h5err)

  ! close the dataset
  call h5dclose_f(dset_id, h5err)

 end subroutine hdf5_read_real8_5D


 !===============================================!
 !       SERIAL ATTRIBUTE OUTPUT ROUTINES
 !===============================================!


 !=======================================================================
 !> Write an integer4 as attribute
 subroutine hdf5_write_int4_attr0D(id, name, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: data                        !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = 1
  rank = 1
  ! Create scalar data space for the attribute.
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  ! Create datatype for the attribute.
  alen = 4
  call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  ! Create dataset attribute.
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  ! Write the attribute data.
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  ! Close the attribute.
  call h5aclose_f(attr_id, h5err)
  ! Terminate access to the data space.
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_int4_attr0D


 !=======================================================================
 !> Write a 1-D array of integer4 as attribute
 subroutine hdf5_write_int4_attr1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array
  integer(kind=4), dimension(n1), intent(in) :: data         !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  rank = 1
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 4
  call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_int4_attr1D


 !=======================================================================
 ! Write a 2-D array of integer4 as attribute
 subroutine hdf5_write_int4_attr2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array
  integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array
  integer(kind=4), dimension(n1,n2), intent(in) :: data      !< attribute value

  integer(hsize_t), dimension(2) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  dim(2) = n2
  rank = 2
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 4
  call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_int4_attr2D



 !=======================================================================
 ! Write a 3-D array of integer4 as attribute
 subroutine hdf5_write_int4_attr3D(id, name, n, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                !< name of the attribute
  integer(kind=4), intent(in) :: n(3)                          !< first dimension of the attribute array
  integer(kind=4), dimension(n(1),n(2),n(3)), intent(in) :: data   !< attribute value

  integer(hsize_t), dimension(3) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n(1)
  dim(2) = n(2)
  dim(3) = n(3)
  rank = 3
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 4
  call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_int4_attr3D

 !=======================================================================
 !> Write a real4 as attribute
 subroutine hdf5_write_real4_attr0D(id, name, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  real(kind=4), intent(in) :: data                           !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = 1
  rank = 1
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 4
  call h5tcopy_f(H5T_NATIVE_REAL, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_real4_attr0D


 !=======================================================================
 !> Write a 1-D array of real4 as attribute
 subroutine hdf5_write_real4_attr1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array
  real(kind=4), dimension(n1), intent(in) :: data            !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  rank = 1
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 4
  call h5tcopy_f(H5T_NATIVE_REAL, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_real4_attr1D


 !=======================================================================
 ! Write a 2-D array of real4 as attribute
 subroutine hdf5_write_real4_attr2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array
  integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array
  real(kind=4), dimension(n1,n2), intent(in) :: data         !< attribute value

  integer(hsize_t), dimension(2) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  dim(2) = n2
  rank = 2
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 4
  call h5tcopy_f(H5T_NATIVE_REAL, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_real4_attr2D


 !=======================================================================
 !> Write a real8 as attribute
 subroutine hdf5_write_real8_attr0D(id, name, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  real(kind=8), intent(in) :: data                           !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = 1
  rank = 1
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 8
  call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_real8_attr0D


 !=======================================================================
 !> Write a 1-D array of real8 as attribute
 subroutine hdf5_write_real8_attr1D(id, name, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array
  real(kind=8), dimension(n1), intent(in) :: data            !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  rank = 1
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 8
  call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_real8_attr1D


 !=======================================================================
 ! Write a 2-D array of real8 as attribute
 subroutine hdf5_write_real8_attr2D(id, name, n1, n2, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array
  integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array
  real(kind=8), dimension(n1,n2), intent(in) :: data         !< attribute value

  integer(hsize_t), dimension(2) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  dim(2) = n2
  rank = 2
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = 8
  call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_real8_attr2D


 !=======================================================================
 !> Write a string as attribute
 subroutine hdf5_write_char_attr(id, name, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                !< name of the attribute
  character(*), intent(in) :: data                           !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes

  dim(1) = 1
  alen = len(data)
  call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_char_attr


 !=======================================================================
 !> Write a 1-D array of character(LEN=H5STRLEN) as attribute
 subroutine hdf5_write_char_attr1D(id, name, len1, n1, data)

  implicit none

  integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
  character(len=*), intent(in) :: name                !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array
  integer,intent(in) :: len1                                 !< length of the array element
  character(len=len1), dimension(n1), intent(in) :: data !< attribute value

  integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
  integer :: rank
  integer(hid_t) :: aspace_id                                ! id of the dataspace
  integer(hid_t) :: atype_id                                 ! id of the dataspace
  integer(hid_t) :: attr_id                                  ! id of the dataspace
  integer(size_t) :: alen                                    ! Attribute length in bytes


  dim(1) = n1
  rank = 1
  call h5screate_simple_f(rank, dim, aspace_id, h5err)

  alen = len1
  call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
  call h5tset_size_f(atype_id, alen, h5err)
  call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
  call h5awrite_f(attr_id, atype_id, data, dim, h5err)
  call h5aclose_f(attr_id, h5err)
  call h5sclose_f(aspace_id, h5err)
  call h5tclose_f(atype_id, h5err)

 end subroutine hdf5_write_char_attr1D


 !===============================================!
 !         ATTRIBUTE INPUT ROUTINES
 !===============================================!

 !=======================================================================
 !> Read an integer4 attribute in a hdf5 file
 subroutine hdf5_read_int4_attr0D(id, name, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id              !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name              !< name of the attribute
  integer(kind=4), intent(inout) :: attr             !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D       !< dimension of the attribute: here = 1
  integer(kind=hid_t) :: attr_id                     !< attribute id
  integer(kind=hid_t) :: atype_id                    !< attribute type id

  dim1D(1) = 1

  call h5aopen_name_f(id, name, attr_id, h5err)
  if(h5err/=0) then
    print *,'Attribute ',name,'not found in the hdf5 file.'
    return
  end if
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_int4_attr0D


 !=======================================================================
 !> Read a 1-D integer4 attribute in a hdf5 file
 subroutine hdf5_read_int4_attr1D(id, name, n1, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                   !< name of the attribute
  integer(kind=4), intent(in) :: n1                       !< dimension of the attribute array
  integer(kind=4), dimension(n1), intent(inout) :: attr   !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D            !< dimension of the attribute
  integer(kind=hid_t) :: attr_id                          !< attribute id
  integer(kind=hid_t) :: atype_id                         !< attribute type id

  dim1D(1) = n1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_int4_attr1D


 !=======================================================================
 !> Read a 2-D integer4 attribute in a hdf5 file
 subroutine hdf5_read_int4_attr2D(id, name, n1, n2, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                      !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array
  integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array
  integer(kind=4), dimension(n1,n2), intent(inout) :: attr   !< attribute value

  integer(kind=hsize_t), dimension(2) :: dim2D               !< dimensions of the attribute
  integer(kind=hid_t) :: attr_id                             !< attribute id
  integer(kind=hid_t) :: atype_id                            !< attribute type id

  dim2D(1) = n1
  dim2D(2) = n2

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_int4_attr2D

 !=======================================================================
 !> Read a 3-D integer4 attribute in a hdf5 file
 subroutine hdf5_read_int4_attr3D(id, name, n, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                      !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                      !< name of the attribute
  integer(kind=4), intent(in) :: n(3)                          !< first dimension of the attribute array
  integer(kind=4), dimension(n(1),n(2),n(3)), intent(inout) :: attr   !< attribute value

  integer(kind=hsize_t), dimension(3) :: dim2D               !< dimensions of the attribute
  integer(kind=hid_t) :: attr_id                             !< attribute id
  integer(kind=hid_t) :: atype_id                            !< attribute type id

  dim2D(1) = n(1)
  dim2D(2) = n(2)
  dim2D(3) = n(3)

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_int4_attr3D


 !=======================================================================
 !> Read a real4 attribute in a hdf5 file
 subroutine hdf5_read_real4_attr0D(id, name, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id            !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name            !< name of the attribute
  real(kind=4), intent(inout) :: attr              !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D     !< dimension of the attribute: here = 1
  integer(kind=hid_t) :: attr_id                   !< attribute id
  integer(kind=hid_t) :: atype_id                  !< attribute type id

  dim1D(1) = 1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_real4_attr0D


 !=======================================================================
 !> Read a 1-D real4 attribute in a hdf5 file
 subroutine hdf5_read_real4_attr1D(id, name, n1, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                 !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                 !< name of the attribute
  integer(kind=4), intent(in) :: n1                     !< dimension of the attribute array
  real(kind=4), dimension(n1), intent(inout) :: attr    !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D          !< dimension of the attribute
  integer(kind=hid_t) :: attr_id                        !< attribute id
  integer(kind=hid_t) :: atype_id                       !< attribute type id

  dim1D(1) = n1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_real4_attr1D


 !=======================================================================
 !> Read a 2-D real4 attribute in a hdf5 file
 subroutine hdf5_read_real4_attr2D(id, name, n1, n2, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                   !< name of the attribute
  integer(kind=4), intent(in) :: n1                       !< first dimension of the attribute array
  integer(kind=4), intent(in) :: n2                       !< second dimension of the attribute array
  real(kind=4), dimension(n1,n2), intent(inout) :: attr   !< attribute value

  integer(kind=hsize_t), dimension(2) :: dim2D            !< dimensions of the attribute
  integer(kind=hid_t) :: attr_id                          !< attribute id
  integer(kind=hid_t) :: atype_id                         !< attribute type id

  dim2D(1) = n1
  dim2D(2) = n2

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_real4_attr2D


 !=======================================================================
 !> Read a real8 attribute in a hdf5 file
 subroutine hdf5_read_real8_attr0D(id, name, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id            !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name            !< name of the attribute
  real(kind=8), intent(inout) :: attr              !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D     !< dimension of the attribute: here = 1
  integer(kind=hid_t) :: attr_id                   !< attribute id
  integer(kind=hid_t) :: atype_id                  !< attribute type id

  dim1D(1) = 1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_real8_attr0D


 !=======================================================================
 !> Read a 1-D real8 attribute in a hdf5 file
 subroutine hdf5_read_real8_attr1D(id, name, n1, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                 !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                 !< name of the attribute
  integer(kind=4), intent(in) :: n1                     !< dimension of the attribute array
  real(kind=8), dimension(n1), intent(inout) :: attr    !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D          !< dimension of the attribute
  integer(kind=hid_t) :: attr_id                        !< attribute id
  integer(kind=hid_t) :: atype_id                       !< attribute type id

  dim1D(1) = n1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_real8_attr1D


 !=======================================================================
 !> Read a 2-D real4 attribute in a hdf5 file
 subroutine hdf5_read_real8_attr2D(id, name, n1, n2, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                   !< name of the attribute
  integer(kind=4), intent(in) :: n1                       !< first dimension of the attribute array
  integer(kind=4), intent(in) :: n2                       !< second dimension of the attribute array
  real(kind=8), dimension(n1,n2), intent(inout) :: attr   !< attribute value

  integer(kind=hsize_t), dimension(2) :: dim2D            !< dimensions of the attribute
  integer(kind=hid_t) :: attr_id                          !< attribute id
  integer(kind=hid_t) :: atype_id                         !< attribute type id

  dim2D(1) = n1
  dim2D(2) = n2

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_real8_attr2D


 !=======================================================================
 !> Read string attribute in a hdf5 file
 subroutine hdf5_read_char_attr(id, name, n1, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                 !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                 !< name of the attribute
  integer(kind=4), intent(in) :: n1                     !< dimension of the attribute string
  character(len=n1), intent(inout) :: attr              !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D          !< dimension of the attribute
  integer(kind=hid_t) :: attr_id                        !< attribute id
  integer(kind=hid_t) :: atype_id                       !< attribute type id

  dim1D(1) = n1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_char_attr


  !=======================================================================
 !> Read a 1-D character attribute in a hdf5 file
 subroutine hdf5_read_char_attr1D(id, name, len1, n1, attr)

  implicit none

  integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
  character(len=*), intent(in) :: name                   !< name of the attribute
  integer(kind=4), intent(in) :: n1                       !< dimension of the attribute array
  integer,intent(in) :: len1                                 !< length of the array element
  character(len=len1), dimension(n1), intent(inout) :: attr !< attribute value

  integer(kind=hsize_t), dimension(1) :: dim1D            !< dimension of the attribute
  integer(kind=hid_t) :: attr_id                          !< attribute id
  integer(kind=hid_t) :: atype_id                         !< attribute type id

  dim1D(1) = n1

  call h5aopen_name_f(id, name, attr_id, h5err)
  call h5aget_type_f(attr_id, atype_id, h5err)
  call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
  call h5aclose_f(attr_id, h5err)

 end subroutine hdf5_read_char_attr1D


end module modhdf5
