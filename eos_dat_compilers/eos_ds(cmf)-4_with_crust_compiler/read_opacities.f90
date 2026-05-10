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
!

! Example routine to read the hdf5 file and output opacity
! for a given point in T,n_B, Y_e of the EoS grid
! The information on the values of temperature, nb and Y_e are read
! from the Compose files eos.t, eos.nb and eos.yq.
! Thus make sure that they are consistent with the underlying EoS of the
! opacity data table

 subroutine read_opacity_hdf5


   use hdf5
   use modhdf5

   implicit none

   integer,parameter :: dp = selected_real_kind(15)

   character(len=500) :: file_read

   integer itmin,itmax,itmax_kappa,jbmin,jbmax,jbmax_kappa,kymin,kymax,kymax_kappa
   integer npts,ndmax,jd,k,pt,it,jb,ky
   integer, allocatable :: nd(:,:,:)

   integer :: n3D(3),n4D(4),n5D(5)
   integer(hid_t) ::  h5id, h5file_read
   real(dp),allocatable :: enumin(:,:,:,:), enumax(:,:,:,:), coeffs(:,:,:,:,:)
   real(dp),allocatable :: tab_temp(:),tab_ye(:),tab_nb(:)
   character(LEN=1),parameter :: as = '"'
   logical :: it_exist

   real(dp) :: kappa_nu,kappa_nubar,alpha,beta,enu

   ! reading T,nb,Ye from Eos data files

   open(20,file='eos.t', status='old')
   read(20,*) itmin
   read(20,*) itmax
   itmax = itmax-itmin + 1
   allocate(tab_temp(itmax))
   do it = 1,itmax
      read(20,*) tab_temp(it)
   end do
   close(20)

   open(20,file='eos.nb', status='old')
   read(20,*) jbmin
   read(20,*) jbmax
   jbmax = jbmax-jbmin + 1
   allocate(tab_nb(jbmax))
   do jb = 1,jbmax
      read(20,*) tab_nb(jb)
   end do
   close(20)

   open(20,file='eos.yq', status='old')
   read(20,*) kymin
   read(20,*) kymax
   kymax = kymax-kymin + 1
   allocate(tab_ye(kymax))
   do ky = 1,kymax
      read(20,*) tab_ye(ky)
   end do
   close(20)

   write(*,*) 'File name for reading opacity data?'
   read(*,*) file_read
   write(*,*) 'Index for temperature?'
   read(*,*) it
   write(*,*) 'Index for density?'
   read(*,*) jb
   write(*,*) 'Index for ye?'
   read(*,*) ky

   write(*,*)

   inquire(file=file_read,exist=it_exist)
   if(.not. it_exist) then
      print '(3a)', "Error : the file '",trim(file_read),"' does not exist"
      return
   endif

   call hdf5_init()

   call hdf5_open_file(file_read,h5file_read,'r')

   ! First for neutrinos

   write(*,*) '**************************************'
   write(*,*) '***           Neutrinos            ***'
   write(*,*) '**************************************'

   ! read dimensions of table arrays
   call hdf5_open_group(h5file_read,'nu',h5id)
   call hdf5_read_attr(h5id,'npts',npts)
   call hdf5_read_attr(h5id,'pts_t',itmax_kappa)
   call hdf5_read_attr(h5id,'pts_nb',jbmax_kappa)
   call hdf5_read_attr(h5id,'pts_ye',kymax_kappa)
   call hdf5_read_attr(h5id,'nd_max',ndmax)

   ! check if EoS grid indices are within boundaries and coincide with EoS
   if(itmax_kappa.ne.itmax) then
      write(*,*) 'Number of points in temperature of EoS and opacity data not equal'
      return
   end if

   if(jbmax_kappa.ne.jbmax) then
      write(*,*) 'Number of points in baryon number density of EoS and opacity data not equal'
      return
   end if

      if(kymax_kappa.ne.kymax) then
      write(*,*) 'Number of points in electron fraction of EoS and opacity data not equal'
      return
   end if



   if((it.lt.1).or.(it.gt.itmax)) then
      write(*,*) 'Temperature index', it, 'out of range: ', 1,itmax
      return
   else
      write(*,*) 'Calculating for temperature', tab_temp(it),' MeV'
   end if
   if((jb.lt.1).or.(jb.gt.jbmax)) then
      write(*,*) 'Density index', jb, 'out of range: ', 1,jbmax
      return
   else
      write(*,*) 'Calculating for density', tab_nb(jb),' fm^{-3}'
   end if
   if((ky.lt.1).or.(ky.gt.kymax)) then
      write(*,*) 'Electron fraction index', ky, 'out of range: ', 1,kymax
      return
   else
      write(*,*) 'Calculating for electron fraction', tab_ye(ky)
   end if

   allocate(nd(itmax,jbmax,kymax))
   allocate(enumin(ndmax,itmax,jbmax,kymax))
   allocate(enumax(ndmax,itmax,jbmax,kymax))
   allocate(coeffs(npts,ndmax,itmax,jbmax,kymax))

   n3d(1) = itmax
   n3d(2) = jbmax
   n3d(3) = kymax
   n4d(2:4) = n3d(1:3)
   n4d(1) = ndmax
   n5d(2:5) = n4d(1:4)
   n5d(1) = npts

   call hdf5_read_data(h5id,'enumax',n4d,enumax)
   call hdf5_read_data(h5id,'enumin',n4d,enumin)
   call hdf5_read_data(h5id,'nd_tny',n3d,nd)
   call hdf5_read_data(h5id,'coeffs',n5d,coeffs)


   call hdf5_close_group(h5id)


   kappa_nu = 0._dp
   write(*,*) 'Which energy (in MeV) between ',&
        real(exp(enumin(1,it,jb,ky))),'and',&
        real(exp(enumax(nd(it,jb,ky),it,jb,ky))),'? '
   read(*,*) enu
   pt = 0
   domainloop_nu: DO jd = 1,nd(it,jb,ky)
      IF(enu.le.exp(enumax(jd,it,jb,ky)).and.enu.ge.exp(enumin(jd,it,jb,ky))) then
         pt = jd
         exit domainloop_nu
      end if
   end DO domainloop_nu
   IF(pt.ne.0) then
      alpha = (enumax(pt,it,jb,ky)-enumin(pt,it,jb,ky))/2
      beta = (enumax(pt,it,jb,ky)+enumin(pt,it,jb,ky))/2

      kappa_nu = coeffs(1,pt,it,jb,ky)
      DO k = 2,npts
         kappa_nu = kappa_nu + coeffs(k,pt,it,jb,ky)*((log(enu)-beta)/alpha)**(k-1)
      end DO
   else
      write(*,*) 'Energy', enu,'out of range'
      return
   end if
   write(*,*) '**************************************'

   write(*,*) 'Opacity for nu_e is',exp(kappa_nu),' 1/cm at energy', enu, ' MeV'

   if(allocated(nd)) deallocate(nd)
   if(allocated(enumin)) deallocate(enumin)
   if(allocated(enumax)) deallocate(enumax)
   if(allocated(coeffs)) deallocate(coeffs)

   ! Then for anti-neutrinos
   write(*,*) '**************************************'
   write(*,*) '***         Anti-neutrinos         ***'
   write(*,*) '**************************************'

   ! read dimensions of table arrays
   call hdf5_open_group(h5file_read,'nu_bar',h5id)
   call hdf5_read_attr(h5id,'npts',npts)
   call hdf5_read_attr(h5id,'pts_t',itmax_kappa)
   call hdf5_read_attr(h5id,'pts_nb',jbmax_kappa)
   call hdf5_read_attr(h5id,'pts_ye',kymax_kappa)
   call hdf5_read_attr(h5id,'nd_max',ndmax)


   ! check if EoS grid indices are within boundaries and coincide with EoS
   if(itmax_kappa.ne.itmax) then
      write(*,*) 'Number of points in temperature of EoS and opacity data not equal'
      return
   end if

   if(jbmax_kappa.ne.jbmax) then
      write(*,*) 'Number of points in baryon number density of EoS and opacity data not equal'
      return
   end if

      if(kymax_kappa.ne.kymax) then
      write(*,*) 'Number of points in electron fraction of EoS and opacity data not equal'
      return
   end if


   if((it.lt.1).or.(it.gt.itmax)) then
      write(*,*) 'Temperature index', it, 'out of range: ', 1,itmax
      return
   else
      write(*,*) 'Calculating for temperature     ', real(tab_temp(it)),' MeV'
   end if
   if((jb.lt.1).or.(jb.gt.jbmax)) then
      write(*,*) 'Density index', jb, 'out of range: ', 1,jbmax
      return
   else
      write(*,*) 'Calculating for density         ', real(tab_nb(jb)),' fm^{-3}'
   end if
   if((ky.lt.1).or.(ky.gt.kymax)) then
      write(*,*) 'Electron fraction index', ky, 'out of range: ', 1,kymax
      return
   else
      write(*,*) 'Calculating for electron fraction', real(tab_ye(ky))
   end if

   if(allocated(tab_temp)) deallocate(tab_temp)
   if(allocated(tab_nb)) deallocate(tab_nb)
   if(allocated(tab_ye)) deallocate(tab_ye)


   allocate(nd(itmax,jbmax,kymax))
   allocate(enumin(ndmax,itmax,jbmax,kymax))
   allocate(enumax(ndmax,itmax,jbmax,kymax))
   allocate(coeffs(npts,ndmax,itmax,jbmax,kymax))

   n3d(1) = itmax
   n3d(2) = jbmax
   n3d(3) = kymax
   n4d(2:4) = n3d(1:3)
   n4d(1) = ndmax
   n5d(2:5) = n4d(1:4)
   n5d(1) = npts


   call hdf5_read_data(h5id,'enumax',n4d,enumax)
   call hdf5_read_data(h5id,'enumin',n4d,enumin)
   call hdf5_read_data(h5id,'nd_tny',n3d,nd)
   call hdf5_read_data(h5id,'coeffs',n5d,coeffs)


   call hdf5_close_group(h5id)


   kappa_nubar = 0._dp
   write(*,*) 'Which energy (in MeV) between',&
        real(exp(enumin(1,it,jb,ky))),'and'&
        ,real(exp(enumax(nd(it,jb,ky),it,jb,ky))),')? '
   read(*,*) enu
   write(*,*) '********************************'
   pt = 0
   domainloop_nubar: DO jd = 1,nd(it,jb,ky)
      IF(enu.le.exp(enumax(jd,it,jb,ky)).and.enu.ge.exp(enumin(jd,it,jb,ky))) then
         pt = jd
         exit domainloop_nubar
      end if
   end DO domainloop_nubar
   IF(pt.ne.0) then
      alpha = (enumax(pt,it,jb,ky)-enumin(pt,it,jb,ky))/2
      beta = (enumax(pt,it,jb,ky)+enumin(pt,it,jb,ky))/2

      kappa_nubar = coeffs(1,pt,it,jb,ky)
      DO k = 2,npts
         kappa_nubar = kappa_nubar + coeffs(k,pt,it,jb,ky)*((log(enu)-beta)/alpha)**(k-1)
      end DO
   else
      write(*,*) 'Energy', enu,'out of range'
      return
   end if
   write(*,*) 'Opacity for nubar_e is',exp(kappa_nubar),' 1/cm at energy', enu, ' MeV'

   if(allocated(nd)) deallocate(nd)
   if(allocated(enumin)) deallocate(enumin)
   if(allocated(enumax)) deallocate(enumax)
   if(allocated(coeffs)) deallocate(coeffs)

   call hdf5_close_file(h5file_read)

 end subroutine read_opacity_hdf5


 program read_opacity_data

   write(*,*) 'Sample program to read opacity data from the table in HDF5 format '
   write(*,*) 'as downloaded from the Compose data base. Please see arxiv:2020.xxx'
   write(*,*) 'for the definition of the different quantities. '
   write(*,*) 'Attention: The information on the values of temperature, nb and Y_e is read'
   write(*,*) 'from the Compose files eos.t, eos.nb and eos.yq.'
   write(*,*) 'Thus make sure that they are present and consistent with the'
   write(*,*) 'underlying EoS of the opacity data table'
   call read_opacity_hdf5
 end program read_opacity_data
