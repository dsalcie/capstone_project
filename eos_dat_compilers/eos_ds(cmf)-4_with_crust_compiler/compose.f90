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


!***********************************************************************
!***********************************************************************
! Fortran90 program from the CompOSE website http://compose.obspm.fr
!    for reading, checking, interpolating and writing EoS tables.
! Required standard input files
!    from website: eos.t, eos.nb, eos.yq, eos.thermo;
!    supplied by user in the standard version
!    (not in the terminal version):
!    eos.parameters, eos_quantities
! Optional additional input files
!    from website: eos.b, eos.compo, eos.micro.
! See manual for details and alternative options.
! CompOSE core team, version 2.16, 2017/12/13
!***********************************************************************
!***********************************************************************
PROGRAM compose
!***********************************************************************
!***********************************************************************
! Stefan Typel for the CompOSE core team, version 1.09, 2017/09/28
! working with CompOSE EoS tables

 implicit none

 integer :: iwr,iterm

! 1: terminal version for input,
! else: standard version with input from files
iterm = 1

! 1: write information and report to terminal,
! else: not
iwr = 1

if (iterm == 1) then
   !+++++++++++++++++++++
   call run_terminal(iwr)
   !+++++++++++++++++++++
else

   ! read and analyze eos data tables
   ! to be called only once
   !+++++++++++++++++++++++
   call init_eos_table(iwr)
   !+++++++++++++++++++++++

   ! define quantities in eos to be interpolated
   !+++++++++++++++++++++++++
   call define_eos_table(iwr)
   !+++++++++++++++++++++++++

   ! generate eos table
   !++++++++++++++++++++++
   call get_eos_table(iwr)
   !++++++++++++++++++++++

end if

! output quantities
! eos_thermo(1): p [MeV fm^-3]
! eos_thermo(2): S [dimensionless]
! eos_thermo(3): mu_b - m_n [MeV]
! eos_thermo(4): mu_q [MeV]
! eos_thermo(5): mu_l [MeV]
! eos_thermo(6): F/m_n-1 [dimensionless]
! eos_thermo(7): E/m_n-1 [dimensionless]
! eos_thermo(8): H/m_n-1 [dimensionless]
! eos_thermo(9): G/m_n-1 [dimensionless]
! eos_thermo(10): dp/dnb|E [MeV]
! eos_thermo(11): dp/dE|nb [fm^-3]
! eos_thermo(12): c_s^2 [dimensionless]
! eos_thermo(13): c_V [dimensionless]
! eos_thermo(14): c_p [dimensionless]
! eos_thermo(15): Gamma [dimensionless]
! eos_thermo(16): alpha_p [MeV^-1]
! eos_thermo(17): beta_V [fm^-3]
! eos_thermo(18): kappa_T [MeV^-1 fm^3]
! eos_thermo(19): kappa_S [MeV^-1 fm^3]
! eos_thermo(20): F [MeV]
! eos_thermo(21): E [MeV]
! eos_thermo(22): H [MeV]
! eos_thermo(23): G [MeV]
! eos_thermo(24): epsilon (energy density) [MeV/fm^3]
! eos_thermo_add(1), ..., eos_thermo(n_add), at least two
! eos_compo_p(1), ..., eos_compo_p(n_p)
! eos_compo_q(1,1), ..., eos_compo_q(n_q,1): Y^av
! eos_compo_q(1,2), ..., eos_compo_q(n_q,2): A^av
! eos_compo_q(1,3), ..., eos_compo_q(n_q,3): Z^av
! eos_compo_q(1,4), ..., eos_compo_q(n_q,4): N^av
! eos_micro(1), ..., eos_micro(n_m)
! eos_thermo_err(1): absolute error in f/n [MeV]
! eos_thermo_err(2): relative error in f/n [dimensionless]
! eos_thermo_err(3): absolute error in e/n [MeV]
! eos_thermo_err(4): relative error in e/n [dimensionless]
! eos_thermo_err(5): absolute error in p/n [MeV]
! eos_thermo_err(6): relative error in p/n [dimensionless]
! eos_thermo_err(7): absolute error in s/n [dimensionless]
! eos_thermo_err(8): relative error in s/n [dimensionless]
!2017/05/22
! eos_df(1): F [MeV]
! eos_df(2): dF/dT []
! eos_df(3): d^2F/(dT^2) [MeV^-1]
! eos_df(4): d^2F/(dT dn_b) [fm^3]
! eos_df(5): d^2F/(dT dY_q) []
! eos_df(6): dF/dn_b [MeV fm^3]
! eos_df(7): d^2F/(dn_b^2) [MeV fm^6]
! eos_df(8): d^2F/(dn_b dY_q) [MeV fm^3]
! eos_df(9): dF/dY_q [MeV]
! eos_df(10): d^2dF/(dY_q^2) [MeV]

!2017/09/28
!stop
!***********************************************************************
end PROGRAM compose
!***********************************************************************
!***********************************************************************
! subroutines
!***********************************************************************
!***********************************************************************
subroutine init_eos_table(iwr)
! Stefan Typel for the CompOSE core team, version 1.06, 2017/11/16
 use omp_lib
 use m_get_tables
 implicit none
 integer :: nbl,iwr,iyq,iinit,init_flg = 0

! maximum number of subtables
 nbl = 10


 call read_eos_4_tables(iwr,nbl,iyq,ii_thermo=0,ii_tynb=0,unit=0,iinit=init_flg)

 call get_diff_rules()

 call init_ipl_rule()

 call get_eos_report(iwr)

end SUBROUTINE init_eos_table



!***********************************************************************
SUBROUTINE write_errors(ierr)
! Stefan Typel for the CompOSE core team, version 1.08, 2017/05/22
USE compose_internal  
implicit none
integer :: ierr,i,imax

error_msg(0) = ierr
if (ierr > 0) then
   write(*,*)
   write(*,*) '!!! warning !!!'
   write(*,*) ierr,'error(s) detected'
   imax = ierr
   if (imax > dim_err) imax = dim_err
   do i=1,imax,1
      write(*,*) 'error(',i,') =',error_msg(i)

      choose_error: select case (error_msg(i))
      case (1)
         write(*,*) 'input file eos.t does not exist'
      case (2)
         write(*,*) 'temperature in eos.t not increasing'
      case (3)
         write(*,*) 'maximum index smaller than minimum index in eos.t'
      case (5)
         write(*,*) 'index for parameter replacement in eos.b out of range'
      case (6)
         write(*,*) 'input file eos.nb does not exist'
      case (7)
         write(*,*) 'baryon number density in eos.nb not increasing'
      case (8)
         write(*,*) 'maximum index smaller than minimum index in eos.nb'
      case (11)
         write(*,*) 'input file eos.yq does not exist'
      case (12)
         write(*,*) 'baryon charge fraction in eos.yq not increasing'
      case (13)
         write(*,*) 'maximum index smaller than minimum index in eos.yq'
      case (15)
         write(*,*) 'more than three independent variables'
      case (16)
         write(*,*) 'input file eos.b does not exist'
      case (17)
         write(*,*) 'magnetic field strength in eos.b not increasing'
      case (18)
         write(*,*) 'maximum index smaller than minimum index in eos.b'
      case (19)
         write(*,*) 'at leat one of the four parameters should be constant'
      case (20)
         write(*,*) 'quantity index out of range in eos.thermo'
      case (21)
         write(*,*) 'error in allocation of memory for tab_thermo'
      case (27)
         write(*,*) 'error in reading masses from file eos.thermo'
      case (28)
         write(*,*) 'error in allocation of memory for tab_para'
      case (32)
         write(*,*) 'error in allocation of memory for r1d'
      case (33)
         write(*,*) 'error in allocation of memory for r2d'
      case (34)
         write(*,*) 'error in allocation of memory for idx_thermo2'
      case (35)
         write(*,*) 'error in allocation of memory for v_thermo'
      case (36)
         write(*,*) 'error in allocation of memory for eos_thermo_add'
      case (37)
         write(*,*) 'error in allocation of memory for idx_add'
      case (38)
         write(*,*) 'error in allocation of memory for qty in read_eos_table_thermo'
      case (41)
         write(*,*) 'temperature out of range'
      case (42)
         write(*,*) 'baryon number density out of range'
      case (43)
         write(*,*) 'baryonic charge fraction out of range'
      case (51)
         write(*,*) 'not all points of interpolation line exist'
      case (52)
         write(*,*) 'not all points of interpolation square exist'
      case (53)
         write(*,*) 'not all points of interpolation cube exist'
      case (60)
         write(*,*) 'input file eos.quantities does not exist'
      case (61)
         write(*,*) 'number of quantities out of range in eos.quantities'
      case (62)
         write(*,*) 'number of additional quantities out of range in eos.quantities'
      case (63)
         write(*,*) 'index of quantities out of range in eos.quantities'
      case (64)
         write(*,*) 'index of additional quantities out of range in eos.quantities'
      case (65)
         write(*,*) 'number of pairs out of range in eos.quantities'
      case (66)
         write(*,*) 'number of quadruples out of range in eos.quantities'
      case (67)
         write(*,*) 'number of microscopic quantities out of range in eos.quantities'
      case (68)
         write(*,*) 'number of error quantities out of range in eos.quantities'
      case (69)
         write(*,*) 'index for output format out of range'
      case (70)
         write(*,*) 'index for interpolation in T out of range'
      case (71)
         write(*,*) 'index for interpolation in n_b out of range'
      case (72)
         write(*,*) 'index for interpolation in Y_q out of range'
      case (73)
         write(*,*) 'index for interpolation in B out of range'
      case (74)
         write(*,*) 'index for independent interpolation out of range'
      case (75)
         write(*,*) 'number of derivative quantities out of range in eos.quantities'
      case (80)
         write(*,*) 'number of table entries too small'
      case (81)
         write(*,*) 'minimum of parameter T < 0'
      case (82)
         write(*,*) 'minimum of parameter T <= 0'
      case (83)
         write(*,*) 'maximum of parameter T < minimum of parameter T'
      case (84)
         write(*,*) 'minimum of parameter n_b <= 0'
      case (85)
         write(*,*) 'maximum of parameter T < minimum of parameter n_b'
      case (86)
         write(*,*) 'minimum of parameter Y_q < 0'
      case (87)
         write(*,*) 'maximum of parameter Y_q < minimum of parameter Y_q'
      case (88)
         write(*,*) 'maximum of parameter Y_q > 1'
      case (89)
         write(*,*) 'minimum of parameter B < 0'
      case (90)
         write(*,*) 'maximum of parameter B < minimum of parameter B'
      case (91)
         write(*,*) 'option for beta equilibrium not available'
      case (101)
         write(*,*) 'error in allocation of memory for iqtyp in read_eos_table_compo'
      case (102)
         write(*,*) 'error in allocation of memory for iqtyq in read_eos_table_compo'
      case (103)
         write(*,*) 'error in allocation of memory for qtyp in read_eos_table_compo'
      case (104)
         write(*,*) 'error in allocation of memory for qtyq in read_eos_table_compo'
      case (105)
         write(*,*) 'error in allocation of memory for idxp_compo in read_eos_table_compo'
      case (106)
         write(*,*) 'error in allocation of memory for idxq_compo in read_eos_table_compo'
      case (107)
         write(*,*) 'error in allocation of memory for tabp_compo in read_eos_table_compo'
      case (108)
         write(*,*) 'error in allocation of memory for tabq_compo in read_eos_table_compo'
      case (109)
         write(*,*) 'error in allocation of memory for iqtym in read_eos_table_micro'
      case (110)
         write(*,*) 'error in allocation of memory for qtym in read_eos_table_micro'
      case (111)
         write(*,*) 'error in allocation of memory for idx_mic in read_eos_table_micro'
      case (112)
         write(*,*) 'error in allocation of memory for tab_mic in read_eos_table_micro'
      case (113)
         write(*,*) 'error in allocation of memory for idx_arg'
      case (114)
         write(*,*) 'error in allocation of memory for idx_argx'
      case (115)
         write(*,*) 'error in allocation of memory for iarg'
      case (120)
         write(*,*) 'number of entries in eos.t not consistent with indices'
      case (121)
         write(*,*) 'number of entries in eos.nb not consistent with indices'
      case (122)
         write(*,*) 'number of entries in eos.yq not consistent with indices'
      case (123)
         write(*,*) 'number of entries in eos.b not consistent with indices'
      case (130)
         write(*,*) 'error in allocation of memory for iphase in get_eos_report'
      case (131)
         write(*,*) 'error in allocation of memory for idx_p in read_eos_table_compo'
      case (132)
         write(*,*) 'error in allocation of memory for idx_q in read_eos_table_compo'
      case (133)
         write(*,*) 'error in allocation of memory for idx_m in read_eos_table_compo'
      case (134)
         write(*,*) 'error in allocation of memory for idx_compo_p in read_eos_table_compo'
      case (135)
         write(*,*) 'error in allocation of memory for idx_compo_q in read_eos_table_compo'
      case (136)
         write(*,*) 'error in allocation of memory for idx_micro in read_eos_table_compo'
      case (137)
         write(*,*) 'error in allocation of memory for eos_compo_p in read_eos_table_compo'
      case (138)
         write(*,*) 'error in allocation of memory for eos_compo_q in read_eos_table_compo'
      case (139)
         write(*,*) 'error in allocation of memory for eos_micro in read_eos_table_compo'
      case (140)
         write(*,*) 'error in allocation of memory for mat in eos_interpol_d3'
      case (141)
         write(*,*) 'error in allocation of memory for mat3 in eos_interpol_d3'
      case (142)
         write(*,*) 'error in allocation of memory for vp_compo in eos_interpol_d3'
      case (143)
         write(*,*) 'error in allocation of memory for vq_compo in eos_interpol_d3'
      case (144)
         write(*,*) 'error in allocation of memory for v_micro in eos_interpol_d3'
      case (145)
         write(*,*) 'error in allocation of memory for eos_q in readeos_table_compo'
      case (146)
         write(*,*) 'error in allocation of memory for mat2 in eos_interpol_d3'
      case (151)
         write(*,*) 'error in allocation of memory for iqtyp in get_eos_report'
      case (152)
         write(*,*) 'error in allocation of memory for iqtyq in get_eos_report'
      case (153)
         write(*,*) 'error in allocation of memory for iqtym in get_eos_report'
      case (154)
         write(*,*) 'np_miss > 0 in get_eos_report'
      case (155)
         write(*,*) 'nq_miss > 0 in get_eos_report'
      case (156)
         write(*,*) 'nm_miss > 0 in get_eos_report'
      case (157)
         write(*,*) 'nphase_miss > 0 in get_eos_report'
      case (160)
         write(*,*) 'determination of beta-equilibrium not possible'
      case (161)
         write(*,*) 'determination of composition quantities not possible without file eos.compo'
      case (162)
         write(*,*) 'determination of microscopic quantities not possible without file eos.micro'
      case (163)
         write(*,*) 'error in indices in eos_interpol_d2'
      case (164)
         write(*,*) 'error in indices in eos_interpol_d1'
      case (165)
         write(*,*) 'error in first and second indices in eos_interpol_d1'
      case (170)
         write(*,*) 'error in allocation of memory for mat in eos_interpol_d2'
      case (171)
         write(*,*) 'error in allocation of memory for mat3 in eos_interpol_d2'
      case (172)
         write(*,*) 'error in allocation of memory for vp_compo in eos_interpol_d2'
      case (173)
         write(*,*) 'error in allocation of memory for vq_compo in eos_interpol_d2'
      case (174)
         write(*,*) 'error in allocation of memory for v_micro in eos_interpol_d2'
      case (175)
         write(*,*) 'number of pairs larger than dim_ip'
      case (176)
         write(*,*) 'number of quadruples larger than dim_iq'
      case (177)
         write(*,*) 'number of microscopic quantities larger than dim_im'
      end select choose_error

   end do
   write(*,*) 'program terminated'
   write(*,*)
   stop 2
end if

return
end SUBROUTINE write_errors
!***********************************************************************
subroutine get_diff_rules()
 ! Stefan Typel for the CompOSE core team, version 1.07, 2016/11/08
 USE compose_internal
 implicit none
 integer :: ip,in,it,iy,ia,ia_min,ia_max,ir,iap,iq,alloc_status,ierr
 integer, allocatable :: iarg(:)
 double precision :: z0,zm,zm2,zp,zp2

 ! choose maximum differentiation rule

 ! used grid points  -4 -3 -2 -1  0  1  2  3  4
 ! case        1            *  *  *  *  *
 ! case        2         *  *  *  *  *
 ! case       -2               *  *  *  *  *
 ! case        3      *  *  *  *  *
 ! case       -3                  *  *  *  *  *
 ! case        4            *  *  *  *
 ! case       -4               *  *  *  *
 ! case        5         *  *  *  *
 ! case       -5                  *  *  *  *
 ! case        6               *  *  *
 ! case        7            *  *  *
 ! case       -7                  *  *  *
 ! case        8               *  *
 ! case       -8                  *  *


 ierr = 0
 allocate(iarg(1:dim_a),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 115
 end if

 allocate(idx_argx(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:3),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 114
 end if

 allocate(r1d(3,1:dim_a,-4:4,-8:8),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 32
 end if

 allocate(r2d(3,1:dim_a,-4:4,-8:8),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 33
 end if

call write_errors(ierr)


 idx_argx(:,:,:,1:3) = idx_arg(:,:,:,1:3)

 do iy=1,dim_idx(3),1
   do in=1,dim_idx(2),1
     do it=1,dim_idx(1),1
       if (all(idx_argx(it,in,iy,:) == idx_ex(1:3)))  idx_arg(it,in,iy,0) = 1
     end do
   end do
 end do

 if (all(dim_idx(1:3) > 0)) idx_arg(1:dim_idx(1),1:dim_idx(2),1:dim_idx(3),1:3) = 0

 ! first index
 ip = 1

 do in=1,dim_idx(2),1
   do iy=1,dim_idx(3),1
     do it=1,(dim_idx(1)-1),1
       if ((idx_arg(it,in,iy,0) == 1).and.&
         (idx_arg(it+1,in,iy,0) == 1)) then
         iarg(it) = 1
       else
         iarg(it) = 0
       end if
     end do
     iarg(dim_idx(1)) = 0
     do it=1,(dim_idx(1)-1),1
       if ((iarg(it) == 1).and.(iarg(it+1) == 1)) then
         iarg(it) = 2
       end if
     end do
     do it=1,(dim_idx(1)-2),1
       if ((iarg(it) == 2).and.(iarg(it+2) >= 1)) then
         iarg(it) = 3
       end if
     end do
     do it=1,(dim_idx(1)-3),1
       if ((iarg(it) == 3).and.(iarg(it+3) >= 1)) then
         iarg(it) = 4
       end if
     end do

     do it=1,dim_idx(1),1
       select case(iarg(it))
       case(4)
         idx_arg(it+4,in,iy,ip) = 3
         idx_arg(it+3,in,iy,ip) = 2
         idx_arg(it+2,in,iy,ip) = 1
         if (idx_arg(it+1,in,iy,ip) == 0) idx_arg(it+1,in,iy,ip) = -2
         if (idx_arg(it  ,in,iy,ip) == 0) idx_arg(it  ,in,iy,ip) = -3
       case(3)
         if (idx_arg(it+3,in,iy,ip) == 0) idx_arg(it+3,in,iy,ip) = 5
         if (idx_arg(it+2,in,iy,ip) == 0) idx_arg(it+2,in,iy,ip) = 4
         if (idx_arg(it+1,in,iy,ip) == 0) idx_arg(it+1,in,iy,ip) = -4
         if (idx_arg(it  ,in,iy,ip) == 0) idx_arg(it  ,in,iy,ip) = -5
       case(2)
         if (idx_arg(it+2,in,iy,ip) == 0) idx_arg(it+2,in,iy,ip) = 7
         if (idx_arg(it+1,in,iy,ip) == 0) idx_arg(it+1,in,iy,ip) = 6
         if (idx_arg(it  ,in,iy,ip) == 0) idx_arg(it  ,in,iy,ip) = -7
       case(1)
         if (idx_arg(it+1,in,iy,ip) == 0) idx_arg(it+1,in,iy,ip) = 8
         if (idx_arg(it  ,in,iy,ip) == 0) idx_arg(it  ,in,iy,ip) = -8
       end select
     end do
   end do
 end do

 ! second index
 ip = 2
 do it=1,dim_idx(1),1
   do iy=1,dim_idx(3),1
     do in=1,(dim_idx(2)-1),1
       if ((idx_arg(it,in,iy,0) == 1).and.&
         (idx_arg(it,in+1,iy,0) == 1)) then
         iarg(in) = 1
       else
         iarg(in) = 0
       end if
     end do
     iarg(dim_idx(2)) = 0
     do in=1,(dim_idx(2)-1),1
       if ((iarg(in) == 1).and.(iarg(in+1) == 1)) then
         iarg(in) = 2
       end if
     end do
     do in=1,(dim_idx(2)-2),1
       if ((iarg(in) == 2).and.(iarg(in+2) >= 1)) then
         iarg(in) = 3
       end if
     end do
     do in=1,(dim_idx(2)-3),1
       if ((iarg(in) == 3).and.(iarg(in+3) >= 1)) then
         iarg(in) = 4
       end if
     end do
     do in=1,dim_idx(2),1
       select case(iarg(in))
       case(4)
         idx_arg(it,in+4,iy,ip) = 3
         idx_arg(it,in+3,iy,ip) = 2
         idx_arg(it,in+2,iy,ip) = 1
         if (idx_arg(it,in+1,iy,ip) == 0)    idx_arg(it,in+1,iy,ip) = -2
         if (idx_arg(it,in  ,iy,ip) == 0)    idx_arg(it,in  ,iy,ip) = -3
       case(3)
         if (idx_arg(it,in+3,iy,ip) == 0)    idx_arg(it,in+3,iy,ip) = 5
         if (idx_arg(it,in+2,iy,ip) == 0)    idx_arg(it,in+2,iy,ip) = 4
         if (idx_arg(it,in+1,iy,ip) == 0)    idx_arg(it,in+1,iy,ip) = -4
         if (idx_arg(it,in  ,iy,ip) == 0)    idx_arg(it,in  ,iy,ip) = -5
       case(2)
         if (idx_arg(it,in+2,iy,ip) == 0)    idx_arg(it,in+2,iy,ip) = 7
         if (idx_arg(it,in+1,iy,ip) == 0)    idx_arg(it,in+1,iy,ip) = 6
         if (idx_arg(it,in  ,iy,ip) == 0)    idx_arg(it,in  ,iy,ip) = -7
       case(1)
         if (idx_arg(it,in+1,iy,ip) == 0)    idx_arg(it,in+1,iy,ip) = 8
         if (idx_arg(it,in  ,iy,ip) == 0)    idx_arg(it,in  ,iy,ip) = -8
       end select
     end do
   end do
 end do

 ! third index
 ip = 3
 do it=1,dim_idx(1),1
   do in=1,dim_idx(2),1
     do iy=1,(dim_idx(3)-1),1
       if ((idx_arg(it,in,iy,0) == 1).and.&
         (idx_arg(it,in,iy+1,0) == 1)) then
         iarg(iy) = 1
       else
         iarg(iy) = 0
       end if
     end do
     iarg(dim_idx(3)) = 0
     do iy=1,(dim_idx(3)-1),1
       if ((iarg(iy) == 1).and.(iarg(iy+1) == 1)) then
         iarg(iy) = 2
       end if
     end do
     do iy=1,(dim_idx(3)-2),1
       if ((iarg(iy) == 2).and.(iarg(iy+2) >= 1)) then
         iarg(iy) = 3
       end if
     end do
     do iy=1,(dim_idx(3)-3),1
       if ((iarg(iy) == 3).and.(iarg(iy+3) >= 1)) then
         iarg(iy) = 4
       end if
     end do
     do iy=1,dim_idx(3),1
       select case(iarg(iy))
       case(4)
         idx_arg(it,in,iy+4,ip) = 3
         idx_arg(it,in,iy+3,ip) = 2
         idx_arg(it,in,iy+2,ip) = 1
         if (idx_arg(it,in,iy+1,ip) == 0)   idx_arg(it,in,iy+1,ip) = -2
         if (idx_arg(it,in,iy  ,ip) == 0)   idx_arg(it,in,iy  ,ip) = -3
       case(3)
         if (idx_arg(it,in,iy+3,ip) == 0)   idx_arg(it,in,iy+3,ip) = 5
         if (idx_arg(it,in,iy+2,ip) == 0)   idx_arg(it,in,iy+2,ip) = 4
         if (idx_arg(it,in,iy+1,ip) == 0)   idx_arg(it,in,iy+1,ip) = -4
         if (idx_arg(it,in,iy  ,ip) == 0)   idx_arg(it,in,iy  ,ip) = -5
       case(2)
         if (idx_arg(it,in,iy+2,ip) == 0)   idx_arg(it,in,iy+2,ip) = 7
         if (idx_arg(it,in,iy+1,ip) == 0)   idx_arg(it,in,iy+1,ip) = 6
         if (idx_arg(it,in,iy  ,ip) == 0)   idx_arg(it,in,iy  ,ip) = -7
       case(1)
         if (idx_arg(it,in,iy+1,ip) == 0)   idx_arg(it,in,iy+1,ip) = 8
         if (idx_arg(it,in,iy  ,ip) == 0)   idx_arg(it,in,iy  ,ip) = -8
       end select
     end do
   end do
 end do

 ! define differentiation rules
 do ip=1,3,1
   ia_min = 1
   ia_max = dim_idx(ip)
   ! initialisation with zero
   do ia=ia_min,ia_max,1
     r1d(ip,ia,-4:4,-8:8) = 0.d00
     r2d(ip,ia,-4:4,-8:8) = 0.d00

     ! five-point formulas
     do iq=-2,2,1
       choose_rule5: select case (iq)
       case (2)
         ir = 3
       case (1)
         ir = 2
       case (0)
         ir = 1
       case (-1)
         ir = -2
       case (-2)
         ir = -3
       end select choose_rule5
       iap = ia-iq
       if ((iap >= (ia_min+2)).and.(iap <= (ia_max-2))) then
         zp2 = tab_para(iap+2,ip)-tab_para(ia,ip)
         zp  = tab_para(iap+1,ip)-tab_para(ia,ip)
         z0  = tab_para(iap  ,ip)-tab_para(ia,ip)
         zm  = tab_para(iap-1,ip)-tab_para(ia,ip)
         zm2 = tab_para(iap-2,ip)-tab_para(ia,ip)
         ! first derivative
         r1d(ip,ia, 2-iq,ir) = &
           -((zp+z0)*zm*zm2+zp*z0*(zm+zm2))/&
           ((zp2-zp)*(zp2-z0)*(zp2-zm)*(zp2-zm2))
         r1d(ip,ia, 1-iq,ir) = &
           -((zp2+z0)*zm*zm2+zp2*z0*(zm+zm2))/&
           ((zp-zp2)*(zp-z0)*(zp-zm)*(zp-zm2))
         r1d(ip,ia,  -iq,ir) = &
           -((zp2+zp)*zm*zm2+zp2*zp*(zm+zm2))/&
           ((z0-zp2)*(z0-zp)*(z0-zm)*(z0-zm2))
         r1d(ip,ia,-1-iq,ir) = &
           -((zp2+zp)*z0*zm2+zp2*zp*(z0+zm2))/&
           ((zm-zp2)*(zm-zp)*(zm-z0)*(zm-zm2))
         r1d(ip,ia,-2-iq,ir) = &
           -((zp2+zp)*z0*zm+zp2*zp*(z0+zm))/&
           ((zm2-zp2)*(zm2-zp)*(zm2-z0)*(zm2-zm))
         ! second derivative
         r2d(ip,ia, 2-iq,ir) = 2.d00*&
           (zp*z0+(zp+z0)*(zm+zm2)+zm*zm2)/&
           ((zp2-zp)*(zp2-z0)*(zp2-zm)*(zp2-zm2))
         r2d(ip,ia, 1-iq,ir) = 2.d00*&
           (zp2*z0+(zp2+z0)*(zm+zm2)+zm*zm2)/&
           ((zp-zp2)*(zp-z0)*(zp-zm)*(zp-zm2))
         r2d(ip,ia,  -iq,ir) = 2.d00*&
           (zp2*zp+(zp2+zp)*(zm+zm2)+zm*zm2)/&
           ((z0-zp2)*(z0-zp)*(z0-zm)*(z0-zm2))
         r2d(ip,ia,-1-iq,ir) = 2.d00*&
           (zp2*zp+(zp2+zp)*(z0+zm2)+z0*zm2)/&
           ((zm-zp2)*(zm-zp)*(zm-z0)*(zm-zm2))
         r2d(ip,ia,-2-iq,ir) = 2.d00*&
           (zp2*zp+(zp2+zp)*(z0+zm)+z0*zm)/&
           ((zm2-zp2)*(zm2-zp)*(zm2-z0)*(zm2-zm))
       end if
     end do

     ! four-point formulas
     do iq=-2,2,1
       choose_rule4: select case (iq)
       case (2)
         ir = 5
       case (1)
         ir = 4
       case (-1)
         ir = -4
       case (-2)
         ir = -5
       end select choose_rule4
       iap = ia-iq
       if (iq > 0) then
         if ((iap >= (ia_min+1)).and.(iap <= (ia_max-2))) then
           zp2 = tab_para(iap+2,ip)-tab_para(ia,ip)
           zp  = tab_para(iap+1,ip)-tab_para(ia,ip)
           z0  = tab_para(iap  ,ip)-tab_para(ia,ip)
           zm  = tab_para(iap-1,ip)-tab_para(ia,ip)
           r1d(ip,ia, 2-iq,ir) =&
             (zp*z0+zp*zm+z0*zm)/&
             ((zp2-zp)*(zp2-z0)*(zp2-zm))
           r1d(ip,ia, 1-iq,ir) =&
             (zp2*z0+zp2*zm+z0*zm)/&
             ((zp-zp2)*(zp-z0)*(zp-zm))
           r1d(ip,ia,  -iq,ir) =&
             (zp2*zp+zp2*zm+zp*zm)/&
             ((z0-zp2)*(z0-zp)*(z0-zm))
           r1d(ip,ia,-1-iq,ir) =&
             (zp2*zp+zp2*z0+zp*z0)/&
             ((zm-zp2)*(zm-zp)*(zm-z0))
           r2d(ip,ia, 2-iq,ir) = 2.d00*&
             (-zp-z0-zm)/&
             ((zp2-zp)*(zp2-z0)*(zp2-zm))
           r2d(ip,ia, 1-iq,ir) = 2.d00*&
             (-zp2-z0-zm)/&
             ((zp-zp2)*(zp-z0)*(zp-zm))
           r2d(ip,ia,  -iq,ir) = 2.d00*&
             (-zp2-zp-zm)/&
             ((z0-zp2)*(z0-zp)*(z0-zm))
           r2d(ip,ia,-1-iq,ir) = 2.d00*&
             (-zp2-zp-z0)/&
             ((zm-zp2)*(zm-zp)*(zm-z0))
         end if
       else if (iq < 0) then
         if ((iap >= (ia_min+2)).and.(iap <= (ia_max-1))) then
           zp  = tab_para(iap+1,ip)-tab_para(ia,ip)
           z0  = tab_para(iap  ,ip)-tab_para(ia,ip)
           zm  = tab_para(iap-1,ip)-tab_para(ia,ip)
           zm2 = tab_para(iap-2,ip)-tab_para(ia,ip)
           r1d(ip,ia, 1-iq,ir) =&
             (z0*zm+z0*zm2+zm*zm2)/&
             ((zp-z0)*(zp-zm)*(zp-zm2))
           r1d(ip,ia,  -iq,ir) =&
             (zp*zm+zp*zm2+zm*zm2)/&
             ((z0-zp)*(z0-zm)*(z0-zm2))
           r1d(ip,ia,-1-iq,ir) =&
             (zp*z0+zp*zm2+z0*zm2)/&
             ((zm-zp)*(zm-z0)*(zm-zm2))
           r1d(ip,ia,-2-iq,ir) =&
             (zp*z0+zp*zm+z0*zm)/&
             ((zm2-zp)*(zm2-z0)*(zm2-zm))
           r2d(ip,ia, 1-iq,ir) = 2.d00*&
             (-z0-zm-zm2)/&
             ((zp-z0)*(zp-zm)*(zp-zm2))
           r2d(ip,ia,  -iq,ir) = 2.d00*&
             (-zp-zm-zm2)/&
             ((z0-zp)*(z0-zm)*(z0-zm2))
           r2d(ip,ia,-1-iq,ir) = 2.d00*&
             (-zp-z0-zm2)/&
             ((zm-zp)*(zm-z0)*(zm-zm2))
           r2d(ip,ia,-2-iq,ir) = 2.d00*&
             (-zp-z0-zm)/&
             ((zm2-zp)*(zm2-z0)*(zm2-zm))
         end if
       end if
     end do

     ! three-point formulas
     do iq=-1,1,1
       choose_rule3: select case (iq)
       case (1)
         ir = 7
       case (0)
         ir = 6
       case (-1)
         ir = -7
       end select choose_rule3
       iap = ia-iq
       if ((iap >= (ia_min+1)).and.(iap <= (ia_max-1))) then
         zp = tab_para(iap+1,ip)-tab_para(ia,ip)
         z0 = tab_para(iap  ,ip)-tab_para(ia,ip)
         zm = tab_para(iap-1,ip)-tab_para(ia,ip)
         r1d(ip,ia, 1-iq,ir) =&
           -(z0+zm)/&
           ((zp-z0)*(zp-zm))
         r1d(ip,ia,  -iq,ir) =&
           -(zp+zm)/&
           ((z0-zp)*(z0-zm))
         r1d(ip,ia,-1-iq,ir) =&
           -(zp+z0)/&
           ((zm-zp)*(zm-z0))
         r2d(ip,ia, 1-iq,ir) = 2.d00/&
           ((zp-z0)*(zp-zm))
         r2d(ip,ia,  -iq,ir) = 2.d00/&
           ((z0-zp)*(z0-zm))
         r2d(ip,ia,-1-iq,ir) = 2.d00/&
           ((zm-zp)*(zm-z0))
       end if
     end do

     ! two-point formulas
     do iq=-1,1,1
       choose_rule2: select case (iq)
       case (1)
         ir = 8
       case (-1)
         ir = -8
       end select choose_rule2
       iap = ia-iq
       if (iq > 0) then
         if ((iap >= (ia_min)).and.(iap <= (ia_max-1))) then
           zp = tab_para(iap+1,ip)-tab_para(ia,ip)
           z0 = tab_para(iap  ,ip)-tab_para(ia,ip)
           r1d(ip,ia, 1-iq,ir) =&
             1.d00/(zp-z0)
           r1d(ip,ia,  -iq,ir) =&
             1.00/(z0-zp)
         end if
       else if (iq < 0) then
         if ((iap >= (ia_min+1)).and.(iap <= (ia_max))) then
           z0 = tab_para(iap  ,ip)-tab_para(ia,ip)
           zm = tab_para(iap-1,ip)-tab_para(ia,ip)
           r1d(ip,ia,  -iq,ir) =&
             1.d00/(z0-zm)
           r1d(ip,ia,-1-iq,ir) =&
             1.d00/(zm-z0)
         end if
       end if

     end do
   end do
 end do

 if (allocated(iarg)) deallocate(iarg)

end SUBROUTINE get_diff_rules
!***********************************************************************
SUBROUTINE init_ipl_rule()
! Stefan Typel for the CompOSE core team, version 1.01, 2016/06/16
USE compose_internal
implicit none
integer :: ir

! third order
do ir=-8,8,1
   ipl_rule(0,3,ir) = ir
   ipl_rule(1,3,ir) = ir
end do

! second order
do ir=-8,8,1
   ipl_rule(0,2,ir) = 6
end do
ipl_rule(0,2, 3) =  7
ipl_rule(0,2,-3) = -7
ipl_rule(0,2, 5) =  7
ipl_rule(0,2,-5) = -7
ipl_rule(0,2, 7) =  7
ipl_rule(0,2,-7) = -7
ipl_rule(0,2, 8) =  8
ipl_rule(0,2,-8) = -8
do ir=-8,8,1
   ipl_rule(1,2,ir) = ipl_rule(0,2,ir)
end do

! first order
do ir=-8,8,1
   ipl_rule(0,1,ir) = -8
   ipl_rule(1,1,ir) = 8
end do
ipl_rule(0,1,-3) = -8
ipl_rule(0,1,-5) = -8
ipl_rule(0,1,-7) = -8
ipl_rule(0,1,-8) = -8
ipl_rule(1,1, 3) =  8
ipl_rule(1,1, 5) =  8
ipl_rule(1,1, 7) =  8
ipl_rule(1,1, 8) =  8

return
end SUBROUTINE init_ipl_rule
!***********************************************************************
subroutine get_eos_report(iwr)
 ! Stefan Typel for the CompOSE core team, version 1.17, 2017/11/23
 use compose_internal
 use omp_lib
 implicit none
 integer iwr,icc
 integer, parameter :: dim_q=10,dim_d=16
 integer, allocatable :: iqtyp(:),iqtyq(:),iqtym(:)
 integer :: alloc_status,ierror,ierr,ip,iq,it,in,iy,iflag,icnt(0:3),ifa,itmp,&
   nqty,iphase(dim_q),nphase,nphase_miss,&
   nqtyp,np_miss,nqtyq,nq_miss,nqtym,nm_miss,&
   distr(2,-1:dim_d),ipl(1:3),idxe(2,3)
 integer :: in_lorene,count_lorene
 double precision :: minimum_enthalpy, enthalpy
 double precision :: nsat,bsat,k,kp,j,l,ksym,tmp,dfa_max(2),&
   dlnt,f1,f2,fn,en,t,nb,yq,b
 logical :: first_in
#ifdef DEBUG
 double precision :: timei,timef
#endif

#ifdef DEBUG
   timei = omp_get_wtime()
#endif

 ! error counter
 ierr = 0

 ! dimension of eos table and interpolation in parameters
 eos_dim = 0
 do ip=1,3,1
   if (dim_idx(ip) > 1) then
     eos_dim = eos_dim+1
     idx_ipl(ip) = 1
   else
     idx_ipl(ip) = 0
   end if
 end do

 ! completeness and consistency check
 icnt(0:3) = 0

 do iy=1,dim_idx(3),1
   do in=1,dim_idx(2),1
     do it=1,dim_idx(1),1
       if (idx_arg(it,in,iy,0) == 1) then
         icnt(0) = icnt(0)+1
       else
         do ip=1,3,1
           if ((idx_ex(ip) == 1).and.(idx_argx(it,in,iy,ip) /= 1)) then
             icnt(ip) = icnt(ip)+1
           end if
         end do
       end if
     end do
   end do
 end do

 if (iwr == 1) then
   write(*,*)
   if (idx_ex(1) == 1) then
     write(*,*) icnt(1),'entry/ies in eos.thermo missing'
   else
     write(*,*) ' no file eos.thermo'
     stop 2
   end if
   if (idx_ex(2) == 1) then
     write(*,*) icnt(2),'entry/ies in eos.compo missing'
   else
     write(*,*) ' no file eos.compo'
   end if
   if (idx_ex(3) == 1) then
     write(*,*) icnt(3),'entry/ies in eos.micro missing'
   else
     write(*,*) ' no file eos.micro'
   end if
   write(*,*)
   write(*,*) icnt(0),'complete entries in eos tables'
   write(*,*)
 end if

  
 ! allocation of arrays
 allocate(iqtyp(np_max),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 151
 end if

 allocate(iqtyq(nq_max),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 152
 end if

 allocate(iqtym(nm_max),stat=alloc_status)
 if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 153
 end if

 ! additional quantities in eos.thermo
 nqty = nadd_max

 ! quantities in eos.compo
 ! phases
 if (dim_q > 0) iphase(1:dim_q) = -999
 nphase = 0
 nphase_miss = 0
 if (idx_ex(2) == 1) then
   do iy=1,dim_idx(3),1
     do in=1,dim_idx(2),1
       do it=1,dim_idx(1),1
         iflag = 0
         ! if(any(idx_arg(it,in,iy,4) == iphase(1:nphase))) iflag = 1
         do ip=1,nphase,1
           if (idx_arg(it,in,iy,4) == iphase(ip)) iflag = 1
         end do
         if (iflag == 0) then
           if (nphase < dim_q) then
             nphase = nphase+1
             iphase(nphase) = idx_arg(it,in,iy,4)
           else
             nphase_miss = nphase_miss+1
           end if
         end if
       end do
     end do
   end do
 end if

 ! particles
 if (np_max > 0) iqtyp(1:np_max) = -1
 nqtyp = 0
 np_miss = 0
 if (idx_ex(2) == 1) then
   do ip=1,np_max,1
     do iy=1,dim_idx(3),1
       do in=1,dim_idx(2),1
         do it=1,dim_idx(1),1
           iflag = 0
           do iq=1,nqtyp,1
             if (idxp_compo(it,in,iy,ip) == iqtyp(iq)) iflag = 1
           end do
           if (iflag == 0) then
             if (idxp_compo(it,in,iy,ip) > -1) then
               if (nqtyp < np_max) then
                 nqtyp = nqtyp+1
                 iqtyp(nqtyp) = idxp_compo(it,in,iy,ip)
               else
                 np_miss = np_miss+1
                 write(*,*) it,in,iy,nqtyp,np_max
                 write(*,*) (iqtyp(iq),iq=1,np_max,1)
                 write(*,*) (idxp_compo(it,in,iy,iq),iq=1,np_max,1)
                 stop
               end if
             end if
           end if
         end do
       end do
     end do
   end do
 end if

 ! particle sets

 if (nq_max > 0) iqtyq(1:nq_max) = -1
 nqtyq = 0
 nq_miss = 0
 if (idx_ex(2) == 1) then
   do ip=1,nq_max,1
     do iy=1,dim_idx(3),1
       do in=1,dim_idx(2),1
         do it=1,dim_idx(1),1
           iflag = 0
           do iq=1,nqtyq,1
             if (idxq_compo(it,in,iy,ip) == iqtyq(iq)) iflag = 1
           end do
           if (iflag == 0) then
             if (idxq_compo(it,in,iy,ip) > -1) then
               if (nqtyq < nq_max) then
                 nqtyq = nqtyq+1
                 iqtyq(nqtyq) = idxq_compo(it,in,iy,ip)
               else
                 nq_miss = nq_miss+1
               end if
             end if
           end if
         end do
       end do
     end do
   end do
 end if

 ! quantities in eos.micro

 if (nm_max > 0) iqtym(1:nm_max) = -1
 nqtym = 0
 nm_miss = 0
 if (idx_ex(3) == 1) then
   do ip=1,nm_max,1
     do iy=1,dim_idx(3),1
       do in=1,dim_idx(2),1
         do it=1,dim_idx(1),1
           iflag = 0
           do iq=1,nqtym,1
             if (idx_mic(it,in,iy,ip) == iqtym(iq)) iflag = 1
           end do
           if (iflag == 0) then
             if (idx_mic(it,in,iy,ip) > -1) then
               if (nqtym < nm_max) then
                 nqtym = nqtym+1
                 iqtym(nqtym) = idx_mic(it,in,iy,ip)
               else
                 nm_miss = nm_miss+1
               end if
             end if
           end if
         end do
       end do
     end do
   end do
 end if

 if (np_miss > 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 154
 end if
 if (nq_miss > 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 155
 end if
 if (nm_miss > 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 156
 end if
 if (nphase_miss > 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 157
 end if

 call write_errors(ierr)

 nsat = 0.d00
 bsat = 0.d00
 k    = 0.d00
 kp   = 0.d00
 j    = 0.d00
 l    = 0.d00
 ksym = 0.d00

 open(unit=23,file='eos.beta',&
   status='unknown',action='write',iostat=ierror)
! create beta-equilibrium files for run with Lorene (check that entalpy is monotonic)
 open(unit=24,file='eos.nb.ns',&
   status='unknown',action='write',iostat=ierror)
 open(unit=25,file='eos.thermo.ns',&
   status='unknown',action='write',iostat=ierror)
 ! check for dependence on nb and yq
 if (inbyq == 1) then
   if (incl_l /= 1) then
     if (irpl == 0) then
       ! nuclear matter parameters at lowest temperature
       !+++++++++++++++++++++++++++++++++++++++++++++
       call get_eos_nmp(nsat,bsat,k,kp,j,l,ksym,ierr)
       !+++++++++++++++++++++++++++++++++++++++++++++
     end if
     write(23,*) '# no data available'
   else
     ! EoS of beta-equilibrated matter at lowest temperature

     ipl(1:3) = 3
     if ((irpl == 0).or.(irpl == 3)) then
       if (irpl == 0) then
         t = tab_para(1,1)
         b = val_rpl
       else
         t = val_rpl
         b = tab_para(1,1)
       end if
       ip = jmap(1)
       if (ip > 0) then
         t = tab_para(1,1)
       else
         t = val_rpl
       end if
       yq = 0.5d00
       iq = jmap(3)


       if (iq > 0) then
         if (dim_idx(iq) > 0) then
           ip = jmap(2)
           if (ip > 0) then
              minimum_enthalpy = 10.d0 ! initialize
              ! look for minimum enthalpy to cut the low density entries (if necessary)
              ! which are non-monotonic
              ! and at the highest densities, too
              in_lorene = 0
              count_lorene = 0
              first_in = .true.
              
              enthalpy_loop: do in=1,dim_idx(ip),1
                 nb = tab_para(in,ip)
                 !++++++++++++++++++++++++++++++
                 call get_eos(t,nb,yq,b,ipl,1,0,eos_thermo)

                 !++++++++++++++++++++++++++++++
                 if (yq > 0.d00) then
                    enthalpy = eos_thermo(7) + 1.d0 + eos_thermo(1)/(m_n*nb)
                    !write(*,*) 'first: nb and enthalpy', nb, enthalpy,count_lorene,yq
                    if(first_in) then
                       if((minimum_enthalpy.ne.10.d0).and.(enthalpy.ge.minimum_enthalpy)) then
                          in_lorene = in - 1
                          count_lorene = count_lorene + 2
                          first_in = .false.
                       end if
                       minimum_enthalpy = enthalpy

                    else
                       if(enthalpy.gt.minimum_enthalpy) then
                          !write(*,*) 'nb and enthalpy', nb, enthalpy,count_lorene,yq
                          count_lorene = count_lorene + 1
                          minimum_enthalpy = enthalpy
                       end if
                    end if

                 end if
              end do enthalpy_loop

              if(in_lorene.gt.0) then
                 write(24,*) 1
                 write(24,*) count_lorene
                 write(25,*) m_n,m_p,1
              end if
              count_lorene = 0
              minimum_enthalpy = 0.d0
              
              do in=1,dim_idx(ip),1
                 nb = tab_para(in,ip)

                 !++++++++++++++++++++++++++++++
                 call get_eos(t,nb,yq,b,ipl,1,0,eos_thermo)
                 !++++++++++++++++++++++++++++++
                 if (yq > 0.d00) then
                    enthalpy = eos_thermo(7) + 1.d0 + eos_thermo(1)/(m_n*nb)
!                    write(*,*) 'enthalpy, density',real(enthalpy),real(nb)
                    if(in.ge.in_lorene) then
                       if(enthalpy.gt.minimum_enthalpy) then
                          count_lorene = count_lorene + 1
                          write(24,*) nb
                          write(25,*) 0,count_lorene,1,eos_thermo(1)/nb,&
                               eos_thermo(2),eos_thermo(3)/m_n,eos_thermo(4)/m_n,&
                               eos_thermo(5)/m_n,eos_thermo(6),eos_thermo(7),2,yq,enthalpy
                          minimum_enthalpy = enthalpy
                       end if
                    end if
                    write(23,*) nb,yq,(eos_thermo(6)+1.d00)*m_n*nb,eos_thermo(1)
               else
                 write(*,*) 'warning: no beta-equilibrium found for ',&
                   'baryon density n_b=',nb,'fm^-3'
               end if
             end do
             write(*,*)
             write(*,*) ' file eos.beta written'
             write(*,*)
           end if
         end if
       else
         write(*,*)
         write(*,*) 'preparation of file eos.beta not supported'
         write(*,*)
       end if
     else
       write(*,*)
       write(*,*) 'preparation of file eos.beta not supported'
       write(*,*)
     end if
   end if
 end if



 close (unit=23)
 close (unit=24)
 close (unit=25)

 ! thermodynamic consistency check and distribution of errors
 ! delta(f/n)/(f/n) and delta(e/n)/(e/n)
 icc = 1
 if (icc == 1) then
   dlnt = dlog(10.d00)
   ifa = 0
   dfa_max(1:2) = 0.d00
   if (dim_d > 0) distr(1:2,0:dim_d) = 0
   if (incl_l == 1) then
     f1 = 1.d00
     f2 = 0.d00
   else
     f1 = 0.d00
     f2 = 1.d00
   end if
   if (idx_ex(1) == 1) then
     do iy=1,dim_idx(3),1
       do in=1,dim_idx(2),1
         do it=1,dim_idx(1),1
           if (idx_arg(it,in,iy,0) == 1) then
             ifa = ifa+1
             fn = m_n*(tab_thermo(it,in,iy,3)&
               +tab_para(iy,3)*(f1*tab_thermo(it,in,iy,5)&
               +f2*tab_thermo(it,in,iy,4)))&
               -tab_thermo(it,in,iy,1)
             en = fn+tab_para(it,1)*tab_thermo(it,in,iy,2)
             fn = m_n*tab_thermo(it,in,iy,6)-fn
             fn = fn/(m_n*(tab_thermo(it,in,iy,6)+1.d00))
             en = m_n*tab_thermo(it,in,iy,7)-en
             en = en/(m_n*(tab_thermo(it,in,iy,7)+1.d00))
             do ip=1,2,1
               if (ip == 1) then
                 tmp = fn
               else
                 tmp = en
               end if
               if (dabs(tmp) > dabs(dfa_max(ip))) then
                 dfa_max(ip) = tmp
                 idxe(ip,1) = it
                 idxe(ip,2) = in
                 idxe(ip,3) = iy
               end if
               tmp = dabs(tmp)
               if (tmp < 1.d00) then
                 if (tmp < 1.d-18) then
                   itmp = dim_d
                 else
                   itmp = int(-dlog(tmp)/dlnt)
                   if (itmp > dim_d) itmp = dim_d
                 end if
               else
                 itmp = -1
               end if
               distr(ip,itmp) = distr(ip,itmp)+1
             end do
           end if
         end do
       end do
     end do
   end if
 end if


 open(unit=21,file='eos.report',&
   status='unknown',action='write',iostat=ierror)

 write(21,*) '# dimension of eos table'
 write(21,*) eos_dim
 do ip=1,4,1
   if (ip == 1)&
     write(21,*) '# number of grid points in temperature T,',&
     'minimum temperature, maximum temperature [MeV]'
   if (ip == 2)&
     write(21,*) '# number of grid points in baryon number density n_b,',&
     'minimum baryon number density, maximum baryon number density [fm^-3]'
   if (ip == 3)&
     write(21,*) '# number of grid points in charge fraction Y_q,',&
     'minimum charge fraction, maximum charge fraction [dimensionless]'
   if (ip == 4)&
     write(21,*) '# number of grid points in magnetic field strength B,',&
     'minimum magnetic field, maximum magnetic field [Gauss]'
   write(21,'(3x,i4,3x,e16.8,3x,e16.8)') dim_p(ip),&
     para_min(ip),para_max(ip)
 end do
 do ip=1,3,1
   if (ip == 1)&
     write(21,*) '# thermodynamical data exist (1) or not (0)'
   if (ip == 2)&
     write(21,*) '# composition data exist (1) or not (0)'
   if (ip == 3)&
     write(21,*) '# microscopic data exist (1) or not (0)'
   write(21,*) idx_ex(ip)
 end do
 write(21,*) '# total number of table entries'
 write(21,*) icnt(0)
 write(21,*) '# number of additional quantities in eos.thermo'
 write(21,*) nqty
 write(21,*) '# number of phases'
 write(21,*) nphase
 write(21,*) '# indices of phases (see data sheet)'
 if (nphase > 0) then
   write(21,*) (iphase(ip),ip=1,nphase,1)
 else
   write(21,*)
 end if
 write(21,*) '# number of particles'
 write(21,*) nqtyp
 write(21,*) '# indices of particles (see manual)'
 if (nqtyp > 0) then
   write(21,*) (iqtyp(ip),ip=1,nqtyp,1)
 else
   write(21,*)
 end if
 write(21,*) '# number of particle sets'
 write(21,*) nqtyq
 write(21,*) '# indices of particle sets (see manual and data sheet)'
 if (nqtyq > 0) then
   write(21,*) (iqtyq(ip),ip=1,nqtyq,1)
 else
   write(21,*)
 end if
 write(21,*) '# number of quantities in eos.micro'
 write(21,*) nqtym
 write(21,*) '# indices of quantities (see manual)'
 if (nqtym > 0) then
   write(21,*) (iqtym(ip),ip=1,nqtym,1)
 else
   write(21,*)
 end if
 write(21,*) '# saturation density of symmetric nuclear matter [fm^-3]'
 write(21,'(f10.5)') nsat
 write(21,*) '# binding energy B_sat per baryon at saturation [MeV]'
 write(21,'(f10.3)') bsat
 write(21,*) '# incompressibility K [MeV]'
 write(21,'(f10.3)') k
 write(21,*) '# skewness K^prime [MeV]'
 write(21,'(f10.3)') kp
 write(21,*) '# symmetry energy J [MeV]'
 write(21,'(f10.3)') j
 write(21,*) '# symmetry energy slope parameter L [MeV]'
 write(21,'(f10.3)') l
 write(21,*) '# symmetry incompressibility K_sym [MeV]'
 write(21,'(f10.3)') ksym

 if (icc == 1) then
   tmp = dble(ifa)
   open(unit=22,file='eos.errdistr',&
     status='unknown',action='write',iostat=ierror)
   if (ifa > 0) then
     ! for xmgrace file
     write(22,*) '@target G0.S0'
     write(22,*) '@type xy'
     do iq=0,dim_d,1
       fn = dble(distr(1,iq))/tmp
       if (fn < 1.d-30) fn = 1.d-30
       write(22,*) 10.d00**(-iq),fn
     end do
     write(22,*) 10.d00**(-dim_d-1),fn
     write(22,*) '&'
     write(22,*) '@target G0.S1'
     write(22,*) '@type xy'
     do iq=0,dim_d,1
       en = dble(distr(2,iq))/tmp
       if (en < 1.d-30) en = 1.d-30
       write(22,*) 10.d00**(-iq),en
     end do
     write(22,*) 10.d00**(-dim_d-1),en
     write(22,*) '&'
   else
     write(*,*) 'error in error distribution'
     stop
   end if
   close(22)
   write(21,*) '# maximum relative error in free energy per baryon at point'
   write(21,'(e16.8,i6,i4,i4)') dfa_max(1),idxe(1,1),idxe(1,2),idxe(1,3)
   write(21,*) '# maximum relative error in internal energy per baryon at point'
   write(21,'(e16.8,i6,i4,i4)') dfa_max(2),idxe(2,1),idxe(2,2),idxe(2,3)
 end if

 close(unit=21)

 call write_errors(ierr)

 write(*,*) ' file eos.report written'

 if (allocated(iqtyp)) deallocate(iqtyp)
 if (allocated(iqtyq)) deallocate(iqtyq)
 if (allocated(iqtym)) deallocate(iqtym)

#ifdef DEBUG
 timef = omp_get_wtime()
 print *
 print *,'TIMING: get_eos_report',timef-timei,'(s)'
#endif

end SUBROUTINE get_eos_report
!***********************************************************************
SUBROUTINE get_eos_nmp(nsat,bsat,k,kp,j,l,ksym,ierr)
! Stefan Typel for the CompOSE core team, version 1.06, 2017/11/23
USE compose_internal
implicit none
integer :: i,idx,ierr,irule,ipl(1:3),inmp
double precision :: t,nb,y,f,fp,nsat,bsat,k,kp,j,l,ksym,&
     nb_min,nb_max,f_min,f_max,d_nb,d_f,y_min,y_max,&
     g1(-2:2),g2(-2:2),g3(-2:2),d_b,d_b2

ierr = 0
inmp = 1

ipl(1:3) = 3

irule = 5

! determination of saturation density
! find intervall around minimum in free energy
t = tab_para(1,1)
y = 0.5d00
arg(1) = t
arg(3) = y
arg2(1) = t
arg2(3) = y

idx = -1
fp = 1.d30
do i=dim_idx(2),1,-1
   nb = tab_para(i,2)
   arg(2) = nb
   arg2(2) = nb
   call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
   f = eos_thermo(6)*m_n
   if (f < fp) then
      fp = f
      idx = i
   end if
 end do



if ((1 < idx).and.(dim_idx(2) > idx)) then
   nb_min = tab_para((idx-1),2)
   arg(2) = nb_min
   arg2(2) = nb_min
   call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
   f_min = eos_thermo(1)
   nb_max = tab_para((idx+1),2)
   arg(2) = nb_max
   arg2(2) = nb_max
   call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
   f_max = eos_thermo(1)
   ! determination of zero in pressure
   nb = 0.d00
   d_nb = 1.d00
   do while (dabs(d_nb) > 1.d-12)
      d_nb = nb
      d_f = f_max-f_min
      nb = 0.5d00*(nb_max+nb_min)
      arg(2) = nb
      arg2(2) = nb
      call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
      f = eos_thermo(1)
      d_nb = d_nb-nb
      if (f_min*f > 0.d00) then
         nb_min = nb
         f_min = f
         else
            nb_max = nb
            f_max = f
         end if
      end do

      ! parameters for symmetric nuclear matter
      nsat = nb
      arg(2) = nb
      arg2(2) = nb
      bsat = -(eos_thermo(6)*m_n+0.5d00*(m_n-m_p))
      k = 9.d00*nsat*eos_thermo_add(1)
      kp = 27.d00*nsat*nsat*eos_thermo_add(2)-6.d00*k
      ! isospin related parameters
      y_min = tab_para(1,3)
      y_max = tab_para(dim_idx(3),3)
      if ((y_min < 0.48d00).and.(y_max > 0.52d00)) then
         g1(0) = eos_thermo(6)*m_n
         g2(0) = eos_thermo(1)
         g3(0) = eos_thermo_add(1)
         d_b = 0.01d00
         d_b2 = d_b*d_b
         if (irule == 3) then
            y = 0.5d00-d_b
            arg(3) = y
            arg2(3) = y
            call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
            g1(-1) = eos_thermo(6)*m_n
            g2(-1) = eos_thermo(1)
            g3(-1) = eos_thermo_add(1)
            y = 0.5d00+d_b
            arg(3) = y
            arg2(3) = y
            call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
            g1(1) = eos_thermo(6)*m_n
            g2(1) = eos_thermo(1)
            g3(1) = eos_thermo_add(1)
            j = 0.125d00*(g1(-1)-2.d00*g1(0)+g1(1))/d_b2
            l = 0.375d00*(g2(-1)-2.d00*g2(0)+g2(1))/(d_b2*nsat)
            ksym = 1.125d00*nsat*(g3(-1)-2.d00*g3(0)+g3(1))/d_b2-3.d00*l
         else
            y = 0.5d00-2.d00*d_b
            arg(3) = y
            arg2(3) = y
            call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
            g1(-2) = eos_thermo(6)*m_n
            g2(-2) = eos_thermo(1)
            g3(-2) = eos_thermo_add(1)
            y = 0.5d00-d_b
            arg(3) = y
            arg2(3) = y
            call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
            g1(-1) = eos_thermo(6)*m_n
            g2(-1) = eos_thermo(1)
            g3(-1) = eos_thermo_add(1)
            y = 0.5d00+d_b
            arg(3) = y
            arg2(3) = y
            call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
            g1(1) = eos_thermo(6)*m_n
            g2(1) = eos_thermo(1)
            g3(1) = eos_thermo_add(1)
            y = 0.5d00+2.d00*d_b
            arg(3) = y
            arg2(3) = y
            call get_eos_sub_new(ipl,ierr,inmp,0,0,eos_thermo)
            g1(2) = eos_thermo(6)*m_n
            g2(2) = eos_thermo(1)
            g3(2) = eos_thermo_add(1)
            j = (-g1(-2)+16.d00*g1(-1)-30.d00*g1(0)+16.d00*g1(1)-g1(2))&
                 /(96.d00*d_b2)
            l = (-g2(-2)+16.d00*g2(-1)-30.d00*g2(0)+16.d00*g2(1)-g2(2))&
                 /(32.d00*d_b2*nsat)
            ksym = 3.d00*nsat*&
                 (-g3(-2)+16.d00*g3(-1)-30.d00*g3(0)+16.d00*g3(1)-g3(2))&
                 /(32.d00*d_b2)-3.d00*l
         endif
      end if
   end if

return
end SUBROUTINE get_eos_nmp
!***********************************************************************
SUBROUTINE define_eos_table(iwr)
! Stefan Typel for the CompOSE core team, version 1.07, 2017/05/22
USE compose_internal
implicit none
integer :: i,j,ierr,iunit,ierror,iwr

n_add = 0
n_p = 0
n_q = 0
n_m = 0

if (idx_ex(2) /= 1) then
   np_max = 0
   nq_max = 0
end if

if (idx_ex(3) /= 1) then
   nm_max = 0
end if

ierr = 0
iunit = 20
open(unit=iunit,file='eos.quantities',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   if (iwr == 1) then
      write(*,*)
      write(*,*) ' reading selection of quantities'
   end if

   read(iunit,*)
   read(iunit,*) n_qty,n_add,n_df
   if ((n_qty < 0).or.(n_qty > dim_qtyt)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 61
   end if
   if ((n_add < 0).or.(n_add > nadd_max)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 62
   end if
!2017/05/22
   if ((n_df < 0).or.(n_df > dim_df)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 75
   end if
   if (ierr == 0) then
      read(iunit,*)
      read(iunit,*) (idx_qty(i),i=1,n_qty,1),(idx_add(j),j=1,n_add,1),&
           (idx_df(j),j=1,n_df,1)
      do i=1,n_qty,1
         if ((idx_qty(i) < 0).or.(idx_qty(i) > dim_qtyt)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 63
         end if
      end do
      do i=1,n_add,1
         if ((idx_add(i) < 0).or.(idx_add(i) > nadd_max)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 64
         end if
      end do

      read(iunit,*)
      read(iunit,*) n_p,n_q
      if (idx_ex(2) /= 1) then
!         n_p = 0
!         n_q = 0
         if ((n_p /= 0).or.(n_q /=0)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 161
         else
            n_p = 0
            n_q = 0
         end if
      end if
      if ((n_p < 0).or.(n_p > np_max)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 65
         write(*,*) 'n_p',n_p
      end if
      if ((n_q < 0).or.(n_q > nq_max)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 66
      end if
      if (ierr == 0) then
         read(iunit,*)
         read(iunit,*) (idx_p(i),i=1,n_p,1),(idx_q(j),j=1,n_q,1)

         read(iunit,*)
         read(iunit,*) n_m
         if (idx_ex(3) /= 1) then
!            n_m = 0
            if (n_m /= 0) then
               ierr = ierr+1
               if (ierr < dim_err) error_msg(ierr) = 162
            else
               n_m = 0
            end if
         end if
         if ((n_m < 0).or.(n_m > nm_max)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 67
         end if
         if (ierr == 0) then
            read(iunit,*)
            read(iunit,*) (idx_m(i),i=1,n_m,1)

            read(iunit,*)
            read(iunit,*) n_err
            if ((n_err < 0).or.(n_err > 8)) then
               ierr = ierr+1
               if (ierr < dim_err) error_msg(ierr) = 68
            end if
            if (ierr == 0) then
               read(iunit,*)
               read(iunit,*) (idx_err(i),i=1,n_err,1)
               if ((idx_err(i) < 1).or.(idx_err(i) > 8)) idx_err(i) = i

               read(iunit,*)
               read(iunit,*) iout
#ifndef have_hdf5
               if ((iout < 1).or.(iout > 2)) then
                  ierr = ierr+1
                  if (ierr < dim_err) error_msg(ierr) = 69
               end if
#endif
            end if
         end if

      end if
      if (n_add > nadd_max) n_add = nadd_max
      if (n_p > np_max) n_p = np_max
      if (n_q > nq_max) n_q = nq_max
      if (n_m > nm_max) n_m = nm_max

   end if
else
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 60
endif
close(unit=iunit)

call write_errors(ierr)


if (iwr == 1) then
   write(*,*)
   write(*,*) ' number of extracted regular quantities    =',n_qty
   write(*,*) ' indices of regular quantities             :',&
        (idx_qty(j),j=1,n_qty,1)
   write(*,*)
   write(*,*) ' number of extracted additional quantities =',n_add
   write(*,*) ' indices of additional quantities          :',&
        (idx_add(j),j=1,n_add,1)
!2017/05/22
   write(*,*)
   write(*,*) ' number of extracted derivative quantities =',n_df
   write(*,*) ' indices of derivative quantities          :',&
        (idx_df(j),j=1,n_df,1)
   write(*,*)
   write(*,*) ' number of pairs for composition           =',n_p
   write(*,*) ' indices of pairs                          :',&
        (idx_p(j),j=1,n_p,1)
   write(*,*)
   write(*,*) ' number of quadruples for composition      =',n_q
   write(*,*) ' indices of quadruples                     :',&
        (idx_q(j),j=1,n_q,1)
   write(*,*)
   write(*,*) ' number of microscopic quantities          =',n_m
   write(*,*) ' indices of microscopic quantities         :',&
        (idx_m(j),j=1,n_m,1)
   write(*,*)
   write(*,*) ' number of error quantities                =',n_err
   write(*,*) ' indices of error quantities               :',&
        (idx_err(j),j=1,n_err,1)
   write(*,*)
   if (iout == 1) then
      write(*,*) ' format of output table                    = ',&
           iout,'(ASCII)'
   else
      write(*,*) ' format of output table                    = ',&
           iout,'(HDF5)'
   end if
end if

return
end SUBROUTINE define_eos_table
!***********************************************************************
subroutine get_eos_table(iwr)
 ! Stefan Typel for the CompOSE core team, version 1.19, 2017/12/13
 use general_var, only: tabulation_schema, pts, tabMin, tabMax
 use compose_internal
 use m_out_to_json
#if defined have_hdf5
 use hdfparameters
#endif

 implicit none
 integer :: i,jj, i_tny,i_beta,ipl_t,ipl_n,ipl_y,itest,iwr,&
   i1,i2,i3,i4,i5,i6,i7,j_t,j_nb,j_yq,&
   i_t,i_nb,i_yq,i_entr,&
   j_b,i_b,ierr,iunit,ierror,iunit2,ipl(3)

 double precision :: d_t,d_nb,d_yq,d_b,t,n,y,b

 if (iwr == 1) then
   write(*,*)
   write(*,*) ' begin generating eos table'
 end if

 itest = 0

 ierr = 0

 if (itest == 0) then

   iunit = 20
   open(unit=iunit,file='eos.parameters',&
     status='old',action='read',iostat=ierror)
   if (ierror == 0) then
     if (iwr ==1) then
       write(*,*)
       write(*,*) ' reading parameters'
     end if

     read(iunit,*)
     read(iunit,*) ipl_t,ipl_n,ipl_y
     if ((ipl_t < 1).or.(ipl_t > 3)) then
       ierr = ierr+1
       if (ierr < dim_err) error_msg(ierr) = 70
     end if
     if ((ipl_n < 1).or.(ipl_n > 3)) then
       ierr = ierr+1
       if (ierr < dim_err) error_msg(ierr) = 71
     end if
     if ((ipl_y < 1).or.(ipl_y > 3)) then
       ierr = ierr+1
       if (ierr < dim_err) error_msg(ierr) = 72
     end if

     ipl(1) = ipl_t
     ipl(2) = ipl_n
     ipl(3) = ipl_y

     read(iunit,*)
     read(iunit,*) i_beta,i_entr
     if (i_beta == 1) then
       if ((imap(3) /= 3).or.(incl_l /= 1).or.(dim_idx(3) < 2)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 160
       end if
     end if
     if (i_entr /= 1) i_entr = 0

     if (iwr == 1) then
       if (i_beta == 1) then
         write(*,*)
         write(*,*) ' EoS table for condition of beta equilibrium'
       end if
       if (i_entr == 1) then
         write(*,*)
         write(*,*) ' EoS table for given entropy per baryon'
       end if
     end if

     call write_errors(ierr)

     ! tabulation scheme is red here : 0 explicit , 1 loop
     read(iunit,*)
     read(iunit,*) tabulation_schema
     read(iunit,*)
     if (tabulation_schema == 0) then
       read(iunit,*) pts%tnyb
       if (pts%tnyb < 1) then
         ierr = ierr +1
         if (ierr < dim_err) error_msg(ierr) = 80
       end if
     else

       if (irpl == 0) then
         read(iunit,*) tabMin%t,tabMin%nb,tabMin%yq
         read(iunit,*) tabMax%t,tabMax%nb,tabMax%yq
         read(iunit,*) pts%t,pts%nb,pts%yq
         read(iunit,*) i_t,i_nb,i_yq
       else
         read(iunit,*) tabMin%t,tabMin%nb,tabMin%yq,tabMin%b
         read(iunit,*) tabMax%t,tabMax%nb,tabMax%yq,tabMax%b
         read(iunit,*) pts%t,pts%nb,pts%yq,pts%b
         read(iunit,*) i_t,i_nb,i_yq,i_b
       end if

       if (i_t == 0) then
       else
         if (tabMin%t <= 0.d00) then
           ierr = ierr +1
           if (ierr < dim_err) error_msg(ierr) = 82
         end if
       end if
       if (tabMax%t < tabMin%t) then
         ierr = ierr +1
         if (ierr < dim_err) error_msg(ierr) = 83
       end if
       if (pts%t < 1) pts%t = 1
       if (i_t == 0) then
         if (pts%t == 1) then
           d_t = 0.d00
         else
           d_t = (tabMax%t-tabMin%t)/dble(pts%t-1)
         end if
       else
         if (pts%t == 1) then
           d_t = 1.d00
         else
           d_t = exp(log(tabMax%t/tabMin%t)/dble(pts%t-1))
         end if
       end if

       if (tabMax%nb < tabMin%nb) then
         ierr = ierr +1
         if (ierr < dim_err) error_msg(ierr) = 85
       end if
       if (pts%nb < 1) pts%nb = 1
       if (i_nb == 0) then
         if (pts%nb == 1) then
           d_nb = 0.d00
         else
           d_nb = (tabMax%nb-tabMin%nb)/dble(pts%nb-1)
         end if
       else
         if (pts%nb == 1) then
           d_nb = 1.d00
         else
           d_nb = exp(log(tabMax%nb/tabMin%nb)/dble(pts%nb-1))
         end if
       end if

       if (i_beta /= 0) then
         tabMin%yq = 0.d00
         tabMax%yq = 0.d00
         d_yq = 0.d00
         pts%yq = 1
       else

         if (tabMax%yq < tabMin%yq) then
           ierr = ierr +1
           if (ierr < dim_err) error_msg(ierr) = 87
         end if
         if (tabMax%yq > 1.d00) then
           ierr = ierr +1
           if (ierr < dim_err) error_msg(ierr) = 88
         end if
         if (pts%yq < 1) pts%yq = 1
         if (i_yq == 0) then
           if (pts%yq == 1) then
             d_yq = 0.d00
           else
             d_yq = (tabMax%yq-tabMin%yq)/dble(pts%yq-1)
           end if
         else
           if (pts%yq == 1) then
             d_yq = 1.d00
           else
             d_yq = exp(log(tabMax%yq/tabMin%yq)/dble(pts%yq-1))
           end if
         end if
       end if

       if (irpl > 0) then

         if (tabMax%b < tabMin%b) then
           ierr = ierr +1
           if (ierr < dim_err) error_msg(ierr) = 90
         end if
         if (pts%b < 1) pts%b = 1
         if (i_b == 0) then
           if (pts%b == 1) then
             d_b = 0.d00
           else
             d_b = (tabMax%b-tabMin%b)/dble(pts%b-1)
           end if
         else
           if (pts%b == 1) then
             d_b = 1.d00
           else
             d_b = exp(log(tabMax%b/tabMin%b)/dble(pts%b-1))
           end if
         end if
       else
         tabMin%b = 0.d00
         tabMax%b = 0.d00
         d_b = 0.d00
         pts%b = 1
         i_b = 0
       end if
     end if

     ! output in ASCII or HDF5 format
     iunit2 = 21
     if (iout == 1) then
       ! ASCII
       open(unit=iunit2,file='eos.table',&
         status='unknown',action='write',iostat=ierror)
     else
       ! HDF5
#if defined have_hdf5
       IF(tabulation_schema == 0) then
         call initialise_hdf5(pts%tnyb,1,1,i_beta,i_entr)
      else
         call initialise_hdf5(pts%nb,pts%t,pts%yq,i_beta,i_entr)
       end IF
#endif
     end if

     if (ierror == 0) then
       if (ierr == 0) then
         if (tabulation_schema == 0) then
           ! use list of parameters from file eos.parameters
           do i_tny=1, pts%tnyb
             if (irpl > 0) then
               read(iunit,*) t,n,y,b
             else
               read(iunit,*) t,n,y
             endif
             !++++++++++++++++++++++++++++++++++++++
             call get_eos(t,n,y,b,ipl,i_beta,i_entr,eos_thermo)
             !++++++++++++++++++++++++++++++++++++++
             if (i_entr > 1) then
               write(*,*) ' no solution found for S =',t,&
                 ' at n_b =',arg2(2),' fm^-3 and Y_q =',arg2(3)
               i_entr = 1
             else
#ifndef HAVE_WEB
               !In web-app version this if has to be commented
               if (y > -1.d00) then
#endif
                 ! no output if no solution of beta equilibrium
                 if (iout == 1) then
                   ! ASCII
                   do jj = 1,4
                     do i=1,n_q,1
                       eos_q(4*(i-1)+jj) = eos_compo_q(i,jj)
                     end do
                   end do
                   if (irpl > 0) then
                     write(iunit2,*) arg2(1),arg2(2),arg2(3),arg2(4),&
                       (eos_thermo(idx_qty(i1)),i1=1,n_qty,1),&
                       (eos_thermo_add(i2),i2=1,n_add,1),&
                       (eos_compo_p(i3),i3=1,n_p,1),&
                       (eos_q(i4),i4=1,4*n_q,1),&
                       (eos_micro(i5),i5=1,n_m,1),&
                       (eos_thermo_err(idx_err(i6)),i6=1,n_err,1)

                   else
                     write(iunit2,*) arg2(1),arg2(2),arg2(3),&
                       (eos_thermo(idx_qty(i1)),i1=1,n_qty,1),&
                       (eos_thermo_add(i2),i2=1,n_add,1),&
                       (eos_df(idx_df(i7)),i7=1,n_df,1),&
                       (eos_compo_p(i3),i3=1,n_p,1),&
                       (eos_q(i4),i4=1,4*n_q,1),&
                       (eos_micro(i5),i5=1,n_m,1),&
                       (eos_thermo_err(idx_err(i6)),i6=1,n_err,1)
                   end if
                 else
                   ! HDF5
#if defined have_hdf5
                   t_hdf5(i_tny) = arg2(1)
                   nb_hdf5(i_tny) = arg2(2)
                   y_q_hdf5(i_tny) = arg2(3)
                   if(n_qty.ne.0) then
                     thermo_hdf5(i_tny,1,1,1:n_qty) = eos_thermo(idx_qty(1:n_qty))
                     index_thermo(1:n_qty) = idx_qty(1:n_qty)
                   end if
                   if(n_add.ne.0) then
                     thermo_hdf5_add(i_tny,1,1,1:n_add) = eos_thermo_add(1:n_add)
                     index_thermo_add(1:n_add) = idx_add(1:n_add)
                   end if
                   if(n_df.ne.0) then
                     deriv_hdf5(i_tny,1,1,1:n_df) = eos_df(idx_df(1:n_df))
                     index_deriv(1:n_df) = idx_df(1:n_df)
                   end if
                   if(n_p.ne.0) then
                     yi_hdf5(i_tny,1,1,1:n_p) = eos_compo_p(1:n_p)
                     index_yi(1:n_p) = idx_p(1:n_p)
                   end if
                   if(n_q.ne.0) then
                     !2017/10/09 ->
                     yav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,1)
                     aav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,2)
                     zav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,3)
                     nav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,4)
                     !2017/10/09 <-
                     index_av(1:n_q) = idx_q(1:n_q)
                   end if
                   if(n_m.ne.0) then
                     micro_hdf5(i_tny,1,1,1:n_m) = eos_micro(1:n_m)
                     index_micro(1:n_m) = idx_m(1:n_m)
                   end if
                   if(n_err.ne.0) then
                     err_hdf5(i_tny,1,1,1:n_err) = eos_thermo_err(idx_err(1:n_err))
                     index_err(1:n_err) = idx_err(1:n_err)
                   end if
#endif
                 end if
#ifndef HAVE_WEB
               end if
#endif
             end if
           end do
         else
           ! use cycle form of parameters from file eos.parameters

           ! !$OMP PARALLEL DO &
           ! !$OMP DEFAULT(SHARED) &
           ! !$OMP PRIVATE(j_t,t,j_nb,n,j_yq,y,j_b,b,arg,arg2)

           do j_t=1,pts%t,1

             if (i_t == 0) then
               t = tabMin%t + d_t * dble(j_t-1)
             else
               t = tabMin%t * (d_t**(j_t-1))
             end if

             do j_nb=1, pts%nb
               if (i_nb == 0) then
                 n = tabMin%nb + d_nb * dble(j_nb-1)
               else
                 n = tabMin%nb * (d_nb**(j_nb-1))
               end if

               do j_yq=1, pts%yq
                 if (i_yq == 0) then
                   y = tabMin%yq + d_yq * dble(j_yq-1)
                 else
                   y = tabMin%yq * (d_yq**(j_yq-1))
                 end if

                 do j_b=1, pts%b
                   if (i_b == 0) then
                     b = tabMin%b + d_b * dble(j_b-1)
                   else
                     b = tabMin%b * (d_b**(j_b-1))
                   end if
#if defined have_hdf5
                   ! magnetic field not yet implemented for HDF5 ou,arg3tput
#endif
                   !++++++++++++++++++++++++++++++++++++++
                   call get_eos(t,n,y,b,ipl,i_beta,i_entr,eos_thermo)

                   !++++++++++++++++++++++++++++++++++++++
                   if (i_entr > 1) then
                     write(*,*) ' no solution found for S =',t,&
                       ' at n_b =',arg(2),' fm^-3 and Y_q =',arg(3)
                     i_entr = 1
                   else
#ifndef HAVE_WEB
                     ! In web-app version this if has to be commented
                     if (y > -1.d00) then
#endif
                       ! no output if no solution of beta equilibrium
                       if (iout == 1) then
                         !ASCII
                         do i=1,n_q,1
                           eos_q(4*(i-1)+1) = eos_compo_q(i,1)
                           eos_q(4*(i-1)+2) = eos_compo_q(i,2)
                           eos_q(4*(i-1)+3) = eos_compo_q(i,3)
                           eos_q(4*(i-1)+4) = eos_compo_q(i,4)
                         end do
                         if (irpl > 0) then
                           write(iunit2,*) arg2(1),arg2(2),arg2(3),arg2(4),&
                             (eos_thermo(idx_qty(i1)),i1=1,n_qty,1),&
                             (eos_thermo_add(i2),i2=1,n_add,1),&
                             (eos_compo_p(i3),i3=1,n_p,1),&
                             (eos_q(i4),i4=1,4*n_q,1),&
                             (eos_micro(i5),i5=1,n_m,1),&
                             (eos_thermo_err(idx_err(i6)),i6=1,n_err,1)
                         else
                           write(iunit2,*) arg2(1),arg2(2),arg2(3),&
                             (eos_thermo(idx_qty(i1)),i1=1,n_qty,1),&
                             (eos_thermo_add(i2),i2=1,n_add,1),&
                             (eos_df(idx_df(i7)),i7=1,n_df,1),&
                             (eos_compo_p(i3),i3=1,n_p,1),&
                             (eos_q(i4),i4=1,4*n_q,1),&
                             (eos_micro(i5),i5=1,n_m,1),&
                             (eos_thermo_err(idx_err(i6)),i6=1,n_err,1)
                           !stop
                         endif
                       else
                         ! HDF5
#if defined have_hdf5
                          !2017/05/23
                          if(i_entr ==0) then
                             t_hdf5(j_t) = arg2(1)
                          else
                             t_hdf5(j_nb + (j_t-1)*pts%nb) = arg2(1)
                          end if
                          nb_hdf5(j_nb) = arg2(2)
                          if(i_beta ==0) then
                             y_q_hdf5(j_yq) = arg2(3)
                          else
                             y_q_hdf5(j_nb) = arg2(3)
                          end if
                          IF(n_qty.ne.0) then
                             thermo_hdf5(j_nb,j_t,j_yq,1:n_qty) = eos_thermo(idx_qty(1:n_qty))
                             index_thermo(1:n_qty) = idx_qty(1:n_qty)
                          end IF
                          IF(n_add.ne.0) then
                             thermo_hdf5_add(j_nb,j_t,j_yq,1:n_add) = eos_thermo_add(1:n_add)
                             index_thermo_add(1:n_add) = idx_add(1:n_add)
                          end IF
                          if(n_df.ne.0) then
                             deriv_hdf5(j_nb,j_t,j_yq,1:n_df) = eos_df(idx_df(1:n_df))
                             index_deriv(1:n_df) = idx_df(1:n_df)
                          end if
                          IF(n_p.ne.0) then
                             yi_hdf5(j_nb,j_t,j_yq,1:n_p) = eos_compo_p(1:n_p)
                           index_yi(1:n_p) = idx_p(1:n_p)
                        end IF
                        IF(n_q.ne.0) then
                           !2017/10/09 ->
                           yav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,1)
                           aav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,2)
                           zav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,3)
                           nav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,4)
                           !2017/70/09 <-
                           index_av(1:n_q) = idx_q(1:n_q)
                         end IF
                         if(n_m.ne.0) then
                           micro_hdf5(j_nb,j_t,j_yq,1:n_m) = eos_micro(1:n_m)
                           index_micro(1:n_m) = idx_m(1:n_m)
                         end if
                         IF(n_err.ne.0) then
                           err_hdf5(j_nb,j_t,j_yq,1:n_err) = eos_thermo_err(idx_err(1:n_err))
                           index_err(1:n_err) = idx_err(1:n_err)
                         end IF

#endif
                       end if
#ifndef HAVE_WEB
                     end if
#endif
                   end if
                 end do
               end do
             end do
           end do
           ! !$OMP END PARALLEL DO
         end if
       end if
     end if
   end if
   if (iout == 1) then
     ! ASCII
     close(unit=iunit2)
   else
     ! HDF5
#if defined have_hdf5
     write(*,*)'call writing HDF5 table'
     call write_hdf5(i_beta,i_entr)
     call close_hdf5
#endif
   end if
   close(unit=iunit)

 else
   ! only for testing
 end if

 call write_errors(ierr)

 if (iwr == 1) then
   i2 = 0
   write(*,*)
   write(*,*) ' end generating EoS table'
   write(*,*)
   write(*,*) ' file eos.table written'
   write(*,*)
   write(*,*) ' The columns of the file eos.table contain the following quantities:'
   write(*,*)
   i2 = i2+1
   write(*,*) i2,' temperature T                                      [MeV]        '
   i2 = i2+1
   write(*,*) i2,' baryon number density n_b                          [fm^-3]      '
   i2 = i2+1
   write(*,*) i2,' hadronic charge fraction Y_q                       []           '
   if (irpl > 0) then
     i2 = i2+1
     write(*,*) i2,' magnetic field strength B                          [G]          '
   endif
   do i1=1,n_qty,1
     i2 = i2+1
     if (idx_qty(i1) == 1)&
       write(*,*) i2,' pressure p                                         [MeV fm^-3]  '
     if (idx_qty(i1) == 2)&
       write(*,*) i2,' entropy per baryon S                               []           '
     if (idx_qty(i1) == 3)&
       write(*,*) i2,' shifted baryon chemical potential mu_b-m_n         [MeV]        '
     if (idx_qty(i1) == 4)&
       write(*,*) i2,' charge chemical potential mu_q                     [MeV]        '
     if (idx_qty(i1) == 5)&
       write(*,*) i2,' lepton chemical potential mu_l                     [MeV]        '
     if (idx_qty(i1) == 6)&
       write(*,*) i2,' scaled free energy per baryon F/m_n-1              []           '
     if (idx_qty(i1) == 7)&
       write(*,*) i2,' scaled internal energy per baryon E/m_n-1          []           '
     if (idx_qty(i1) == 8)&
       write(*,*) i2,' scaled enthalpy energy per baryon H/m_n-1          []           '
     if (idx_qty(i1) == 9)&
       write(*,*) i2,' scaled free enthalpy per baryon G/m_n-1            []           '
     if (idx_qty(i1) == 10)&
       write(*,*) i2,' derivative dp/dn_b|E                               [MeV]        '
     if (idx_qty(i1) == 11)&
       write(*,*) i2,' derivative dp/dE|n_b                                [fm^-3]      '
     if (idx_qty(i1) == 12)&
       write(*,*) i2,' square of speed of sound (c_s)^2                   []           '
     if (idx_qty(i1) == 13)&
       write(*,*) i2,' specific heat capacity at constant volume c_V      []           '
     if (idx_qty(i1) == 14)&
       write(*,*) i2,' specific heat capacity at constant pressure c_p    []           '
     if (idx_qty(i1) == 15)&
       write(*,*) i2,' adiabatic index Gamma                              []           '
     if (idx_qty(i1) == 16)&
       write(*,*) i2,' expansion coefficient at constant pressure alpha_p [MeV^-1]     '
     if (idx_qty(i1) == 17)&
       write(*,*) i2,' tension coefficient at constant volume beta_V      [fm^-3]      '
     if (idx_qty(i1) == 18)&
       write(*,*) i2,' isothermal compressibility kappa_T                 [MeV^-1 fm^3]'
     if (idx_qty(i1) == 19)&
       write(*,*) i2,' isentropic compressibility kappa_S                 [MeV^-1 fm^3]'
     if (idx_qty(i1) == 20)&
       write(*,*) i2,' free energy per baryon F                           [MeV]        '
     if (idx_qty(i1) == 21)&
       write(*,*) i2,' internal energy per baryon E                       [MeV]        '
     if (idx_qty(i1) == 22)&
       write(*,*) i2,' enthalpy per baryon H                              [MeV]        '
     if (idx_qty(i1) == 23)&
       write(*,*) i2,' free enthalpy per baryon G                         [MeV]        '
     if (idx_qty(i1) == 24)&
       write(*,*) i2,' energy density epsilon                       [MeV/fm^3]        '
   end do
   if (n_add > 0) then
     do i1=1,n_add,1
       i2 = i2+1
       write(*,*) i2,' additional thermodynamic quantity #                    ',idx_add(i1)
     end do
   end if
   !2017/05/22
   if (n_df > 0) then
     do i1=1,n_df,1
       i2 = i2+1
       write(*,*) i2,' derivative quantity #                                  ',idx_df(i1)
     end do
   end if
   if (n_p > 0) then
     do i1=1,n_p,1
       i2 = i2+1
       write(*,*) i2,' number fraction Y of particle with index               ',idx_p(i1)
     end do
   end if
   if (n_q > 0) then
     do i1=1,n_q,1
       i2 = i2+1
       write(*,*) i2,' total number fraction Y of particle set with index     ',idx_q(i1)
       i2 = i2+1
       write(*,*) i2,' average mass number A_av of particle set with index    ',idx_q(i1)
       i2 = i2+1
       write(*,*) i2,' average proton number Z_av of particle set with index  ',idx_q(i1)
       i2 = i2+1
       write(*,*) i2,' average neutron number N_av of particle set with index ',idx_q(i1)
     end do
   end if
   if (n_m > 0) then
     do i1=1,n_m,1
       i2 = i2+1
       write(*,*) i2,' microscopic quantity with index                        ',idx_m(i1)
     end do
   end if
   if (n_err > 0) then
     do i1=1,n_err,1
       i2 = i2+1
       write(*,*) i2,' error quantity with index                              ',idx_err(i1)
     end do
   end if
   write(*,*)
 end if

 ! write output on json file
 call out_to_json_write(iwr,irpl,i_beta, &
   &                    idx_qty, idx_add, idx_df, idx_p, idx_q, idx_m, idx_err,&
   &                    n_qty,   n_add,    n_df,  n_p,   n_q,   n_m,&
   &                    n_err)


 if (allocated(tab_para)) deallocate(tab_para)
 if (allocated(idx_arg)) deallocate(idx_arg)
 if (allocated(tab_thermo)) deallocate(tab_thermo)
 if (allocated(idx_thermo)) deallocate(idx_thermo)
 if (allocated(v_thermo)) deallocate(v_thermo)
 if (allocated(eos_thermo_add)) deallocate(eos_thermo_add)
 if (allocated(idx_add)) deallocate(idx_add)
 if (allocated(idxp_compo)) deallocate(idxp_compo)
 if (allocated(idxq_compo)) deallocate(idxq_compo)
 if (allocated(tabp_compo)) deallocate(tabp_compo)
 if (allocated(tabq_compo)) deallocate(tabq_compo)
 if (allocated(idx_p)) deallocate(idx_p)
 if (allocated(idx_q)) deallocate(idx_q)
 if (allocated(idx_compo_p)) deallocate(idx_compo_p)
 if (allocated(idx_compo_q)) deallocate(idx_compo_q)
 if (allocated(eos_compo_p)) deallocate(eos_compo_p)
 if (allocated(eos_compo_q)) deallocate(eos_compo_q)
 if (allocated(eos_q)) deallocate(eos_q)
 if (allocated(idx_mic)) deallocate(idx_mic)
 if (allocated(tab_mic)) deallocate(tab_mic)
 if (allocated(idx_m)) deallocate(idx_m)
 if (allocated(idx_micro)) deallocate(idx_micro)
 if (allocated(eos_micro)) deallocate(eos_micro)
 if (allocated(idx_argx)) deallocate(idx_argx)
 if (allocated(r1d)) deallocate(r1d)
 if (allocated(r2d)) deallocate(r2d)

end SUBROUTINE get_eos_table
!***********************************************************************
subroutine get_eos(t,n,y,b,ipl,i_beta,i_entr,eos_thermo)
! Stefan Typel for the CompOSE core team, version  1.13, 2018/09/07
 use compose_internal, only : arg, idx_argx, idx_ipl, imap, jmap, arg2,&
   & tab_thermo, tab_para, dim_idx, eos_dim, val_rpl, dim_qtyt
 implicit none
 integer,intent(in) :: i_beta,ipl(3)
 integer :: i_entr
 integer :: ierr
 double precision,intent(inout) :: t, n, y, b
 double precision,intent(out) :: eos_thermo(dim_qtyt)

 integer ::  i1, itmp
 double precision :: s,s_min,s_max,ds,t_min,t_max,sm,sp

 !eos_thermo = 0.d0

 ! error counter
 ierr = 0

 arg(1) = t
 arg(2) = n
 arg(3) = y

 if (eos_dim < 3) then
   do i1 = 1,3,1
     if (idx_ipl(i1) == 0) then
       arg(imap(i1)) = tab_para(1,i1)
     end if
   end do
   t = arg(1)
   n = arg(2)
   y = arg(3)
 end if

 if (jmap(1) == 0) then
   t = val_rpl
 else
   t = arg(jmap(1))
 endif
 if (jmap(2) == 0) then
   n = val_rpl
 else
   n = arg(jmap(2))
 endif
 if (jmap(3) == 0) then
   y = val_rpl
 else
   y = arg(jmap(3))
 endif
 if (jmap(4) == 0) then
   b = val_rpl
 else
   b = arg(jmap(4))
 endif

 arg2(1) = t
 arg2(2) = n
 arg2(3) = y
 arg2(4) = b

 if (i_entr == 1) then

   ! entropy
   s = arg(1)
   !write(*,*) ' s =',s
   !write(*,*) ' t_min =', tab_para(1,1)
   !write(*,*) ' t_max =', tab_para(dim_idx(1),1)

   ! boundaries

   t_min = tab_para(1,1)*1.00000001d00
   arg(1) = t_min
   arg2(1) = arg(1)
   !+++++++++++++++++++++++++++++++++++++++++++

   call get_eos_sub_s(y,ipl,ierr,i_beta,i_entr,eos_thermo)

   !+++++++++++++++++++++++++++++++++++++++++++
   s_min = eos_thermo(2)
   if (s == s_min) return
  
   t_max = tab_para(dim_idx(1),1)*0.99999999d00
   arg(1) = t_max
   arg2(1) = arg(1)
   !+++++++++++++++++++++++++++++++++++++++++++
   call get_eos_sub_s(y,ipl,ierr,i_beta,i_entr,eos_thermo)
   !+++++++++++++++++++++++++++++++++++++++++++
   s_max = eos_thermo(2)
   if (s == s_max) return

   !write(*,*) ' s_min s_max',s_min,s_max
   if ((s < s_max).and.(s > s_min)) then

     s_max = s_max-s
     s_min = s_min-s
     ds = 1.d00
     i1 = 0
     do while (dabs(ds) > 1.d-12*s)
       i1 = i1+1
       if (i1 == 50) then
         i_entr = 2
         return
       end if

       if (dabs(s_max-s_min) < s) then
         arg(1) = (s_max*t_min-s_min*t_max)/(s_max-s_min)
       else
         arg(1) = 0.5*(t_min+t_max)
       end if
       arg2(1) = arg(1)
       !+++++++++++++++++++++++++++++++++++++++++++
       call get_eos_sub_s(y,ipl,ierr,i_beta,i_entr,eos_thermo)
       !+++++++++++++++++++++++++++++++++++++++++++
       ds = eos_thermo(2)-s
       !write(*,*) i1,arg(1),eos_thermo(2),ds
       if (ds > 0.d00) then
         t_max = arg(1)
         s_max = ds
       else if (ds < 0.d00) then
         t_min = arg(1)
         s_min = ds
       else
         return
       end if
     end do
   else
     i_entr = 2
   end if
 else

   !+++++++++++++++++++++++++++++++++++++++++++
   call get_eos_sub_s(y,ipl,ierr,i_beta,i_entr,eos_thermo)
   !+++++++++++++++++++++++++++++++++++++++++++

 end if
 call write_errors(ierr)

end subroutine get_eos
!***********************************************************************
subroutine get_eos_sub_s(y,ipl,ierr,i_beta,i_entr,eos_thermo)
! Stefan Typel for the CompOSE core team, version 1.00, 2017/11/23
use compose_internal, only : arg2, dim_qtyt
implicit none
integer,intent(in) :: i_beta, i_entr, ipl(3)
integer,intent(out) :: ierr
double precision,intent(in) :: y
double precision,intent(out) :: eos_thermo(dim_qtyt)

if (i_beta == 1) then
   ! determination of beta equilibrium
   !+++++++++++++++++++++++++++++++++++
   call get_eos_beta(y,ipl,ierr,i_entr, eos_thermo)
   !+++++++++++++++++++++++++++++++++++
   arg2(3) = y
else
   !++++++++++++++++++++++++++++++++++++
   call get_eos_sub_new(ipl,ierr,0,0,i_entr, eos_thermo)
   !++++++++++++++++++++++++++++++++++++
end if

return
end SUBROUTINE get_eos_sub_s
!***********************************************************************
subroutine get_eos_beta(y,ipl,ierr,i_entr, eos_thermo)
 ! Stefan Typel for the CompOSE core team, version 1.08, 2017/11/23
 use compose_internal, only:  arg, tab_para, dim_idx, dim_qtyt
 use omp_lib
 implicit none
 integer :: ierr,i_entr,ipl(3),itmp,i3
 double precision ,intent(out):: eos_thermo(dim_qtyt)
 double precision :: y,y_min,y_max,f_min,f_max,f,d_f,d_y,fm,fp
#ifdef DEBUG
 double precision :: timei,timef
#endif

 ! error counter
 ierr = 0

 ! boundaries
 arg(3) = tab_para(1,3)*1.00000001d00
 y_min = arg(3)

 !++++++++++++++++++++++++++++++++++++
 call get_eos_sub_new(ipl,ierr,0,1,i_entr,eos_thermo) 
 !++++++++++++++++++++++++++++++++++++

 f_min = eos_thermo(5)

 arg(3) = tab_para(dim_idx(3),3)*0.99999999d00
 y_max = arg(3)

 !++++++++++++++++++++++++++++++++++++
 call get_eos_sub_new(ipl,ierr,0,1,i_entr,eos_thermo)
 !++++++++++++++++++++++++++++++++++++
 f_max = eos_thermo(5)

 if ((f_min*f_max > 0.d00).or.(y_min > y_max)) then
    ! no beta equilibrium
    if(y_min>y_max) then
       y = -2.d00
    else
       y = y_min ! return the value at the lowest Y entry of the table
       arg(3) = y
       !++++++++++++++++++++++++++++++++++++
       call get_eos_sub_new(ipl,ierr,0,1,i_entr,eos_thermo)
       !++++++++++++++++++++++++++++++++++++
    end if

   return
 else
  
#ifdef DEBUG
   timei = omp_get_wtime()
#endif

   ! determination of zero

   d_y = 1.d00
   do while (dabs(d_y) > 1.d-12)
     d_y = y
     d_f = f_max-f_min
     if (dabs(d_f) > 1.d-08) then
       y = (f_max*y_min-f_min*y_max)/d_f
     else
       y = 0.5d00*(y_max+y_min)
     end if
     arg(3) = y
     !++++++++++++++++++++++++++++++++++++
     call get_eos_sub_new(ipl,ierr,0,1,i_entr,eos_thermo)
     !++++++++++++++++++++++++++++++++++++
     f = eos_thermo(5)
     d_y = d_y-y
     if (f_min*f > 0.d00) then
       y_min = y
       f_min = f
     else
       y_max = y
       f_max = f
     end if
   end do
 endif

#ifdef DEBUG
 timef = omp_get_wtime()
 ! print '(2x,a,f12.5,a,i3,a)','TIME for read tables : ',&
 !   &    timef-timei,  '(s) with ', 1, ' threads'
#endif

 !++++++++++++++++++++++++++++++++++++
 call get_eos_sub_new(ipl,ierr,0,0,i_entr,eos_thermo)
 !++++++++++++++++++++++++++++++++++++
end SUBROUTINE get_eos_beta
!***********************************************************************
SUBROUTINE get_eos_sub_s2(ipl,ierr,inmp,ibeta,i_entr)
 ! Stefan Typel for the CompOSE core team, version 1.00, 2017/11/23
 use compose_internal, only : dim_qtyt,arg, arg2, tab_para, dim_idx
 implicit none
 integer,intent(in) :: inmp, ibeta, ipl(3)
 integer,intent(out) :: ierr
 integer :: i_entr,i1
 double precision :: s,t_min,t_max,s_min,s_max,ds, eos_thermo(dim_qtyt)

 if (i_entr == 0) then
   !+++++++++++++++++++++++++++++++++++++++++++
   call get_eos_sub_new(ipl,ierr,inmp,ibeta,i_entr,eos_thermo)
   !+++++++++++++++++++++++++++++++++++++++++++
 else
   ! entropy
   s = arg(1)
   !write(*,*) ' s =',s
   !write(*,*) ' t_min =', tab_para(1,1)
   !write(*,*) ' t_max =', tab_para(dim_idx(1),1)

   ! boundaries

   t_min = tab_para(1,1)*1.00000001d00
   arg(1) = t_min
   arg2(1) = arg(1)
   !+++++++++++++++++++++++++++++++++++++++++++
   call get_eos_sub_new(ipl,ierr,inmp,ibeta,i_entr,eos_thermo)
   !+++++++++++++++++++++++++++++++++++++++++++
   s_min = eos_thermo(2)
   if (s == s_min) return

   t_max = tab_para(dim_idx(1),1)*0.99999999d00
   arg(1) = t_max
   arg2(1) = arg(1)
   !+++++++++++++++++++++++++++++++++++++++++++
   call get_eos_sub_new(ipl,ierr,inmp,ibeta,i_entr,eos_thermo)
   !+++++++++++++++++++++++++++++++++++++++++++
   s_max = eos_thermo(2)
   if (s == s_max) return

   !write(*,*) ' s_min s_max',s_min,s_max
   if ((s < s_max).and.(s > s_min)) then
     s_max = s_max-s
     s_min = s_min-s
     ds = 1.d00
     i1 = 0
     do while (dabs(ds) > 1.d-10*s)
       i1 = i1+1
       if (i1 == 40) then
         i_entr = 2
         return
       end if

       if (dabs(s_max-s_min) < s) then
         arg(1) = (s_max*t_min-s_min*t_max)/(s_max-s_min)
       else
         arg(1) = 0.5*(t_min+t_max)
       end if
       arg2(1) = arg(1)
       !+++++++++++++++++++++++++++++++++++++++++++
       call get_eos_sub_new(ipl,ierr,inmp,ibeta,i_entr,eos_thermo)
       !+++++++++++++++++++++++++++++++++++++++++++
       ds = eos_thermo(2)-s
       !write(*,*) i1,arg(1),eos_thermo(2),ds
       if (ds > 0.d00) then
         t_max = arg(1)
         s_max = ds
       else
         if (ds < 0.d00) then
           t_min = arg(1)
           s_min = ds
         else
           return
         end if
       end if
     end do
   else
     i_entr = 2
   end if
 end if

 return
end SUBROUTINE get_eos_sub_s2
!***********************************************************************
subroutine get_eos_sub(ipl,ierr,inmp,ibeta,i_entr)
 ! Stefan Typel for the CompOSE core team, version 1.09, 2017/11/23
 USE compose_internal
 implicit none
 integer :: ierr,inmp,ibeta,i_entr,m(3),ip,ic,i1,i2,i3,ipl(3)
 double precision :: q(3)

 ! default values for interpolation rules
 do ip=1,3,1
   if ((ipl(ip) < 1).or.(ipl(ip) > 3)) ipl(ip) = 3
 end do

 ! grid parameters
 do ip = 1,3,1
   if (idx_ipl(ip) == 1) then
     !++++++++++++++++++++++++++++++++++
     call get_eos_grid_para(ip,m,q,ierr)
     !++++++++++++++++++++++++++++++++++
   else
     m(ip) = 1
     q(ip) = 0.d00
   end if
 end do

 ! test for existence of corner points
 ic = 1
 ! 3 dim
 if (all(idx_ipl == [1,1,1])) then
   do i1=m(1),m(1)+1,1
     do i2=m(2),m(2)+1,1
       do i3=m(3),m(3)+1,1
         ic = ic*idx_arg(i1,i2,i3,0)
       end do
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 53
   end if


   ! 2 dim
 else if (all(idx_ipl == [0,1,1])) then
   i1 = m(1)
   do i2=m(2),m(2)+1,1
     do i3=m(3),m(3)+1,1
       ic = ic*idx_arg(i1,i2,i3,0)
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 52
   end if

 else if (all(idx_ipl == [1,0,1])) then
   do i1=m(1),m(1)+1,1
     i2 = m(2)
     do i3=m(3),m(3)+1,1
       ic = ic*idx_arg(i1,i2,i3,0)
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 52
   end if

 else if (all(idx_ipl == [1,1,0])) then
   do i1=m(1),m(1)+1,1
     do i2=m(2),m(2)+1,1
       i3 = m(3)
       ic = ic*idx_arg(i1,i2,i3,0)
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 52
   end if


   ! 1 dim
 else if (all(idx_ipl == [0,0,1])) then
   i1 = m(1)
   i2 = m(2)
   do i3=m(3),m(3)+1,1
     ic = ic*idx_arg(i1,i2,i3,0)
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 51
   end if
 else if (all(idx_ipl == [0,1,0])) then
   i1 = m(1)
   do i2=m(2),m(2)+1,1
     i3 = m(3)
     ic = ic*idx_arg(i1,i2,i3,0)
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 51
   end if

 else if (all(idx_ipl == [1,0,0])) then
   do i1=m(1),m(1)+1,1
     i2 = m(2)
     i3 = m(3)
     ic = ic*idx_arg(i1,i2,i3,0)
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 51
   end if
 end if

 if (i_entr == 0) then
   call write_errors(ierr)
 end if

 if (ierr == 0) then
   !++++++++++++++++++++++++++++++++++++
   call eos_interpol(m,q,ipl,inmp,ibeta,eos_thermo)
   !++++++++++++++++++++++++++++++++++++
 end if


end subroutine get_eos_sub

!***********************************************************************
subroutine get_eos_sub_new(ipl,ierr,inmp,ibeta,i_entr, thermo_out)
 ! Stefan Typel for the CompOSE core team, version 1.09, 2017/11/23
 use compose_internal, only : idx_arg, dim_qtyt, idx_ipl,error_msg, dim_err
 implicit none

 double precision,intent(out) :: thermo_out(dim_qtyt)
 integer,intent(out) :: ierr
 integer,intent(in) :: inmp,ibeta,i_entr
 integer :: ip,ic,i1,i2,i3,ipl(3),m(3)
 double precision :: q(3)

 ! default values for interpolation rules
 do ip=1,3,1
   if ((ipl(ip) < 1).or.(ipl(ip) > 3)) ipl(ip) = 3
 end do
 ierr = 0

 ! grid parameters
 do ip = 1,3,1
   if (idx_ipl(ip) == 1) then
     !++++++++++++++++++++++++++++++++++
     call get_eos_grid_para(ip,m,q,ierr)
     !++++++++++++++++++++++++++++++++++
   else
     m(ip) = 1
     q(ip) = 0.d00
   end if
 end do

 ! test for existence of corner points
 ic = 1
 ! 3 dim
 if (all(idx_ipl == [1,1,1])) then
   do i1=m(1),m(1)+1,1
     do i2=m(2),m(2)+1,1
       do i3=m(3),m(3)+1,1
         ic = ic*idx_arg(i1,i2,i3,0)
       end do
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 53
   end if


   ! 2 dim
 else if (all(idx_ipl == [0,1,1])) then
   i1 = m(1)
   do i2=m(2),m(2)+1,1
     do i3=m(3),m(3)+1,1
       ic = ic*idx_arg(i1,i2,i3,0)
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 52
   end if

 else if (all(idx_ipl == [1,0,1])) then
   do i1=m(1),m(1)+1,1
     i2 = m(2)
     do i3=m(3),m(3)+1,1
       ic = ic*idx_arg(i1,i2,i3,0)
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 52
   end if

 else if (all(idx_ipl == [1,1,0])) then
   do i1=m(1),m(1)+1,1
     do i2=m(2),m(2)+1,1
       i3 = m(3)
       ic = ic*idx_arg(i1,i2,i3,0)
     end do
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 52
   end if


   ! 1 dim
 else if (all(idx_ipl == [0,0,1])) then
   i1 = m(1)
   i2 = m(2)
   do i3=m(3),m(3)+1,1
     ic = ic*idx_arg(i1,i2,i3,0)
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 51
   end if
 else if (all(idx_ipl == [0,1,0])) then
   i1 = m(1)
   do i2=m(2),m(2)+1,1
     i3 = m(3)
     ic = ic*idx_arg(i1,i2,i3,0)
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 51
   end if

 else if (all(idx_ipl == [1,0,0])) then
   do i1=m(1),m(1)+1,1
     i2 = m(2)
     i3 = m(3)
     ic = ic*idx_arg(i1,i2,i3,0)
   end do
   if (ic /= 1) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 51
   end if
 end if

 if (i_entr == 0) then
   call write_errors(ierr)
 end if

 if (ierr == 0) then
   !++++++++++++++++++++++++++++++++++++
   call eos_interpol(m,q,ipl,inmp,ibeta,thermo_out)
   !++++++++++++++++++++++++++++++++++++
 end if


end subroutine get_eos_sub_new
!***********************************************************************
SUBROUTINE get_eos_grid_para(ip,m,q,ierr)
! Stefan Typel for the CompOSE core team, version 1.03, 2016/06/20
USE compose_internal
implicit none
integer :: ierr,ip,ia,im,m(3)
double precision q(3)

im = 1
do ia=1,dim_idx(ip),1
   if (arg(ip) >= (tab_para(ia,ip))) im = ia
end do
if (im > (dim_idx(ip)-1)) im = dim_idx(ip)-1
m(ip) = im
q(ip) = (arg(ip)-tab_para(im,ip))/(tab_para(im+1,ip)-tab_para(im,ip))
if ((q(ip) < -1.d-06).or.((q(ip)-1.d00) > 1.d-06)) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 40+ip
end if

return
end SUBROUTINE get_eos_grid_para
!***********************************************************************
subroutine eos_interpol(m,q,ipl,inmp,ibeta,eosthermo_out)
! Stefan Typel for the CompOSE core team, version 1.05, 2016/06/30
USE compose_internal
implicit none
double precision,intent(out) :: eosthermo_out(dim_qtyt)
integer :: m(3),ipl(3),inmp,ibeta,dim_ipl
double precision :: q(3)

if (irpl > 0) then
   write(*,*) 'case with B field not considered completely'
end if

dim_ipl = nall_max
if (np_max > dim_ipl) dim_ipl = np_max
if (nq_max > dim_ipl) dim_ipl = nq_max
if (nm_max > dim_ipl) dim_ipl = nm_max

select case (eos_dim)
case(3)
  !+++++++++++++++++++++++++++++++++++++++++++++++
  call eos_interpol_d3(m,q,ipl,inmp,ibeta,dim_ipl,eosthermo_out)
case(2)
  !+++++++++++++++++++++++++++++++++++++++++++++++
  call eos_interpol_d2(m,q,ipl,inmp,ibeta,dim_ipl,eosthermo_out)
case(1)
  !+++++++++++++++++++++++++++++++++++++++++++++++
  call eos_interpol_d1(m,q,ipl,inmp,ibeta,dim_ipl,eosthermo_out)

end select

end SUBROUTINE eos_interpol
!***********************************************************************
subroutine eos_interpol_d1(m,q,ipl,inmp,ibeta,dim_ipl,eosthermo_out)
! Stefan Typel for the CompOSE core team, version 1.07, 2017/05/23
USE compose_internal
implicit none
 double precision,intent(out) :: eosthermo_out(dim_qtyt)
integer, parameter :: dim_igp=100
integer :: m(3),inmp,ibeta,i1,i2,i3,i1p,i2p,i3p,iq,iq2,is,dim_ipl,&
     igp,ngp,igp3,ngp3,igp_r,igp3_r,j1,j2,j3,alloc_status,ierr,&
     idx1(dim_igp),idx2(dim_igp),idx3(10),ir(-4:5),irp,ipl(3),ik(3)
double precision :: q(3),vec(-4:5),vec3(-4:5,3),tmp,dh(0:2)
double precision, allocatable :: mat(:,:,:),mat3(:,:,:,:),&
     vp_compo(:),vq_compo(:,:),v_micro(:)

ierr = 0

allocate(mat(1:dim_ipl,-4:5,-4:5),stat=alloc_status)
if (alloc_status /= 0) then
   write(*,*) 'alloc_status = ',alloc_status
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 140
end if
mat = 0.d0

allocate(mat3(dim_ipl,-4:5,-4:5,3),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 141
end if
mat3 = 0.d0

allocate(vp_compo(np_max),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 142
end if

allocate(vq_compo(nq_max,3),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 143
end if

allocate(v_micro(nm_max),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 144
end if

if (dim_idx(1) > 1) then
  ik = [2,3,1]
else if (dim_idx(2) > 1) then
  ik = [1,3,2]
else if (dim_idx(3) > 1) then
  ik = [1,2,3]
else
  ierr = ierr+1
  if (ierr < dim_err) error_msg(ierr) = 164
end if

! list of grid points and reference point in indices 1 and 2
igp = 0
igp_r = 0
i1 = 0
i1p = i1+m(ik(1))
if ((i1p >= 1).and.(i1p <= dim_idx(ik(1)))) then
   i2 = 0
   i2p = i2+m(ik(2))
   if ((i2p >= 1).and.(i2p <= dim_idx(ik(2)))) then
      igp = igp+1
      idx1(igp) = i1
      idx2(igp) = i2
      igp_r = igp
   else
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 165
   end if
end if
ngp = igp

call write_errors(ierr)

! list of grid points and reference point in index 3
if (idx_ipl(ik(3)) == 1) then
   igp3 = 0
   igp3_r = 0
   do i3=-4,5,1 ! Y_q
      i3p = i3+m(ik(3))
      if ((i3p >= 1).and.(i3p <= dim_idx(ik(3)))) then
         igp3 = igp3+1
         idx3(igp3) = i3
         if (q(ik(3)) <= 0.5d00) then
            if (i3 == 0) igp3_r = igp3
         else
            if (i3 == 1) igp3_r = igp3
         end if
      end if
   end do
   ngp3 = igp3
end if

! list of grid points and reference point in index 3
if (idx_ipl(ik(3)) == 1) then
   igp3 = 0
   igp3_r = 0
   do i3=-4,5,1 ! Y_q
      i3p = i3+m(ik(3))
      if ((i3p >= 1).and.(i3p <= dim_idx(ik(3)))) then
         igp3 = igp3+1
         idx3(igp3) = i3
         if (q(ik(3)) <= 0.5d00) then
            if (i3 == 0) igp3_r = igp3
         else
            if (i3 == 1) igp3_r = igp3
         end if
      end if
   end do
   ngp3 = igp3
end if

! thermodynamic quantities

if (ibeta == 1) then
   ! for determination of beta-equilibrium
   if (nall_max > 0) idx_thermo(1:nall_max) = 0
   idx_thermo(5) = 1
else
   do iq=1,nall_max,1
      idx_thermo(iq) = iq
   end do
end if

!2017/05/23
eos_df(1:10) = 0.d00
! interpolation in index 3
!2016/10/28
if (nall_max > 0) v_thermo(1:nall_max,0:4) = 0.d00
do iq=1,nall_max,1
   if (idx_thermo(iq) > 0) then
      do igp=1,ngp,1
         i1p = idx1(igp)+m(ik(1))
         i2p = idx2(igp)+m(ik(2))
         do i3=-4,5,1
            vec(i3) = 0.d00
         end do
         do igp3=1,ngp3,1
            i3p = idx3(igp3)+m(ik(3))
            if (dim_idx(1) > 1) then
              j1 = i3p
              j2 = i1p
              j3 = i2p
            else if (dim_idx(2) > 1) then
              j1 = i1p
              j2 = i3p
              j3 = i2p
            else if (dim_idx(3) > 1) then
              j1 = i1p
              j2 = i2p
              j3 = i3p
            end if
            vec(idx3(igp3)) = tab_thermo(j1,j2,j3,iq)
            irp = idx_arg(j1,j2,j3,ik(3))
            ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
         end do
         call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
         v_thermo(iq,0) = dh(0)
         v_thermo(iq,1) = dh(1)
         v_thermo(iq,2) = dh(2)
!2017/05/23
! eos_df(1): F [MeV]
! eos_df(2): dF/dT []
! eos_df(3): d^2F/(dT^2) [MeV^-1]
! eos_df(4): d^2F/(dT dn_b) [fm^3]
! eos_df(5): d^2F/(dT dY_q) []
! eos_df(6): dF/dn_b [MeV fm^3]
! eos_df(7): d^2F/(dn_b^2) [MeV fm^6]
! eos_df(8): d^2F/(dn_b dY_q) [MeV fm^3]
! eos_df(9): dF/dY_q [MeV]
! eos_df(10): d^2dF/(dY_q^2) [MeV]
         if (iq == 6) then
            eos_df(1) = (dh(0)+1.d00)*m_n
            if (dim_idx(1) > 1) then
               eos_df(2) = dh(1)*m_n
               eos_df(3) = dh(2)*m_n
            end if
            if (dim_idx(2) > 1) then
               eos_df(6) = dh(1)*m_n
               eos_df(7) = dh(2)*m_n
            end if
            if (dim_idx(3) > 1) then
               eos_df(9) = dh(1)*m_n
               eos_df(10) = dh(2)*m_n
            end if
         end if
      end do
   end if

end do

! mu_l lepton number chemical potential
eosthermo_out(5) = v_thermo(5,0)*m_n

if (ibeta /= 1) then

   ! p pressure
   eosthermo_out(1) = v_thermo(1,0)*arg2(2)
   ! S entropy ber baryon
   eosthermo_out(2) = v_thermo(2,0)
   ! mu_b-m_n baryon number chemical potential
   eosthermo_out(3) = v_thermo(3,0)*m_n
   ! mu_q charge chemical potential
   eosthermo_out(4) = v_thermo(4,0)*m_n

   ! F/m_n-1 free energy per baryon
   eosthermo_out(6) = v_thermo(6,0)
   ! E/m_n-1 internal energy per baryon
   eosthermo_out(7) = v_thermo(7,0)
   !eosthermo_out(7) = v_thermo(6,0)+x(1)*eosthermo_out(2)/m_n

   if (inmp == 1) then
      eos_thermo_add(1) = 0.d00
      eos_thermo_add(2) = 0.d00
      return
   endif

   tmp = v_thermo(1,0)/m_n
   ! H/m_n-1  enthalpy per baryon
   eosthermo_out(8) = eosthermo_out(7)+tmp
   ! G/m_n-1 free enthalpy per baryon
   eosthermo_out(9) = eosthermo_out(6)+tmp

   ! F = f/n
   eosthermo_out(20) = (eosthermo_out(6)+1.d00)*m_n
   ! E = e/n
   eosthermo_out(21) = (eosthermo_out(7)+1.d00)*m_n
   ! H = h/n
   eosthermo_out(22) = (eosthermo_out(8)+1.d00)*m_n
   ! G = g/n
   eosthermo_out(23) = (eosthermo_out(9)+1.d00)*m_n
   ! epsilon = e
   eosthermo_out(24) = eosthermo_out(21)*arg2(2)

   ! error quantities

   ! delta (f/n) = f/n+p/n-(mu_b+y_q*(mu_q+mu_e), mu_e = mu_l-mu_q

   if (incl_l == 1) then
      ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_l)
      eos_thermo_err(1) = eosthermo_out(6)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(5))
   else
      ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_q)
      eos_thermo_err(1) = eosthermo_out(6)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(4))
   end if

   ! (delta f/n)/(f/n)
   if (eosthermo_out(6) /= -1.d00) then
      eos_thermo_err(2) = eos_thermo_err(1)/((eosthermo_out(6)+1.d00)*m_n)
   else
      eos_thermo_err(2) = 0.d00
   end if

   if (incl_l == 0) then
      ! delta (e/n) = f/n+p/n-(mu_b+y_q*mu_q)-Ts
      eos_thermo_err(3) = eosthermo_out(7)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(4))&
           -arg2(1)*eosthermo_out(2)
   else
      ! delta (e/n) = e/n+p/n-(mu_b+y_q*mu_l)-Ts
      eos_thermo_err(3) = eosthermo_out(7)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(5))&
           -arg2(1)*eosthermo_out(2)
   end if

   ! (delta e/n)/(e/n)
   if (eosthermo_out(7) /= -1.d00) then
      eos_thermo_err(4) = eos_thermo_err(3)/((eosthermo_out(7)+1.d00)*m_n)
   else
      eos_thermo_err(4) = 0.d00
   end if

   ! delta (p/n) = p/n-n*d(f/n)/dn
   eos_thermo_err(5) = 0.d00
   eos_thermo_err(6) = 0.d00

   ! delta (s/n) = s/n+d(f/n)/dt
   eos_thermo_err(7) = 0.d00
   eos_thermo_err(8) = 0.d00

   ! quantities depending on derivatives

   eosthermo_out(10:19) = 0.d00

   if (irpl == 0) then
      if (ik(3) == 1) then
         ! beta_v
         eosthermo_out(17) = v_thermo(1,1)*arg2(2)
         ! c_v
         eosthermo_out(13) = v_thermo(2,1)*arg2(1)
         ! dp/de|n
         eosthermo_out(11) = arg2(1)*v_thermo(2,1)
         if (eosthermo_out(11) /= 0.d00) then
            eosthermo_out(11) = arg2(2)*v_thermo(1,1)/eosthermo_out(11)
         else
            eosthermo_out(11) = 0.d00
         end if
      end if
      if (ik(3) == 2) then
         ! kappa_t
         eosthermo_out(18) = arg2(2)*(arg2(2)*v_thermo(1,1)+v_thermo(1,0))
         if (eosthermo_out(18) /= 0.d00) then
            eosthermo_out(18) = 1.d00/eosthermo_out(18)
         else
            eosthermo_out(18) = 0.d00
         end if
         if (arg2(1) == 0.d00) then
            ! gamma
            eosthermo_out(15) = 1.d00
            ! kappa_s
            if (eosthermo_out(15) /= 0.d00) then
               eosthermo_out(19) = eosthermo_out(18)/eosthermo_out(15)
            else
               eosthermo_out(19) = 0.d00
            end if
            ! cs^2
            eosthermo_out(12) = eosthermo_out(19)*(eosthermo_out(8)+1.d00)*m_n*arg2(2)
            if (eosthermo_out(12) /= 0.d00) then
               eosthermo_out(12) = 1.d00/eosthermo_out(12)
            end if
         end if
      end if
   end if

   ! additional quantities
   do iq=1,nadd_max,1
      eos_thermo_add(iq) = v_thermo(iq+dim_reg,0)
   end do

   ! compositional quantities

   if (np_max > 0) then
      idx_compo_p(1:np_max) = 0
      eos_compo_p(1:np_max) = 0.d00
   end if
   if (nq_max > 0) then
      idx_compo_q(1:nq_max) = 0
      eos_compo_q(1:nq_max,1:4) = 0.d00
   end if

   if (idx_ex(2) == 1) then

      ! pairs

      ! interpolation in index 3
      do iq=1,np_max,1
         do igp=1,ngp,1
            i1p = idx1(igp)+m(ik(1))
            i2p = idx2(igp)+m(ik(2))
            do i3=-4,5,1
               vec(i3) = 0.d00
            end do
            do igp3=1,ngp3,1
              i3p = idx3(igp3)+m(ik(3))
              if (dim_idx(1) > 1) then
                j1 = i3p
                j2 = i1p
                j3 = i2p
              else if (dim_idx(2) > 1) then
                j1 = i1p
                j2 = i3p
                j3 = i2p
              else if (dim_idx(3) > 1) then
                j1 = i1p
                j2 = i2p
                j3 = i3p
              end if
               vec(idx3(igp3)) = 0.d00
               do iq2=1,np_max,1
                  if (idxp_compo(j1,j2,j3,iq2) == idx_p(iq)) then
                     vec(idx3(igp3)) = tabp_compo(j1,j2,j3,iq2)
                  end if
               end do
               irp = idx_arg(j1,j2,j3,ik(3))
               ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
            end do
            call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
            vp_compo(iq) = dh(0)
         end do
      end do

      do iq=1,np_max,1
         !         if (vp_compo(iq) > 0.d00) then
         eos_compo_p(iq) = vp_compo(iq)
         !         end if
         idx_compo_p(iq) = idx_p(iq)
      end do

      ! quadruples

      ! interpolation in index 3
      do iq=1,nq_max,1
         do igp=1,ngp,1
            i1p = idx1(igp)+m(ik(1))
            i2p = idx2(igp)+m(ik(2))
            vec(-4:5) = 0.d00
            vec3(-4:5,1:3) = 0.d00
            do igp3=1,ngp3,1
              i3p = idx3(igp3)+m(ik(3))
              if (dim_idx(1) > 1) then
                j1 = i3p
                j2 = i1p
                j3 = i2p
              else if (dim_idx(2) > 1) then
                j1 = i1p
                j2 = i3p
                j3 = i2p
              else if (dim_idx(3) > 1) then
                j1 = i1p
                j2 = i2p
                j3 = i3p
              end if
              vec3(idx3(igp3),1:3) = 0.d00
              do iq2=1,nq_max,1
                if (idxq_compo(j1,j2,j3,iq2) == idx_q(iq)) then
                  do is=1,3,1
                    vec3(idx3(igp3),is) = tabq_compo(j1,j2,j3,3*(iq2-1)+is)
                  end do
                end if
              end do
              irp = idx_arg(j1,j2,j3,ik(3))
              ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
            end do
            do is=1,3,1
               do i3=-4,5,1
                  vec(i3) = vec3(i3,is)
               end do
               call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
               vq_compo(iq,is) = dh(0)
            end do
         end do
      end do

      do iq=1,nq_max,1
         if ((vq_compo(iq,3) > 0.d00).and.(vq_compo(iq,2) > 0.d00).and.&
              (vq_compo(iq,1) > 0.d00)) then
            eos_compo_q(iq,1) = vq_compo(iq,3)
            eos_compo_q(iq,2) = vq_compo(iq,1)
            eos_compo_q(iq,3) = vq_compo(iq,2)
            eos_compo_q(iq,4) = vq_compo(iq,1)-vq_compo(iq,2)
         else
            do is=1,4,1
               eos_compo_q(iq,is) = 0.d00
            end do
         end if
         idx_compo_q(iq) = idx_q(iq)
      end do
   end if

   ! microscopic quantities

   if (nm_max > 0) then
      idx_micro(1:nm_max) = 0
      eos_micro(1:nm_max) = 0.d00
   end if

   if (idx_ex(3) == 1) then

      ! interpolation in index 3
      do iq=1,nm_max,1
         do igp=1,ngp,1
            i1p = idx1(igp)+m(ik(1))
            i2p = idx2(igp)+m(ik(2))
            do i3=-4,5,1
               vec(i3) = 0.d00
            end do
            do igp3=1,ngp3,1
               i3p = idx3(igp3)+m(ik(3))
               if (dim_idx(1) > 1) then
                j1 = i3p
                j2 = i1p
                j3 = i2p
              else if (dim_idx(2) > 1) then
                j1 = i1p
                j2 = i3p
                j3 = i2p
              else if (dim_idx(3) > 1) then
                j1 = i1p
                j2 = i2p
                j3 = i3p
              end if

               vec(idx3(igp3)) = 0.d00
               do iq2=1,nm_max,1
                  if (idx_mic(j1,j2,j3,iq2) == idx_m(iq)) then
                     vec(idx3(igp3)) = tab_mic(j1,j2,j3,iq2)
                  end if
               end do
               irp = idx_arg(j1,j2,j3,ik(3))
               ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
            end do
            call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
            v_micro(iq) = dh(0)
         end do
      end do

      do iq=1,nm_max,1
         eos_micro(iq) = v_micro(iq)
         idx_micro(iq) = idx_m(iq)
      end do
   end if

end if

if (allocated(mat)) deallocate(mat)
if (allocated(mat3)) deallocate(mat3)
if (allocated(vp_compo)) deallocate(vp_compo)
if (allocated(vq_compo)) deallocate(vq_compo)
if (allocated(v_micro)) deallocate(v_micro)

end SUBROUTINE eos_interpol_d1
!***********************************************************************
subroutine eos_interpol_d3(m,q,ipl,inmp,ibeta,dim_ipl, eosthermo_out)
 ! Stefan Typel for the CompOSE core team, version 1.17, 2017/05/22
 use compose_internal
 use omp_lib
 implicit none
 logical, save :: first = .true.
 integer,intent(in) :: m(3),inmp,ibeta,dim_ipl
 double precision,intent(in) :: q(3)
 double precision,intent(out) :: eosthermo_out(dim_qtyt)
 integer, parameter :: dim_igp=100
 integer :: i1,i2,i3,i1p,i2p,i3p,iq,iq2,is,&
   igp,ngp,igp3,ngp3,igp_r,igp3_r,j1,j2,j3,alloc_status,ierr,&
   idx1(dim_igp),idx2(dim_igp),idx3(10),ir(-4:5),irp,ipl(3),ik(3)
 double precision :: vec(-4:5),vec3(-4:5,3),&
   qx,qy,dg(0:2,0:2),tmp,dh(0:2)
 double precision, allocatable,save :: mat(:,:,:),mat3(:,:,:,:),&
   vp_compo(:),vq_compo(:,:),v_micro(:)
 double precision :: mat2(-4:5,-4:5,0:2)
 integer :: nthreads
#ifdef DEBUG
 double precision :: timei,timef
#endif


 ierr = 0


 if(first) then

   allocate(mat(1:dim_ipl,-4:5,-4:5),stat=alloc_status)
   if (alloc_status /= 0) then
     write(*,*) 'alloc_status = ',alloc_status
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 140
   end if
   mat = 0.d0

   allocate(mat3(dim_ipl,-4:5,-4:5,3),stat=alloc_status)
   if (alloc_status /= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 141
   end if
   mat3 = 0.d0

   !2017/05/22
   allocate(vp_compo(np_max),stat=alloc_status)
   if (alloc_status /= 0) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 142
   end if

   allocate(vq_compo(nq_max,3),stat=alloc_status)
   if (alloc_status /= 0) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 143
   end if

   allocate(v_micro(nm_max),stat=alloc_status)
   if (alloc_status /= 0) then
     ierr = ierr+1
     if (ierr < dim_err) error_msg(ierr) = 144
   end if

   call write_errors(ierr)

   first = .false.
 end if

 ! initialisation of mat2
 mat2 = 0.d0

 ! interpolation for three dimensional table
 ! standard choice: two-dimensional interpolation in ik(1) and ik(2)
 !                  one-dimensional interpolation in ik(3)

 ! irpl == 0 or 3
 select case(irpl)
 case(0)
   ik = [ (i1,  i1 = 1,3, 1 )]
 case(1)
   ik = [ 2, 3, 1]
 case(2)
   ik = [ 1, 3, 2]
 case default
   ik = [ (i1,  i1 = 1,3, 1 )]
 end select

 ! list of grid points and reference point in indices 1 and 2
 igp = 0
 igp_r = 0
 do i1=-4,5,1 ! first index
   i1p = i1+m(ik(1))
   if ((i1p >= 1).and.(i1p <= dim_idx(ik(1)))) then
     do i2=-4,5,1 ! second index
       i2p = i2+m(ik(2))
       if ((i2p >= 1).and.(i2p <= dim_idx(ik(2)))) then
         igp = igp+1
         idx1(igp) = i1
         idx2(igp) = i2
         if (q(ik(1)) <= 0.5d00) then
           if (q(ik(2)) <= 0.5d00) then
             if ((i1 == 0).and.(i2 == 0)) igp_r = igp
           else
             if ((i1 == 0).and.(i2 == 1)) igp_r = igp
           end if
         else
           if (q(ik(2)) <= 0.5d00) then
             if ((i1 == 1).and.(i2 == 0)) igp_r = igp
           else
             if ((i1 == 1).and.(i2 == 1)) igp_r = igp
           end if
         end if
       end if
     end do
   end if
 end do
 ngp = igp
 ! list of grid points and reference point in index 3
 if (idx_ipl(ik(3)) == 1) then
   igp3 = 0
   igp3_r = 0
   do i3=-4,5,1 ! Y_q
     i3p = i3+m(ik(3))
     if ((i3p >= 1).and.(i3p <= dim_idx(ik(3)))) then
       igp3 = igp3+1
       idx3(igp3) = i3
       if (q(ik(3)) <= 0.5d00) then
         if (i3 == 0) igp3_r = igp3
       else
         if (i3 == 1) igp3_r = igp3
       end if
     end if
   end do
   ngp3 = igp3
 end if


 ! for interpolation in first and second index
 call get_idx_arg2(m,ik)
 call get_diff_rules2(m,ipl,ik)

 ! thermodynamic quantities

 if (ibeta == 1) then
   ! for determination of beta-equilibrium

   if (nall_max > 0) idx_thermo(1:nall_max) = 0
   idx_thermo(5) = 1
 else
   do iq=1,nall_max,1
     idx_thermo(iq) = iq
   end do
 end if

 ! interpolation in index 3
 !2016/10/28
 if (nall_max > 0) mat(1:nall_max,-4:5,-4:5) = 0.d00

 qx = q(ik(1))
 qy = q(ik(2))

! #ifdef DEBUG
!  timei = omp_get_wtime()
! #endif

 !! !$OMP PARALLEL &
 !! !$OMP DEFAULT(SHARED)

 !! !$OMP DO &
 !! !$OMP PRIVATE(iq,igp,i1p,i2p,vec,igp3,i3p,j1,j2,j3,irp,dh,dg)
 do iq=1,nall_max,1
   if (idx_thermo(iq) > 0) then
     do igp=1,ngp,1
       i1p = idx1(igp)+m(ik(1))
       i2p = idx2(igp)+m(ik(2))
       vec(-4:5) = 0.d00
       do igp3=1,ngp3,1
         i3p = idx3(igp3)+m(ik(3))
         select case(irpl)
         case(1)
           j1 = i3p
           j2 = i1p
           j3 = i2p
         case(2)
           j1 = i1p
           j2 = i3p
           j3 = i2p
         case default
           j1 = i1p
           j2 = i2p
           j3 = i3p
         end select

         vec(idx3(igp3)) = tab_thermo(j1,j2,j3,iq)
         irp = idx_arg(j1,j2,j3,ik(3))
         if  (idx3(igp3) > 0) then
           ir(idx3(igp3)) = ipl_rule(1,ipl(ik(3)),irp)
         else
           ir(idx3(igp3)) = ipl_rule(0,ipl(ik(3)),irp)
           ! ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
         endif
       end do
       call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
       mat(iq,idx1(igp),idx2(igp)) = dh(0)
       !2017/05/22
       if (iq == 6) then
          mat2(idx1(igp),idx2(igp),0:2) = dh(0:2)
       end if
     end do

     ! interpolation in index 1 and 2
     call make_interp_xy(ipl,ik,qx,qy,dg,1,mat(iq,-4:5,-4:5))
     v_thermo(iq,0) = dg(0,0)
     v_thermo(iq,1) = dg(0,1)
     v_thermo(iq,2) = dg(1,0)
     v_thermo(iq,3) = dg(0,2)
     v_thermo(iq,4) = dg(2,0)
   else
     v_thermo(iq,0:4) = 0.d00
   end if
 end do
 !! !$OMP END DO


 !2017/05/22
 !! !$OMP WORKSHARE
 eos_df(1:10) = 0.d00
 !! !$OMP END WORKSHARE

 if (irpl == 0) then
   if (idx_thermo(6) > 0) then
     !! !$OMP DO PRIVATE(is,dg)
     do is=0,2
        call make_interp_xy(ipl,ik,qx,qy,dg,1,mat2(-4:5,-4:5,is))
       ! call get_derivatives(ipl,ik)
       ! call get_coefficients()
       ! call get_interpol_xy(qx,qy,dg,1)
       ! eos_df(1): F [MeV]
       ! eos_df(2): dF/dT [MeV]
       ! eos_df(3): d^2F/(dT^2) []
       ! eos_df(4): d^2F/(dT dn_b) [fm^3]
       ! eos_df(5): d^2F/(dT dY_q) []
       ! eos_df(6): dF/dn_b [MeV fm^3]
       ! eos_df(7): d^2F/(dn_b^2) [MeV fm^6]
       ! eos_df(8): d^2F/(dn_b dY_q) [MeV fm^3]
       ! eos_df(9): dF/dY_q [MeV]
       ! eos_df(10): d^2dF/(dY_q^2) [MeV]
       select case (is)
       case(0)
         eos_df(1) = (dg(0,0)+1.d00)*m_n
         eos_df(2) = dg(1,0)*m_n
         eos_df(3) = dg(2,0)*m_n
         eos_df(4) = dg(1,1)*m_n
         eos_df(6) = dg(0,1)*m_n
         eos_df(7) = dg(0,2)*m_n
       case(1)
         eos_df(5) = dg(1,0)*m_n
         eos_df(8) = dg(0,1)*m_n
         eos_df(9) = dg(0,0)*m_n
       case(2)
         eos_df(10) = dg(0,0)*m_n
       end select
     end do
     !! !$OMP END DO
   end if
 end if
 !! !$OMP SINGLE
 ! nthreads = omp_get_num_threads()

 ! mu_l lepton number chemical potential
 eosthermo_out(5) = v_thermo(5,0)*m_n
 !! !$OMP END SINGLE


 if (ibeta /= 1) then

   !! !$OMP SINGLE
   ! p pressure
   eosthermo_out(1) = v_thermo(1,0)*arg2(2)
   ! S entropy ber baryon
   eosthermo_out(2) = v_thermo(2,0)
   ! mu_b-m_n baryon number chemical potential
   eosthermo_out(3) = v_thermo(3,0)*m_n
   ! mu_q charge chemical potential
   eosthermo_out(4) = v_thermo(4,0)*m_n

   ! F/m_n-1 free energy per baryon
   eosthermo_out(6) = v_thermo(6,0)
   ! E/m_n-1 internal energy per baryon
   eosthermo_out(7) = v_thermo(7,0)
   !eosthermo_out(7) = v_thermo(6,0)+x(1)*eosthermo_out(2)/m_n

   tmp = v_thermo(1,0)/m_n
   ! H/m_n-1  enthalpy per baryon
   eosthermo_out(8) = eosthermo_out(7)+tmp
   ! G/m_n-1 free enthalpy per baryon
   eosthermo_out(9) = eosthermo_out(6)+tmp

   ! F = f/n
   eosthermo_out(20) = (eosthermo_out(6)+1.d00)*m_n
   ! E = e/n
   eosthermo_out(21) = (eosthermo_out(7)+1.d00)*m_n
   ! H = h/n
   eosthermo_out(22) = (eosthermo_out(8)+1.d00)*m_n
   ! G = g/n
   eosthermo_out(23) = (eosthermo_out(9)+1.d00)*m_n
   ! epsilon = e
   eosthermo_out(24) = eosthermo_out(21)*arg2(2)
   !! !$OMP END SINGLE

   if (inmp == 1) then
     !! !$OMP SINGLE
     eos_thermo_add(1) = v_thermo(1,1)
     eos_thermo_add(2) = v_thermo(1,3)
     !! !$OMP END SINGLE
   else
     !! !$OMP SINGLE

     ! error quantities

     ! delta (f/n) = f/n+p/n-(mu_b+y_q*(mu_q+mu_e), mu_e = mu_l-mu_q

     if (incl_l == 1) then
       ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_l)
       eos_thermo_err(1) = eosthermo_out(6)*m_n+v_thermo(1,0)&
         -(eosthermo_out(3)+arg2(3)*eosthermo_out(5))
     else
       ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_q)
       eos_thermo_err(1) = eosthermo_out(6)*m_n+v_thermo(1,0)&
         -(eosthermo_out(3)+arg2(3)*eosthermo_out(4))
     end if

     ! (delta f/n)/(f/n)
     if (eosthermo_out(6) /= -1.d00) then
       eos_thermo_err(2) = eos_thermo_err(1)/((eosthermo_out(6)+1.d00)*m_n)
     else
       eos_thermo_err(2) = 0.d00
     end if

     if (incl_l == 0) then
       ! delta (e/n) = f/n+p/n-(mu_b+y_q*mu_q)-Ts
       eos_thermo_err(3) = eosthermo_out(7)*m_n+v_thermo(1,0)&
         -(eosthermo_out(3)+arg2(3)*eosthermo_out(4))&
         -arg2(1)*eosthermo_out(2)
     else
       ! delta (e/n) = e/n+p/n-(mu_b+y_q*mu_l)-Ts
       eos_thermo_err(3) = eosthermo_out(7)*m_n+v_thermo(1,0)&
         -(eosthermo_out(3)+arg2(3)*eosthermo_out(5))&
         -arg2(1)*eosthermo_out(2)
     end if

     ! (delta e/n)/(e/n)
     if (eosthermo_out(7) /= -1.d00) then
       eos_thermo_err(4) = eos_thermo_err(3)/((eosthermo_out(7)+1.d00)*m_n)
     else
       eos_thermo_err(4) = 0.d00
     end if

     ! delta (p/n) = p/n-n*d(f/n)/dn
     eos_thermo_err(5) = v_thermo(1,0)-arg2(2)*v_thermo(6,1)*m_n
     if (v_thermo(1,0) /= 0.d00) then
       eos_thermo_err(6) = eos_thermo_err(5)/v_thermo(1,0)
     else
       eos_thermo_err(6) = 0.d00
     end if

     ! delta (s/n) = s/n+d(f/n)/dt
     eos_thermo_err(7) = eosthermo_out(2)+v_thermo(6,2)*m_n
     if (eosthermo_out(2) /= 0.d00) then
       eos_thermo_err(8) = eos_thermo_err(7)/eosthermo_out(2)
     else
       eos_thermo_err(8) = 0.d00
     end if

     ! quantities depending on derivatives

     eosthermo_out(10:19) = 0.d00
     if ((irpl == 0).or.(irpl == 3)) then
       ! beta_v
       eosthermo_out(17) = v_thermo(1,2)*arg2(2)
       ! kappa_t
       eosthermo_out(18) = arg2(2)*(arg2(2)*v_thermo(1,1)+v_thermo(1,0))
       if (eosthermo_out(18) /= 0.d00) then
         eosthermo_out(18) = 1.d00/eosthermo_out(18)
       else
         eosthermo_out(18) = 0.d00
       end if
       ! alpha_p
       eosthermo_out(16) = eosthermo_out(17)*eosthermo_out(18)
       ! c_v
       eosthermo_out(13) = v_thermo(2,2)*arg2(1)
       ! c_p
       eosthermo_out(14) = eosthermo_out(13)+eosthermo_out(16)*eosthermo_out(17)*arg2(1)/arg2(2)
       ! gamma
       if (eosthermo_out(13) /= 0.d00) then
         eosthermo_out(15) = eosthermo_out(14)/eosthermo_out(13)
       else
         eosthermo_out(15) = 1.d00
       end if
       ! kappa_s
       if (eosthermo_out(15) /= 0.d00) then
         eosthermo_out(19) = eosthermo_out(18)/eosthermo_out(15)
       else
         eosthermo_out(19) = 0.d00
       end if
       ! dp/de|n
       eosthermo_out(11) = arg2(1)*v_thermo(2,2)
       if (eosthermo_out(11) /= 0.d00) then
         eosthermo_out(11) = arg2(2)*v_thermo(1,2)/eosthermo_out(11)
       else
         eosthermo_out(11) = 0.d00
       end if
       ! cs^2 and dp/dn|e
       eosthermo_out(12) = eosthermo_out(19)*(eosthermo_out(8)+1.d00)*m_n*arg2(2)
       if (eosthermo_out(12) /= 0.d00) then
         eosthermo_out(12) = 1.d00/eosthermo_out(12)
         eosthermo_out(10) = (eosthermo_out(8)+1.d00)*m_n*eosthermo_out(12)&
           -v_thermo(1,0)*eosthermo_out(11)/arg2(2)
       else
         eosthermo_out(12) = 0.d00
         eosthermo_out(10) = 0.d00
       end if
     end if

     ! additional quantities
     do iq=1,nadd_max,1
       eos_thermo_add(iq) = v_thermo(iq+dim_reg,0)
     end do

     ! compositional quantities

     if (np_max > 0) then
       idx_compo_p(1:np_max) = 0
       eos_compo_p(1:np_max) = 0.d00
     end if
     if (nq_max > 0) then
       idx_compo_q(1:nq_max) = 0
       eos_compo_q(1:nq_max,1:4) = 0.d00
     end if

! #ifdef DEBUG
!      timei = omp_get_wtime()
! #endif


     !! !$OMP END SINGLE


     if (idx_ex(2) == 1) then

       !! !$OMP SINGLE
       ! interpolation in index 1 and 2
       qx = q(ik(1))
       qy = q(ik(2))

       mat(:,-4:5,-4:5) = 0.d00
       !! !$OMP END SINGLE

       ! interpolation in index 3
       !! !$OMP DO COLLAPSE(2) &
       !! !$OMP PRIVATE(iq,igp,i1p,i2p,vec,igp3,i3p,j1,j2,j3,irp,dh,dg)
       do iq=1,np_max,1
         do igp=1,ngp,1
           i1p = idx1(igp)+m(ik(1))
           i2p = idx2(igp)+m(ik(2))

           vec(-4:5) = 0.d00
           do igp3=1,ngp3,1
             i3p = idx3(igp3)+m(ik(3))
             vec(idx3(igp3)) = 0.d00
             select case(irpl)
             case(1)
               j1 = i3p
               j2 = i1p
               j3 = i2p
             case(2)
               j1 = i1p
               j2 = i3p
               j3 = i2p
             case default
               j1 = i1p
               j2 = i2p
               j3 = i3p
             end select
             do iq2=1,np_max,1
               if (idxp_compo(j1,j2,j3,iq2) == idx_p(iq)) then
                 vec(idx3(igp3)) = tabp_compo(j1,j2,j3,iq2)
               end if
             end do
             irp = idx_arg(j1,j2,j3,ik(3))
             ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
           end do
           call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
           mat(iq,idx1(igp),idx2(igp)) = dh(0)
         end do


         call make_interp_xy(ipl,ik,qx,qy,dg,0,mat(iq,-4:5,-4:5))

         vp_compo(iq) = dg(0,0)

         !         if (vp_compo(iq) > 0.d00) then
         eos_compo_p(iq) = vp_compo(iq)
         !         end if
         idx_compo_p(iq) = idx_p(iq)

         ! quadruples
         mat(iq,-4:5,-4:5) = 0.d00
         mat3(iq,-4:5,-4:5,1:3) = 0.d00

       end do
       !! !$OMP END DO


       ! interpolation in index 3
       !! !$OMP DO COLLAPSE(2) &
       !! !$OMP PRIVATE(vec3,iq,igp,i1p,i2p,vec,iq2)
       do iq=1,nq_max,1
         do igp=1,ngp,1
           i1p = idx1(igp)+m(ik(1))
           i2p = idx2(igp)+m(ik(2))
           vec(-4:5) = 0.d00
           vec3(-4:5,1:3) = 0.d00
           do igp3=1,ngp3,1
             i3p = idx3(igp3)+m(ik(3))
             vec3(idx3(igp3),1:3) = 0.d00
             select case(irpl)
             case(1)
               j1 = i3p
               j2 = i1p
               j3 = i2p
             case(2)
               j1 = i1p
               j2 = i3p
               j3 = i2p
             case default
               j1 = i1p
               j2 = i2p
               j3 = i3p
             end select
             do iq2=1,nq_max,1
               if (idxq_compo(j1,j2,j3,iq2) == idx_q(iq)) then
                 do is=1,3,1
                   vec3(idx3(igp3),is) =&
                     tabq_compo(j1,j2,j3,3*(iq2-1)+is)
                 end do
               end if
             end do
             irp = idx_arg(j1,j2,j3,ik(3))
             ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
           end do
           do is=1,3,1
             vec(-4:5) = vec3(-4:5,is)
             call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
             mat3(iq,idx1(igp),idx2(igp),is) = dh(0)
           end do
         end do

         ! interpolation in index 1 and 2
         do is=1,3,1
           call make_interp_xy(ipl,ik,qx,qy,dg,1,mat3(iq,-4:5,-4:5,is))
           vq_compo(iq,is) = dg(0,0)
         end do

         if ((vq_compo(iq,3) > 0.d00).and.(vq_compo(iq,2) > 0.d00).and.&
           (vq_compo(iq,1) > 0.d00)) then
           eos_compo_q(iq,1) = vq_compo(iq,3)
           eos_compo_q(iq,2) = vq_compo(iq,1)
           eos_compo_q(iq,3) = vq_compo(iq,2)
           eos_compo_q(iq,4) = vq_compo(iq,1)-vq_compo(iq,2)
         else
           eos_compo_q(iq,1:4) = 0.d00
         end if
         idx_compo_q(iq) = idx_q(iq)
       end do
       !! !$OMP END DO
     end if

     !! !$OMP WORKSHARE
     ! microscopic quantities
     idx_micro(1:nm_max) = 0
     eos_micro(1:nm_max) = 0.d00
     !! !$OMP END WORKSHARE

     if (idx_ex(3) == 1) then
       !! !$OMP SINGLE
       qx = q(ik(1))
       qy = q(ik(2))

       mat(1,-4:5,-4:5) = 0.d00
       !! !$OMP END SINGLE

       ! interpolation in index 3
       !! !$OMP DO COLLAPSE(2) &
       !! !$OMP private(vec3,iq,igp,i1p,i2p,i3p,vec,iq2,j1,j2,j3)
       do iq=1,nm_max,1
         do igp=1,ngp,1
           i1p = idx1(igp)+m(ik(1))
           i2p = idx2(igp)+m(ik(2))

           vec(-4:5) = 0.d00
           do igp3=1,ngp3,1
             i3p = idx3(igp3)+m(ik(3))
             vec(idx3(igp3)) = 0.d00
             select case(irpl)
             case(1)
               j1 = i3p
               j2 = i1p
               j3 = i2p
             case(2)
               j1 = i1p
               j2 = i3p
               j3 = i2p
             case default
               j1 = i1p
               j2 = i2p
               j3 = i3p
             end select
             do iq2=1,nm_max,1
               if (idx_mic(j1,j2,j3,iq2) == idx_m(iq)) then
                 vec(idx3(igp3)) = tab_mic(j1,j2,j3,iq2)
               end if
             end do
             irp = idx_arg(j1,j2,j3,ik(3))
             ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(ik(3)))
           end do
           call get_interpol_x(m,q,ir,vec,dh,ipl,ik(3))
           mat(iq,idx1(igp),idx2(igp)) = dh(0)
         end do
         ! end do

         ! interpolation in index 1 and 2
         ! do iq=1,nm_max,1
         call make_interp_xy(ipl,ik,qx,qy,dg,0,mat(iq,-4:5,-4:5))
         v_micro(iq) = dg(0,0)
         eos_micro(iq) = v_micro(iq)
         idx_micro(iq) = idx_m(iq)
         mat(iq,-4:5,-4:5) = 0.d00
       end do
       !! !$OMP END  DO


     end if


   end if
 end if
 !! !$OMP END PARALLEL
! #ifdef DEBUG
!  timef = omp_get_wtime()
!  print '(2x,a,f12.5,a,i3,a)','TIME line 3900 : ',&
!    &    timef-timei,  '(s) with ', nthreads, ' threads'
! #endif

end SUBROUTINE eos_interpol_d3
!***********************************************************************
subroutine eos_interpol_d2(m,q,ipl,inmp,ibeta,dim_ipl,eosthermo_out)
! Stefan Typel for the CompOSE core team, version 1.18, 2017/05/23
USE compose_internal
implicit none
integer, parameter :: dim_igp=100
double precision,intent(out) :: eosthermo_out(dim_qtyt)
integer :: m(3),inmp,ibeta,i1,i2,i1p,i2p,iq,iq2,is,dim_ipl,&
     igp,ngp,igp_r,j1,j2,j3,alloc_status,ierr,&
     idx1(dim_igp),idx2(dim_igp),ipl(3),ik(3)
double precision :: q(3),qx,qy,dg(0:2,0:2),tmp
double precision, allocatable :: mat(:,:,:),mat3(:,:,:,:),&
     vp_compo(:),vq_compo(:,:),v_micro(:)

ierr = 0

allocate(mat(1:dim_ipl,-4:5,-4:5),stat=alloc_status)
if (alloc_status /= 0) then
   write(*,*) 'alloc_status = ',alloc_status
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 170
end if
mat = 0.d0

allocate(mat3(dim_ipl,-4:5,-4:5,3),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 171
end if
mat3 = 0.d0

allocate(vp_compo(np_max),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 172
end if

allocate(vq_compo(nq_max,3),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 173
end if

allocate(v_micro(nm_max),stat=alloc_status)
if (alloc_status /= 0) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 174
end if

call write_errors(ierr)

! interpolation for two dimensional table
! standard choice: two-dimensional interpolation in ik(1) and ik(2)

if (dim_idx(1) == 1) then
  ik(1:3) = [2,3,1]
else if (dim_idx(2) == 1) then
  ik(1:3) = [1,3,2]
else if (dim_idx(3) == 1) then
  ik(1:3) = [1,2,3]
else
  ierr = ierr+1
  if (ierr < dim_err) error_msg(ierr) = 163
end if

call write_errors(ierr)

! list of grid points and reference point in indices 1 and 2
igp = 0
igp_r = 0
do i1=-4,5,1 ! first index
   i1p = i1+m(ik(1))
   if ((i1p >= 1).and.(i1p <= dim_idx(ik(1)))) then
      do i2=-4,5,1 ! second index
         i2p = i2+m(ik(2))
         if ((i2p >= 1).and.(i2p <= dim_idx(ik(2)))) then
            igp = igp+1
            idx1(igp) = i1
            idx2(igp) = i2
            if (q(ik(1)) <= 0.5d00) then
               if (q(ik(2)) <= 0.5d00) then
                  if ((i1 == 0).and.(i2 == 0)) igp_r = igp
               else
                  if ((i1 == 0).and.(i2 == 1)) igp_r = igp
               end if
            else
               if (q(ik(2)) <= 0.5d00) then
                  if ((i1 == 1).and.(i2 == 0)) igp_r = igp
               else
                  if ((i1 == 1).and.(i2 == 1)) igp_r = igp
               end if
            end if
         end if
      end do
   end if
end do
ngp = igp

! for interpolation in first and second index
call get_idx_arg2(m,ik)
call get_diff_rules2(m,ipl,ik)

! thermodynamic quantities

if (ibeta == 1) then
   ! for determination of beta-equilibrium
   idx_thermo(1:nall_max) = 0
   idx_thermo(5) = 1
else
   do iq=1,nall_max,1
      idx_thermo(iq) = iq
   end do
end if

!2016/10/28
mat(1:nall_max,-4:5,-4:5) = 0.d00
do iq=1,nall_max,1
   if (idx_thermo(iq) > 0) then
      do igp=1,ngp,1
         i1p = idx1(igp)+m(ik(1))
         i2p = idx2(igp)+m(ik(2))
         if (ik(3) == 1) then
            j1 = 1
            j2 = i1p
            j3 = i2p
         else if (ik(2) == 1) then
           j1 = i1p
           j2 = 1
           j3 = i2p
         else
           j1 = i1p
           j2 = i2p
           j3 = 1
         end if
         mat(iq,idx1(igp),idx2(igp)) = tab_thermo(j1,j2,j3,iq)
      end do
   end if
end do

! interpolation in index 1 and 2
qx = q(ik(1))
qy = q(ik(2))

!2017/05/22
eos_df(1:10) = 0.d00
do iq=1,nall_max,1
   if (idx_thermo(iq) > 0) then
      call make_interp_xy(ipl,ik,qx,qy,dg,1,mat(iq,-4:5,-4:5))

      v_thermo(iq,0) = dg(0,0)
      v_thermo(iq,1) = dg(0,1)
      v_thermo(iq,2) = dg(1,0)
      v_thermo(iq,3) = dg(0,2)
      v_thermo(iq,4) = dg(2,0)
!2017/05/22
! eos_df(1): F [MeV]
! eos_df(2): dF/dT []
! eos_df(3): d^2F/(dT^2) [MeV^-1]
! eos_df(4): d^2F/(dT dn_b) [fm^3]
! eos_df(5): d^2F/(dT dY_q) []
! eos_df(6): dF/dn_b [MeV fm^3]
! eos_df(7): d^2F/(dn_b^2) [MeV fm^6]
! eos_df(8): d^2F/(dn_b dY_q) [MeV fm^3]
! eos_df(9): dF/dY_q [MeV]
! eos_df(10): d^2dF/(dY_q^2) [MeV]
      if (iq == 6) then
         eos_df(1) = (dg(0,0)+1.d00)*m_n
         if (ik(1) == 1) then
            eos_df(2) = dg(1,0)*m_n
            eos_df(3) = dg(2,0)*m_n
            eos_df(4) = dg(1,1)*m_n
            eos_df(6) = dg(0,1)*m_n
            eos_df(7) = dg(0,2)*m_n
         end if
         if (ik(2) == 1) then
            eos_df(2) = dg(1,0)*m_n
            eos_df(3) = dg(2,0)*m_n
            eos_df(5) = dg(1,1)*m_n
            eos_df(9) = dg(0,1)*m_n
            eos_df(10) = dg(0,2)*m_n
         end if
         if (ik(3) == 1) then
            eos_df(6) = dg(1,0)*m_n
            eos_df(7) = dg(2,0)*m_n
            eos_df(8) = dg(1,1)*m_n
            eos_df(9) = dg(0,1)*m_n
            eos_df(10) = dg(0,2)*m_n
         end if
      end if
   else
      v_thermo(iq,0:4) = 0.d00
   end if
end do

! mu_l lepton number chemical potential
eosthermo_out(5) = v_thermo(5,0)*m_n

if (ibeta /= 1) then

   ! p pressure
   eosthermo_out(1) = v_thermo(1,0)*arg2(2)
   ! S entropy ber baryon
   eosthermo_out(2) = v_thermo(2,0)
   ! mu_b-m_n baryon number chemical potential
   eosthermo_out(3) = v_thermo(3,0)*m_n
   ! mu_q charge chemical potential
   eosthermo_out(4) = v_thermo(4,0)*m_n

   ! F/m_n-1 free energy per baryon
   eosthermo_out(6) = v_thermo(6,0)
   ! E/m_n-1 internal energy per baryon
   eosthermo_out(7) = v_thermo(7,0)
   !eosthermo_out(7) = v_thermo(6,0)+x(1)*eosthermo_out(2)/m_n

   if (inmp == 1) then
!      eos_thermo_add(1) = v_thermo(1,1)
!      eos_thermo_add(2) = v_thermo(1,3)
     eos_thermo_add(1:2) = 0.d00
     call cpu_time(time1)
     print '("Time = ",f15.8," seconds.")',time1-time0
      return
   endif

   tmp = v_thermo(1,0)/m_n
   ! H/m_n-1  enthalpy per baryon
   eosthermo_out(8) = eosthermo_out(7)+tmp
   ! G/m_n-1 free enthalpy per baryon
   eosthermo_out(9) = eosthermo_out(6)+tmp

   ! F = f/n
   eosthermo_out(20) = (eosthermo_out(6)+1.d00)*m_n
   ! E = e/n
   eosthermo_out(21) = (eosthermo_out(7)+1.d00)*m_n
   ! H = h/n
   eosthermo_out(22) = (eosthermo_out(8)+1.d00)*m_n
   ! G = g/n
   eosthermo_out(23) = (eosthermo_out(9)+1.d00)*m_n
   ! epsilon = e
   eosthermo_out(24) = eosthermo_out(21)*arg2(2)

   ! error quantities

   ! delta (f/n) = f/n+p/n-(mu_b+y_q*(mu_q+mu_e), mu_e = mu_l-mu_q

   if (incl_l == 1) then
      ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_l)
      eos_thermo_err(1) = eosthermo_out(6)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(5))
   else
      ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_q)
      eos_thermo_err(1) = eosthermo_out(6)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(4))
   end if

   ! (delta f/n)/(f/n)
   if (eosthermo_out(6) /= -1.d00) then
      eos_thermo_err(2) = eos_thermo_err(1)/((eosthermo_out(6)+1.d00)*m_n)
   else
      eos_thermo_err(2) = 0.d00
   end if

   if (incl_l == 0) then
      ! delta (e/n) = f/n+p/n-(mu_b+y_q*mu_q)-Ts
      eos_thermo_err(3) = eosthermo_out(7)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(4))&
           -arg2(1)*eosthermo_out(2)
   else
      ! delta (e/n) = e/n+p/n-(mu_b+y_q*mu_l)-Ts
      eos_thermo_err(3) = eosthermo_out(7)*m_n+v_thermo(1,0)&
           -(eosthermo_out(3)+arg2(3)*eosthermo_out(5))&
           -arg2(1)*eosthermo_out(2)
   end if

   ! (delta e/n)/(e/n)
   if (eosthermo_out(7) /= -1.d00) then
      eos_thermo_err(4) = eos_thermo_err(3)/((eosthermo_out(7)+1.d00)*m_n)
   else
      eos_thermo_err(4) = 0.d00
   end if

   ! delta (p/n) = p/n-n*d(f/n)/dn
   eos_thermo_err(5) = v_thermo(1,0)-arg2(2)*v_thermo(6,1)*m_n
   if (v_thermo(1,0) /= 0.d00) then
      eos_thermo_err(6) = eos_thermo_err(5)/v_thermo(1,0)
   else
      eos_thermo_err(6) = 0.d00
   end if

   ! delta (s/n) = s/n+d(f/n)/dt
   eos_thermo_err(7) = eosthermo_out(2)+v_thermo(6,2)*m_n
   if (eosthermo_out(2) /= 0.d00) then
      eos_thermo_err(8) = eos_thermo_err(7)/eosthermo_out(2)
   else
      eos_thermo_err(8) = 0.d00
   end if

   ! quantities depending on derivatives

   eosthermo_out(10:19) = 0.d00
   !if (((irpl == 0).or.(irpl == 3)).and.(idx_ipl(3) /= 1)) then
   if (irpl == 0) then
      if (dim_idx(3) == 1) then
         ! beta_v
         eosthermo_out(17) = v_thermo(1,2)*arg2(2)
         ! kappa_t
         eosthermo_out(18) = arg2(2)*(arg2(2)*v_thermo(1,1)+v_thermo(1,0))
         if (eosthermo_out(18) /= 0.d00) then
            eosthermo_out(18) = 1.d00/eosthermo_out(18)
         else
            eosthermo_out(18) = 0.d00
         end if
         ! alpha_p
         eosthermo_out(16) = eosthermo_out(17)*eosthermo_out(18)
         ! c_v
         eosthermo_out(13) = v_thermo(2,2)*arg2(1)
         ! c_p
         eosthermo_out(14) = eosthermo_out(13)+eosthermo_out(16)*eosthermo_out(17)*arg2(1)/arg2(2)
         ! gamma
         if (eosthermo_out(13) /= 0.d00) then
            eosthermo_out(15) = eosthermo_out(14)/eosthermo_out(13)
         else
            eosthermo_out(15) = 1.d00
         end if
         ! kappa_s
         if (eosthermo_out(15) /= 0.d00) then
            eosthermo_out(19) = eosthermo_out(18)/eosthermo_out(15)
         else
            eosthermo_out(19) = 0.d00
         end if
         ! dp/de|n
         eosthermo_out(11) = arg2(1)*v_thermo(2,2)
         if (eosthermo_out(11) /= 0.d00) then
            eosthermo_out(11) = arg2(2)*v_thermo(1,2)/eosthermo_out(11)
         else
            eosthermo_out(11) = 0.d00
         end if
         ! cs^2 and dp/dn|e
         eosthermo_out(12) = eosthermo_out(19)*(eosthermo_out(8)+1.d00)*m_n*arg2(2)
         if (eosthermo_out(12) /= 0.d00) then
            eosthermo_out(12) = 1.d00/eosthermo_out(12)
            eosthermo_out(10) = (eosthermo_out(8)+1.d00)*m_n*eosthermo_out(12)&
            -v_thermo(1,0)*eosthermo_out(11)/arg2(2)
         else
            eosthermo_out(12) = 0.d00
            eosthermo_out(10) = 0.d00
         end if
      end if
      if (dim_idx(2) == 1) then
         ! beta_v
         eosthermo_out(17) = v_thermo(1,2)*arg2(2)
         ! c_v
         eosthermo_out(13) = v_thermo(2,2)*arg2(1)
         ! dp/de|n
         eosthermo_out(11) = arg2(1)*v_thermo(2,2)
         if (eosthermo_out(11) /= 0.d00) then
            eosthermo_out(11) = arg2(2)*v_thermo(1,2)/eosthermo_out(11)
         else
            eosthermo_out(11) = 0.d00
         end if
      end if
      if (dim_idx(1) == 1) then
         ! kappa_t
         eosthermo_out(18) = arg2(2)*(arg2(2)*v_thermo(1,2)+v_thermo(1,0))
         if (eosthermo_out(18) /= 0.d00) then
            eosthermo_out(18) = 1.d00/eosthermo_out(18)
         else
            eosthermo_out(18) = 0.d00
         end if
      end if
   end if

   ! additional quantities
   do iq=1,nadd_max,1
      eos_thermo_add(iq) = v_thermo(iq+dim_reg,0)
   end do

   ! compositional quantities

   if (np_max > 0) then
      idx_compo_p(1:np_max) = 0
      eos_compo_p(1:np_max) = 0.d00
   end if
   if (nq_max > 0) then
      idx_compo_q(1:nq_max) = 0
      eos_compo_q(1:nq_max,1:4) = 0.d00
   end if

   if (idx_ex(2) == 1) then

      ! pairs
!2016/10/28
      mat(1:np_max,-4:5,-4:5) = 0.d00
      do iq=1,np_max,1

         do igp=1,ngp,1
            i1p = idx1(igp)+m(ik(1))
            i2p = idx2(igp)+m(ik(2))
            if (ik(3) == 1) then
               j1 = 1
               j2 = i1p
               j3 = i2p
            else
               if (ik(2) == 1) then
                  j1 = i1p
                  j2 = 1
                  j3 = i2p
               else
                  j1 = i1p
                  j2 = i2p
                  j3 = 1
               end if
            end if
            do iq2=1,np_max,1
               if (idxp_compo(j1,j2,j3,iq2) == idx_p(iq)) then
                  mat(iq,idx1(igp),idx2(igp)) = tabp_compo(j1,j2,j3,iq2)
               end if
            end do
         end do
      end do


      ! interpolation in index 1 and 2
      qx = q(ik(1))
      qy = q(ik(2))
      do iq=1,np_max,1
         call make_interp_xy(ipl,ik,qx,qy,dg,0,mat(iq,-4:5,-4:5))
         vp_compo(iq) = dg(0,0)
      end do

      do iq=1,np_max,1
         !         if (vp_compo(iq) > 0.d00) then
         eos_compo_p(iq) = vp_compo(iq)
         !         end if
         idx_compo_p(iq) = idx_p(iq)
      end do


      ! quadruples
!2016/10/28
      mat(1:nq_max,-4:5,-4:5) = 0.d00
      mat3(1:nq_max,-4:5,-4:5,1:3) = 0.d00
      do iq=1,nq_max,1

         do igp=1,ngp,1
            i1p = idx1(igp)+m(ik(1))
            i2p = idx2(igp)+m(ik(2))
            if (ik(3) == 1) then
               j1 = 1
               j2 = i1p
               j3 = i2p
            else
               if (ik(2) == 1) then
                  j1 = i1p
                  j2 = 1
                  j3 = i2p
               else
                  j1 = i1p
                  j2 = i2p
                  j3 = 1
               end if
            end if

            do iq2=1,nq_max,1
               if (idxq_compo(j1,j2,j3,iq2) == idx_q(iq)) then
                  do is=1,3,1
                     mat3(iq,idx1(igp),idx2(igp),is) =&
                     tabq_compo(j1,j2,j3,3*(iq2-1)+is)
                  end do
               end if
            end do
         end do
      end do


      ! interpolation in index 1 and 2
      qx = q(ik(1))
      qy = q(ik(2))
      do iq=1,nq_max,1
         do is=1,3,1
            call make_interp_xy(ipl,ik,qx,qy,dg,0,mat3(iq,-4:5,-4:5,is))
            vq_compo(iq,is) = dg(0,0)
         end do
      end do

      do iq=1,nq_max,1
         if ((vq_compo(iq,3) > 0.d00).and.(vq_compo(iq,2) > 0.d00).and.&
              (vq_compo(iq,1) > 0.d00)) then
            eos_compo_q(iq,1) = vq_compo(iq,3)
            eos_compo_q(iq,2) = vq_compo(iq,1)
            eos_compo_q(iq,3) = vq_compo(iq,2)
            eos_compo_q(iq,4) = vq_compo(iq,1)-vq_compo(iq,2)
         else
            eos_compo_q(iq,1:4) = 0.d00
         end if
         idx_compo_q(iq) = idx_q(iq)
      end do

   end if

   ! microscopic quantities

   idx_micro(1:nm_max) = 0
   eos_micro(1:nm_max) = 0.d00

   if (idx_ex(3) == 1) then

!2016/10/28
      mat(1:nm_max,-4:5,-4:5) = 0.d00
      do iq=1,nm_max,1

         do igp=1,ngp,1
            i1p = idx1(igp)+m(ik(1))
            i2p = idx2(igp)+m(ik(2))
            if (ik(3) == 1) then
               j1 = 1
               j2 = i1p
               j3 = i2p
            else
               if (ik(2) == 1) then
                  j1 = i1p
                  j2 = 1
                  j3 = i2p
               else
                  j1 = i1p
                  j2 = i2p
                  j3 = 1
               end if
            end if

            do iq2=1,nm_max,1
               if (idx_mic(j1,j2,j3,iq2) == idx_m(iq)) then
                  mat(iq,idx1(igp),idx2(igp)) = tab_mic(j1,j2,j3,iq2)
               end if
            end do
         end do
      end do

      ! interpolation in index 1 and 2
      qx = q(ik(1))
      qy = q(ik(2))
      do iq=1,nm_max,1
         call make_interp_xy(ipl,ik,qx,qy,dg,0,mat(iq,-4:5,-4:5))
         v_micro(iq) = dg(0,0)
      end do

      do iq=1,nm_max,1
         eos_micro(iq) = v_micro(iq)
         idx_micro(iq) = idx_m(iq)
      end do
   end if

end if

if (allocated(mat)) deallocate(mat)
if (allocated(mat3)) deallocate(mat3)
if (allocated(vp_compo)) deallocate(vp_compo)
if (allocated(vq_compo)) deallocate(vq_compo)
if (allocated(v_micro)) deallocate(v_micro)

return
end SUBROUTINE eos_interpol_d2
!***********************************************************************
subroutine get_idx_arg2(m,ik)
 ! Stefan Typel for the CompOSE core team, version 1.05, 2016/10/28
 use compose_internal
 implicit none
 integer,intent(in) :: m(3),ik(3)
 integer :: i1,i2,i1p,i2p,varg(-4:5)

 ! initialisation

 idx_arg2(-4:5,-4:5,0:2) = 0

 if (eos_dim == 3) then
   if (idx_ipl(ik(3)) == 1) then
     do i1=-4,5,1 ! T
       i1p = i1+m(ik(1))
       if ((i1p >= 1).and.(i1p <= dim_idx(ik(1)))) then
         do i2=-4,5,1 ! n_b
           i2p = i2+m(ik(2))
           if ((i2p >= 1).and.(i2p <= dim_idx(ik(2)))) then
             if (irpl == 1) then
               if ((idx_arg(m(ik(3)),i1p,i2p,0) == 1).and.&
                 (idx_arg(m(ik(3))+1,i1p,i2p,0)== 1)) then
                 idx_arg2(i1,i2,0) = 1
               end if
             else
               if (irpl == 2) then
                 if ((idx_arg(i1p,m(ik(3)),i2p,0) == 1).and.&
                   (idx_arg(i1p,m(ik(3))+1,i2p,0)== 1)) then
                   idx_arg2(i1,i2,0) = 1
                 end if
               else
                 if ((idx_arg(i1p,i2p,m(ik(3)),0) == 1).and.&
                   (idx_arg(i1p,i2p,m(ik(3))+1,0)== 1)) then
                   idx_arg2(i1,i2,0) = 1
                 end if
               end if
             end if
           end if
         end do
       end if
     end do
   else
     do i1=-4,5,1 ! T
       i1p = i1+m(ik(1))
       if ((i1p >= 1).and.(i1p <= dim_idx(ik(1)))) then
         do i2=-4,5,1 ! n_b
           i2p = i2+m(ik(2))
           if ((i2p >= 1).and.(i2p <= dim_idx(ik(2)))) then
             if (irpl == 1) then
               if (idx_arg(m(ik(3)),i1p,i2p,0) == 1) then
                 idx_arg2(i1,i2,0) = 1
               end if
             else
               if (irpl == 2) then
                 if (idx_arg(i1p,m(ik(3)),i2p,0) == 1) then
                   idx_arg2(i1,i2,0) = 1
                 end if
               else
                 if (idx_arg(i1p,i2p,m(ik(3)),0) == 1) then
                   idx_arg2(i1,i2,0) = 1
                 end if
               end if
             end if
           end if
         end do
       end if
     end do
   end if
 else
   do i1=-4,5,1 ! T
     i1p = i1+m(ik(1))
     if ((i1p >= 1).and.(i1p <= dim_idx(ik(1)))) then
       do i2=-4,5,1 ! n_b
         i2p = i2+m(ik(2))
         if ((i2p >= 1).and.(i2p <= dim_idx(ik(2)))) then
           if (ik(3) == 1) then
             if (idx_arg(m(ik(3)),i1p,i2p,0) == 1) then
               idx_arg2(i1,i2,0) = 1
             end if
           else
             if (ik(3) == 2) then
               if (idx_arg(i1p,m(ik(3)),i2p,0) == 1) then
                 idx_arg2(i1,i2,0) = 1
               end if
             else
               if (idx_arg(i1p,i2p,m(ik(3)),0) == 1) then
                 idx_arg2(i1,i2,0) = 1
               end if
             end if
           end if
         end if
       end do
     end if
   end do
 end if

 ! x
 do i2=-4,5,1
   do i1=-4,4,1
     if ((idx_arg2(i1,i2,0) == 1).and.&
       (idx_arg2(i1+1,i2,0) == 1)) then
       varg(i1) = 1
     else
       varg(i1) = 0
     end if
     varg(5) = 0
   end do
   do i1=-4,4,1
     if ((varg(i1) == 1).and.(varg(i1+1) >= 1)) then
       varg(i1) = 2
     end if
   end do
   do i1=-4,3,1
     if ((varg(i1) == 2).and.(varg(i1+2) >= 1)) then
       varg(i1) = 3
     end if
   end do
   do i1=-4,2,1
     if ((varg(i1) == 3).and.(varg(i1+3) >= 1)) then
       varg(i1) = 4
     end if
   end do
   do i1=-4,5,1
     select case(varg(i1))
     case(4)
       idx_arg2(i1+4,i2,1) = 3
       idx_arg2(i1+3,i2,1) = 2
       idx_arg2(i1+2,i2,1) = 1
       if (idx_arg2(i1+1,i2,1) == 0) idx_arg2(i1+1,i2,1) = -2
       if (idx_arg2(i1  ,i2,1) == 0) idx_arg2(i1  ,i2,1) = -3
     case(3)
       if (idx_arg2(i1+3,i2,1) == 0) idx_arg2(i1+3,i2,1) =  5
       if (idx_arg2(i1+2,i2,1) == 0) idx_arg2(i1+2,i2,1) =  4
       if (idx_arg2(i1+1,i2,1) == 0) idx_arg2(i1+1,i2,1) = -4
       if (idx_arg2(i1  ,i2,1) == 0) idx_arg2(i1  ,i2,1) = -5
     case(2)
       if (idx_arg2(i1+2,i2,1) == 0) idx_arg2(i1+2,i2,1) =  7
       if (idx_arg2(i1+1,i2,1) == 0) idx_arg2(i1+1,i2,1) =  6
       if (idx_arg2(i1  ,i2,1) == 0) idx_arg2(i1  ,i2,1) = -7
     case(1)
       if (idx_arg2(i1+1,i2,1) == 0) idx_arg2(i1+1,i2,1) =  8
       if (idx_arg2(i1  ,i2,1) == 0) idx_arg2(i1  ,i2,1) = -8
     end select
   end do
 end do

 ! y
 do i1=-4,5,1
   do i2=-4,4,1
     if ((idx_arg2(i1,i2,0) == 1).and.&
       (idx_arg2(i1,i2+1,0) == 1)) then
       varg(i2) = 1
     else
       varg(i2) = 0
     end if
     varg(5) = 0
   end do
   do i2=-4,4,1
     if ((varg(i2) == 1).and.(varg(i2+1) >= 1)) then
       varg(i2) = 2
     end if
   end do
   do i2=-4,3,1
     if ((varg(i2) == 2).and.(varg(i2+2) >= 1)) then
       varg(i2) = 3
     end if
   end do
   do i2=-4,2,1
     if ((varg(i2) == 3).and.(varg(i2+3) >= 1)) then
       varg(i2) = 4
     end if
   end do
   do i2=-4,5,1
     select case(varg(i2))
     case(4)
       idx_arg2(i1,i2+4,2) = 3
       idx_arg2(i1,i2+3,2) = 2
       idx_arg2(i1,i2+2,2) = 1
       if (idx_arg2(i1,i2+1,2) == 0) idx_arg2(i1,i2+1,2) = -2
       if (idx_arg2(i1,i2  ,2) == 0) idx_arg2(i1,i2  ,2) = -3
     case(3)
       if (idx_arg2(i1,i2+3,2) == 0) idx_arg2(i1,i2+3,2) =  5
       if (idx_arg2(i1,i2+2,2) == 0) idx_arg2(i1,i2+2,2) =  4
       if (idx_arg2(i1,i2+1,2) == 0) idx_arg2(i1,i2+1,2) = -4
       if (idx_arg2(i1,i2  ,2) == 0) idx_arg2(i1,i2  ,2) = -5

     case(2)
       if (idx_arg2(i1,i2+2,2) == 0) idx_arg2(i1,i2+2,2) =  7
       if (idx_arg2(i1,i2+1,2) == 0) idx_arg2(i1,i2+1,2) =  6
       if (idx_arg2(i1,i2  ,2) == 0) idx_arg2(i1,i2  ,2) = -7
     case(1)
       if (idx_arg2(i1,i2+1,2) == 0) idx_arg2(i1,i2+1,2) =  8
       if (idx_arg2(i1,i2  ,2) == 0) idx_arg2(i1,i2+1,2) = -8
     end select
   end do
 end do

end subroutine get_idx_arg2
!***********************************************************************
SUBROUTINE get_diff_rules2(m,ipl,ik)
! Stefan Typel for the CompOSE core team, version 1.07, 2016/10/28
USE compose_internal
implicit none
integer :: m(3),i1,i2,ix,iy,ixp,iyp,i1p,i2p,iq,ir,irp,ipl(3),ik(3)


d1x(-4:5,-4:5,-4:4) = 0.d00
d2x(-4:5,-4:5,-4:4) = 0.d00
d1y(-4:5,-4:5,-4:4) = 0.d00
d2y(-4:5,-4:5,-4:4) = 0.d00


ix = ik(1)
ixp = ix+4

iy = ik(2)
iyp = iy+4

! x
do i1=-4,5,1
   i1p = i1+m(ix)
   if ((i1p >= 1).and.(i1p <= dim_idx(ix))) then
      do i2=-4,5,1
         irp = idx_arg2(i1,i2,1)
         ir = get_ipl_rule(irp,i1,ipl(ix))
         d1x(i1,i2,-4:4) = r1d(ix,i1p,-4:4,ir)
         d2x(i1,i2,-4:4) = r2d(ix,i1p,-4:4,ir)
      end do
   end if
end do
dx = tab_para(m(ix),ixp)
dx2 = dx*dx

! y
do i2=-4,5,1
   i2p = i2+m(iy)
   if ((i2p >= 1).and.(i2p <= dim_idx(iy))) then
      do i1=-4,5,1
         irp = idx_arg2(i1,i2,2)
         ir = get_ipl_rule(irp,i2,ipl(iy))
         d1y(i1,i2,-4:4) = r1d(iy,i2p,-4:4,ir)
         d2y(i1,i2,-4:4) = r2d(iy,i2p,-4:4,ir)
      end do
   end if
end do
dy = tab_para(m(iy),iyp)
dy2 = dy*dy

return
end SUBROUTINE get_diff_rules2


!***********************************************************************
subroutine get_interpol_x(m,q,ir,f,dh,ipl,idx)
 ! Stefan Typel for the CompOSE core team, version 1.04, 2016/10/28
 use compose_internal,only : r1d, r2d, tab_para
 implicit none
 double precision,intent(in) :: q(3),f(-4:5)
 double precision,intent(out) :: dh(0:2)
 integer :: m(3),i2,i2p,iq,ir(-4:5),ipl(3),idx
 double precision ::   tmp0,tmp1,tmp2,tmp3,fx(0:1),fxx(0:1),ffc(0:5)

 tmp0 = tab_para(m(idx),idx+4)
 ! derivatives
 !2016/10/28
 fx(0:1) = 0.d00
 fxx(0:1) = 0.d00
 do i2=0,1,1
   i2p = i2+m(idx)
   do iq=-4,4,1
     fx(i2)  = fx(i2) +f(i2+iq)*r1d(idx,i2p,iq,ir(i2))
     fxx(i2) = fxx(i2)+f(i2+iq)*r2d(idx,i2p,iq,ir(i2))
   end do
   ! rescaling
   fx(i2) = fx(i2)*tmp0
   fxx(i2) = fxx(i2)*tmp0*tmp0
 end do

 if (ipl(idx) == 2) then
   tmp1 = 6.d00*(f(1)-f(0))
   fxx(0) = tmp1-2.d00*fx(1)-4.d00*fx(0)
   fxx(1) = 4.d00*fx(1)+2.d00*fx(0)-tmp1
 else if (ipl(idx) == 1) then
   tmp1 = f(1)-f(0)
   fx(0) = tmp1
   fx(1) = tmp1
   fxx(0) = 0.d00
   fxx(1) = 0.d00
 end if

 ! coefficients
 ffc(0) = f(0)
 ffc(1) = fx(0)
 ffc(2) = 0.5d00*fxx(0)
 tmp1 = f(1)-f(0)-fx(0)-0.5d00*fxx(0)
 tmp2 = fx(1)-fx(0)-fxx(0)
 tmp3 = fxx(1)-fxx(0)
 ffc(3) =  10.d00*tmp1-4.d00*tmp2+0.5d00*tmp3
 ffc(4) = -15.d00*tmp1+7.d00*tmp2-tmp3
 ffc(5) =   6.d00*tmp1-3.d00*tmp2+0.5d00*tmp3

 ! interpolation
 ! functions
 dh(0) = ffc(0)+q(idx)*(ffc(1)+q(idx)*(ffc(2)&
   +q(idx)*(ffc(3)+q(idx)*(ffc(4)+q(idx)*ffc(5)))))
 ! first derivative
 dh(1) = ffc(1)+q(idx)*(2.d00*ffc(2)&
   +q(idx)*(3.d00*ffc(3)+q(idx)*(4.d00*ffc(4)+q(idx)*5.d00*ffc(5))))
 ! second derivative
 dh(2) = 2.d00*ffc(2)&
   +q(idx)*(6.d00*ffc(3)+q(idx)*(12.d00*ffc(4)+q(idx)*20.d00*ffc(5)))
 ! rescaling
 dh(1) = dh(1)/tmp0
 dh(2) = dh(2)/(tmp0*tmp0)


end SUBROUTINE get_interpol_x
!***********************************************************************


subroutine make_interp_xy(ipl,ik,qx,qy,dg,order,df00)
 use compose_internal, only : d1x, d1y , d2x, d2y,dx,dx2,dy2,dy
 double precision,intent(in) ::  df00(-4:5, -4:5), qx,qy
 double precision,intent(out) :: dg(0:2,0:2)
 integer,intent(in) :: ipl(3),ik(3),order

 integer :: i1,i2,i1p,i2p,iq,ix,iy
 double precision :: df(0:3, 0:3, -4:5, -4:5), fc(0:5, 0:5)
 double precision :: xx(-2:5),yy(-2:5)

 double precision :: tmp

 df = 0.d0
 df(0,0,-4:5,-4:5) = df00

 !2017/10/09 (sic!)
 df(1:2,0:2,-4:5,-4:5) = 0.d00
 df(0:2,1:2,-4:5,-4:5) = 0.d00

 do i2=-4,5,1
   do i1=-4,5,1
     ! first derivatives
     ! x
     do iq=-4,4,1
       i1p = i1+iq
       if ((i1p >= -4).and.(i1p <= 5)) then
          df(1,0,i1,i2) = df(1,0,i1,i2)+df(0,0,i1p,i2)*d1x(i1,i2,iq)
       end if
     end do
     ! y
     do iq=-4,4,1
       i2p = i2+iq
       if ((i2p >= -4).and.(i2p <= 5)) then
         df(0,1,i1,i2) = df(0,1,i1,i2)+df(0,0,i1,i2p)*d1y(i1,i2,iq)
       end if
     end do
     ! second derivatives
     ! xx
     do iq=-4,4,1
       i1p = i1+iq
       if ((i1p >= -4).and.(i1p <= 5)) then
         df(2,0,i1,i2) = df(2,0,i1,i2)+df(0,0,i1p,i2)*d2x(i1,i2,iq)
       end if
     end do
     ! yy
     do iq=-4,4,1
       i2p = i2+iq
       if ((i2p >= -4).and.(i2p <= 5)) then
         df(0,2,i1,i2) = df(0,2,i1,i2)+df(0,0,i1,i2p)*d2y(i1,i2,iq)
       end if
     end do
   end do
 end do

 ! mixed derivatives
 do i2=0,1,1
   do i1=0,1,1
     ! second derivatives
     ! xy
     do iq=-4,4,1
       i1p = i1+iq
       if ((i1p >= -4).and.(i1p <= 5)) then
         df(1,1,i1,i2) = df(1,1,i1,i2)+df(0,1,i1p,i2)*d1x(i1,i2,iq)
       end if
     end do
     do iq=-4,4,1
       i2p = i2+iq
       if ((i2p >= -4).and.(i2p <= 5)) then
         df(1,1,i1,i2) = df(1,1,i1,i2)+df(1,0,i1,i2p)*d1y(i1,i2,iq)
       end if
     end do
     df(1,1,i1,i2) = 0.5d00*df(1,1,i1,i2)
     ! third derivatives
     ! xxy
     do iq=-4,4,1
       i2p = i2+iq
       if ((i2p >= -4).and.(i2p <= 5)) then
         df(2,1,i1,i2) = df(2,1,i1,i2)+df(2,0,i1,i2p)*d1y(i1,i2,iq)
       end if
     end do
     ! xyy
     do iq=-4,4,1
       i1p = i1+iq
       if ((i1p >= -4).and.(i1p <= 5)) then
         df(1,2,i1,i2) = df(1,2,i1,i2)+df(0,2,i1p,i2)*d1x(i1,i2,iq)
       end if
     end do
     ! fourth derivative
     ! xxyy
     do iq=-4,4,1
       i1p = i1+iq
       if ((i1p >= -4).and.(i1p <= 5)) then
         df(2,2,i1,i2) = df(2,2,i1,i2)+df(0,2,i1p,i2)*d2x(i1,i2,iq)
       end if
     end do
     do iq=-4,4,1
       i2p = i2+iq
       if ((i2p >= -4).and.(i2p <= 5)) then
         df(2,2,i1,i2) = df(2,2,i1,i2)+df(2,0,i1,i2p)*d2y(i1,i2,iq)
       end if
     end do
     df(2,2,i1,i2) = 0.5d00*df(2,2,i1,i2)
   end do
 end do


 ! rescaling
 if (1 == 1) then
   do i2=0,1,1
     do i1=0,1,1
       df(1,0,i1,i2) = df(1,0,i1,i2)*dx
       df(0,1,i1,i2) = df(0,1,i1,i2)*dy
       df(2,0,i1,i2) = df(2,0,i1,i2)*dx2
       df(0,2,i1,i2) = df(0,2,i1,i2)*dy2
       df(1,1,i1,i2) = df(1,1,i1,i2)*dx*dy
       df(2,1,i1,i2) = df(2,1,i1,i2)*dx2*dy
       df(1,2,i1,i2) = df(1,2,i1,i2)*dx*dy2
       df(2,2,i1,i2) = df(2,2,i1,i2)*dx2*dy2
     end do
   end do
 end if

 ! print *,dx,dy,dx2,dy2,dx*dy,dx2,dy2

 !        print *,df(1,0,:,:)
 !        print *,df(0,1,:,:)
 !        print *,df(2,0,:,:)
 !        print *,df(0,2,:,:)
 !        print *,df(1,1,:,:)
 !        print *,df(2,1,:,:)
 !        print *,df(1,2,:,:)
 !        print *,df(2,2,:,:)

 !        stop



 !???
 if (0 == 1) then
   if (ipl(ik(2)) == 2) then
     do i2=0,2,1
       do i1=0,1,1
         tmp = 6.d00*(df(i2,0,i1,1)-df(i2,0,i1,0))
         df(i2,2,i1,0) = tmp-2.d00*df(i2,1,i1,1)-4.d00*df(i2,1,i1,0)
         df(i2,2,i1,1) = 4.d00*df(i2,1,i1,1)+2.d00*df(i2,1,i1,0)-tmp
       end do
     end do
   end if

   if (ipl(ik(1)) == 2) then
     do i1=0,1,1
       do i2=0,2,1
         tmp = 6.d00*(df(0,i2,1,i1)-df(0,i2,0,i1))
         df(2,i2,0,i1) = tmp-2.d00*df(1,i2,1,i1)-4.d00*df(1,i2,0,i1)
         df(2,i2,1,i1) = 4.d00*df(1,i2,1,i1)+2.d00*df(1,i2,0,i1)-tmp
       end do
     end do
   end if

   if (ipl(ik(2)) == 1) then
     !2016/10/28
     df(0:2,2,0:1,0) = 0.d00
     df(0:2,2,0:1,1) = 0.d00
     do i1=0,1,1
       do i2=0,2,1
         tmp = df(i2,0,i1,1)-df(i2,0,i1,0)
         df(i2,1,i1,0) = tmp
         df(i2,1,i1,1) = tmp
         !         df(i2,2,i1,0) = 0.d00
         !         df(i2,2,i1,1) = 0.d00
       end do
     end do
   end if

   if (ipl(ik(1)) == 1) then
     !2016/10/28
     df(2,0:2,0,0:1) = 0.d00
     df(2,0:2,1,0:1) = 0.d00
     do i1=0,1,1
       do i2=0,2,1
         tmp = df(0,i2,1,i1)-df(0,i2,0,i1)
         df(1,i2,0,i1) = tmp
         df(1,i2,1,i1) = tmp
       end do
     end do
   end if
 end if


 fc(0,0) = df(0,0,0,0)
 fc(0,1) = df(0,1,0,0)
 fc(0,2) = 0.5d00*df(0,2,0,0)
 fc(0,3) = -10.d00*(df(0,0,0,0)-df(0,0,0,1))&
   -6.d00*df(0,1,0,0)-4.d00*df(0,1,0,1)&
   -1.5d00*df(0,2,0,0)+0.5d00*df(0,2,0,1)
 fc(0,4) =  15.d00*(df(0,0,0,0)-df(0,0,0,1))&
   +8.d00*df(0,1,0,0)+7.d00*df(0,1,0,1)&
   +1.5d00*df(0,2,0,0)-1.d00*df(0,2,0,1)
 fc(0,5) =  -6.d00*(df(0,0,0,0)-df(0,0,0,1))&
   -3.d00*(df(0,1,0,0)+df(0,1,0,1))&
   -0.5d00*(df(0,2,0,0)-df(0,2,0,1))
 fc(1,0) = df(1,0,0,0)
 fc(1,1) = df(1,1,0,0)
 fc(1,2) = 0.5d00*df(1,2,0,0)
 fc(1,3) = -10.d00*(df(1,0,0,0)-df(1,0,0,1))&
   -6.d00*df(1,1,0,0)-4.d00*df(1,1,0,1)&
   -1.5d00*df(1,2,0,0)+0.5d00*df(1,2,0,1)
 fc(1,4) =  15.d00*(df(1,0,0,0)-df(1,0,0,1))&
   +8.d00*df(1,1,0,0)+7.d00*df(1,1,0,1)&
   +1.5d00*df(1,2,0,0)-1.d00*df(1,2,0,1)
 fc(1,5) =  -6.d00*(df(1,0,0,0)-df(1,0,0,1))&
   -3.d00*(df(1,1,0,0)+df(1,1,0,1))&
   -0.5d00*(df(1,2,0,0)-df(1,2,0,1))
 fc(2,0) = 0.5d00*df(2,0,0,0)
 fc(2,1) = 0.5d00*df(2,1,0,0)
 fc(2,2) = 0.25d00*df(2,2,0,0)
 fc(2,3) = -5.d00*(df(2,0,0,0)-df(2,0,0,1))&
   -3.d00*df(2,1,0,0)-2.d00*df(2,1,0,1)&
   -0.75d00*df(2,2,0,0)+0.25d00*df(2,2,0,1)
 fc(2,4) =  7.5d00*(df(2,0,0,0)-df(2,0,0,1))&
   +4.d00*df(2,1,0,0)+3.5d00*df(2,1,0,1)&
   +0.75d00*df(2,2,0,0)-0.5d00*df(2,2,0,1)
 fc(2,5) = -3.d00*(df(2,0,0,0)-df(2,0,0,1))&
   -1.5d00*(df(2,1,0,0)+df(2,1,0,1))&
   -0.25d00*(df(2,2,0,0)-df(2,2,0,1))
 fc(3,0) =  -10.d00*(df(0,0,0,0)-df(0,0,1,0))&
   -6.d00*df(1,0,0,0)-4.d00*df(1,0,1,0)&
   -1.5d00*df(2,0,0,0)+0.5d00*df(2,0,1,0)
 fc(3,1) = -10.d00*(df(0,1,0,0)-df(0,1,1,0))&
   -6.d00*df(1,1,0,0)-4.d00*df(1,1,1,0)&
   -1.5d00*df(2,1,0,0)+0.5d00*df(2,1,1,0)
 fc(3,2) = -5.d00*(df(0,2,0,0)-df(0,2,1,0))&
   -3.d00*df(1,2,0,0)-2.d00*df(1,2,1,0)&
   -0.75d00*df(2,2,0,0)+0.25d00*df(2,2,1,0)
 fc(3,3) =  100.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   +60.d00*(df(1,0,0,0)-df(1,0,0,1))+40.d00*(df(1,0,1,0)-df(1,0,1,1))&
   +60.d00*(df(0,1,0,0)-df(0,1,1,0))+40.d00*(df(0,1,0,1)-df(0,1,1,1))&
   +15.d00*(df(2,0,0,0)-df(2,0,0,1))-5.d00*(df(2,0,1,0)-df(2,0,1,1))&
   +36.d00*df(1,1,0,0)+24.d00*(df(1,1,1,0)+df(1,1,0,1))+16.d00*df(1,1,1,1)&
   +15.d00*(df(0,2,0,0)-df(0,2,1,0))-5.d00*(df(0,2,0,1)-df(0,2,1,1))&
   +9.d00*(df(2,1,0,0)+df(1,2,0,0))-3.d00*(df(2,1,1,0)+df(1,2,0,1))&
   +6.d00*(df(2,1,0,1)+df(1,2,1,0))-2.d00*(df(2,1,1,1)+df(1,2,1,1))&
   +2.25d00*df(2,2,0,0)-0.75d00*(df(2,2,1,0)+df(2,2,0,1))&
   +0.25d00*df(2,2,1,1)
 fc(3,4) = -150.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   -90.d00*(df(1,0,0,0)-df(1,0,0,1))-60.d00*(df(1,0,1,0)-df(1,0,1,1))&
   -80.d00*(df(0,1,0,0)-df(0,1,1,0))-70.d00*(df(0,1,0,1)-df(0,1,1,1))&
   -22.5d00*(df(2,0,0,0)-df(2,0,0,1))+7.5d00*(df(2,0,1,0)-df(2,0,1,1))&
   -48.d00*df(1,1,0,0)-32.d00*df(1,1,1,0)&
   -42.d00*df(1,1,0,1)-28.d00*df(1,1,1,1)&
   -15.d00*(df(0,2,0,0)-df(0,2,1,0))+10.d00*(df(0,2,0,1)-df(0,2,1,1))&
   -12.d00*df(2,1,0,0)+4.d00*df(2,1,1,0)&
   -10.5d00*df(2,1,0,1)+3.5d00*df(2,1,1,1)&
   -9.d00*df(1,2,0,0)-6.d00*(df(1,2,1,0)-df(1,2,0,1))+4.d00*df(1,2,1,1)&
   -2.25d00*df(2,2,0,0)+0.75d00*df(2,2,1,0)&
   +1.5d00*df(2,2,0,1)-0.5d00*df(2,2,1,1)
 fc(3,5) =   60.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   +36.d00*(df(1,0,0,0)-df(1,0,0,1))+24.d00*(df(1,0,1,0)-df(1,0,1,1))&
   +30.d00*(df(0,1,0,0)-df(0,1,1,0)+df(0,1,0,1)-df(0,1,1,1))&
   +9.d00*(df(2,0,0,0)-df(2,0,0,1))-3.d00*(df(2,0,1,0)-df(2,0,1,1))&
   +18.d00*(df(1,1,0,0)+df(1,1,0,1))+12.d00*(df(1,1,1,0)+df(1,1,1,1))&
   +5.d00*(df(0,2,0,0)-df(0,2,1,0)-df(0,2,0,1)+df(0,2,1,1))&
   +4.5d00*(df(2,1,0,0)+df(2,1,0,1))-1.5d00*(df(2,1,1,0)+df(2,1,1,1))&
   +3.d00*(df(1,2,0,0)-df(1,2,0,1))+2.d00*(df(1,2,1,0)-df(1,2,1,1))&
   +0.75d00*(df(2,2,0,0)-df(2,2,0,1))&
   -0.25d00*(df(2,2,1,0)-df(2,2,1,1))
 fc(4,0) =   15.d00*(df(0,0,0,0)-df(0,0,1,0))&
   +8.d00*df(1,0,0,0)+7.d00*df(1,0,1,0)&
   +1.5d00*df(2,0,0,0)-1.d00*df(2,0,1,0)
 fc(4,1) = 15.d00*(df(0,1,0,0)-df(0,1,1,0))&
   +8.d00*df(1,1,0,0)+7.d00*df(1,1,1,0)&
   +1.5d00*df(2,1,0,0)-1.d00*df(2,1,1,0)
 fc(4,2) = 7.5d00*(df(0,2,0,0)-df(0,2,1,0))&
   +4.d00*df(1,2,0,0)+3.5d00*df(1,2,1,0)&
   +0.75d00*df(2,2,0,0)-0.5d00*df(2,2,1,0)
 fc(4,3) = -150.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   -80.d00*(df(1,0,0,0)-df(1,0,0,1))-70.d00*(df(1,0,1,0)-df(1,0,1,1))&
   -90.d00*(df(0,1,0,0)-df(0,1,1,0))-60.d00*(df(0,1,0,1)-df(0,1,1,1))&
   -15.d00*(df(2,0,0,0)-df(2,0,0,1))+10.d00*(df(2,0,1,0)-df(2,0,1,1))&
   -48.d00*df(1,1,0,0)-42.d00*df(1,1,1,0)&
   -32.d00*df(1,1,0,1)-28.d00*df(1,1,1,1)&
   -22.5d00*(df(0,2,0,0)-df(0,2,1,0))+7.5d00*(df(0,2,0,1)-df(0,2,1,1))&
   -9.d00*df(2,1,0,0)+6.d00*(df(2,1,1,0)-df(2,1,0,1))+4.d00*df(2,1,1,1)&
   -12.d00*df(1,2,0,0)-10.5d00*df(1,2,1,0)&
   +4.d00*df(1,2,0,1)+3.5d00*df(1,2,1,1)&
   -2.25d00*df(2,2,0,0)+1.5d00*df(2,2,1,0)&
   +0.75d00*df(2,2,0,1)-0.5d00*df(2,2,1,1)
 fc(4,4) =  225.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   +120.d00*(df(1,0,0,0)-df(1,0,0,1))+105.d00*(df(1,0,1,0)-df(1,0,1,1))&
   +120.d00*(df(0,1,0,0)-df(0,1,1,0))+105.d00*(df(0,1,0,1)-df(0,1,1,1))&
   +22.5d00*(df(2,0,0,0)-df(2,0,0,1))-15.d00*(df(2,0,1,0)-df(2,0,1,1))&
   +64.d00*df(1,1,0,0)+56.d00*(df(1,1,1,0)+df(1,1,0,1))+49.d00*df(1,1,1,1)&
   +22.5d00*(df(0,2,0,0)-df(0,2,1,0))-15.d00*(df(0,2,0,1)-df(0,2,1,1))&
   +12.d00*(df(2,1,0,0)+df(1,2,0,0))-8.d00*(df(2,1,1,0)+df(1,2,0,1))&
   +10.5d00*(df(2,1,0,1)+df(1,2,1,0))-7.d00*(df(2,1,1,1)+df(1,2,1,1))&
   +2.25d00*df(2,2,0,0)-1.5d00*(df(2,2,1,0)+df(2,2,0,1))&
   +1.d00*df(2,2,1,1)
 fc(4,5) =  -90.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   -48.d00*(df(1,0,0,0)-df(1,0,0,1))-42.d00*(df(1,0,1,0)-df(1,0,1,1))&
   -45.d00*(df(0,1,0,0)-df(0,1,1,0)+df(0,1,0,1)-df(0,1,1,1))&
   -9.d00*(df(2,0,0,0)-df(2,0,0,1))+6.d00*(df(2,0,1,0)-df(2,0,1,1))&
   -24.d00*(df(1,1,0,0)+df(1,1,0,1))-21.d00*(df(1,1,1,0)+df(1,1,1,1))&
   -7.5d00*(df(0,2,0,0)-df(0,2,1,0)-df(0,2,0,1)+df(0,2,1,1))&
   -4.5d00*(df(2,1,0,0)+df(2,1,0,1))+3.d00*(df(2,1,1,0)+df(2,1,1,1))&
   -4.d00*df(1,2,0,0)-3.5d00*(df(1,2,1,0)-df(1,2,1,1))+4.d00*df(1,2,0,1)&
   -0.75d00*(df(2,2,0,0)-df(2,2,0,1))&
   +0.5d00*(df(2,2,1,0)-df(2,2,1,1))
 fc(5,0) =   -6.d00*(df(0,0,0,0)-df(0,0,1,0))&
   -3.d00*(df(1,0,0,0)+df(1,0,1,0))&
   -0.5d00*(df(2,0,0,0)-df(2,0,1,0))
 fc(5,1) = -6.d00*(df(0,1,0,0)-df(0,1,1,0))&
   -3.d00*(df(1,1,0,0)+df(1,1,1,0))&
   -0.5d00*(df(2,1,0,0)-df(2,1,1,0))
 fc(5,2) = -3.d00*df(0,2,0,0)+3.d00*df(0,2,1,0)&
   -1.5d00*df(1,2,0,0)-1.5d00*df(1,2,1,0)&
   -0.25d00*df(2,2,0,0)+0.25d00*df(2,2,1,0)
 fc(5,3) =   60.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   +30.d00*(df(1,0,0,0)+df(1,0,1,0)-df(1,0,0,1)-df(1,0,1,1))&
   +36.d00*(df(0,1,0,0)-df(0,1,1,0))+24.d00*(df(0,1,0,1)-df(0,1,1,1))&
   +5.d00*(df(2,0,0,0)-df(2,0,1,0)-df(2,0,0,1)+df(2,0,1,1))&
   +18.d00*(df(1,1,0,0)+df(1,1,1,0))+12.d00*(df(1,1,0,1)+df(1,1,1,1))&
   +9.d00*(df(0,2,0,0)-df(0,2,1,0))-3.d00*(df(0,2,0,1)-df(0,2,1,1))&
   +3.d00*(df(2,1,0,0)-df(2,1,1,0))+2.d00*(df(2,1,0,1)-df(2,1,1,1))&
   +4.5d00*df(1,2,0,0)+4.5d00*df(1,2,1,0)&
   -1.5d00*(df(1,2,0,1)+df(1,2,1,1))&
   +0.75d00*(df(2,2,0,0)-df(2,2,1,0))&
   -0.25d00*(df(2,2,0,1)-df(2,2,1,1))
 fc(5,4) =  -90.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   -45.d00*(df(1,0,0,0)+df(1,0,1,0)-df(1,0,0,1)-df(1,0,1,1))&
   -48.d00*(df(0,1,0,0)-df(0,1,1,0))-42.d00*(df(0,1,0,1)-df(0,1,1,1))&
   -7.5d00*(df(2,0,0,0)-df(2,0,1,0)-df(2,0,0,1)+df(2,0,1,1))&
   -24.d00*(df(1,1,0,0)+df(1,1,1,0))-21.d00*(df(1,1,0,1)+df(1,1,1,1))&
   -9.d00*(df(0,2,0,0)-df(0,2,1,0))+6.d00*(df(0,2,0,1)-df(0,2,1,1))&
   -4.d00*(df(2,1,0,0)-df(2,1,1,0))-3.5d00*(df(2,1,0,1)-df(2,1,1,1))&
   -4.5d00*(df(1,2,0,0)+df(1,2,1,0))+3.d00*(df(1,2,0,1)+df(1,2,1,1))&
   -0.75d00*(df(2,2,0,0)-df(2,2,1,0))&
   +0.5d00*(df(2,2,0,1)-df(2,2,1,1))
 fc(5,5) =   36.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
   +18.d00*(df(1,0,0,0)+df(1,0,1,0)-df(1,0,0,1)-df(1,0,1,1))&
   +18.d00*(df(0,1,0,0)-df(0,1,1,0)+df(0,1,0,1)-df(0,1,1,1))&
   +3.d00*(df(2,0,0,0)-df(2,0,1,0)-df(2,0,0,1)+df(2,0,1,1))&
   +9.d00*(df(1,1,0,0)+df(1,1,1,0)+df(1,1,0,1)+df(1,1,1,1))&
   +3.d00*(df(0,2,0,0)-df(0,2,1,0)-df(0,2,0,1)+df(0,2,1,1))&
   +1.5d00*(df(2,1,0,0)-df(2,1,1,0)+df(2,1,0,1)-df(2,1,1,1))&
   +1.5d00*(df(1,2,0,0)+df(1,2,1,0)-df(1,2,0,1)-df(1,2,1,1))&
   +0.25d00*(df(2,2,0,0)-df(2,2,1,0)-df(2,2,0,1)+df(2,2,1,1))




 xx(-2:-1) = 0.d00
 yy(-2:-1) = 0.d00
 xx(0) = 1.d00
 yy(0) = 1.d00
 do ix=1,5,1
   xx(ix) = xx(ix-1)*qx
   yy(ix) = yy(ix-1)*qy
 end do

 dg(0:2,0:2) = 0.d00

 if (order == 0) then
   ! function
   do ix=0,5,1
     do iy=0,5,1
       dg(0,0) = dg(0,0)+fc(ix,iy)*xx(ix)*yy(iy)
     end do
   end do
 else
   ! function and derivatives
   do ix=0,5,1
     do iy=0,5,1
       dg(0,0) = dg(0,0)+fc(ix,iy)*xx(ix)*yy(iy)
       dg(1,0) = dg(1,0)+fc(ix,iy)*xx(ix-1)*yy(iy)*dble(ix)
       dg(0,1) = dg(0,1)+fc(ix,iy)*xx(ix)*yy(iy-1)*dble(iy)
       dg(2,0) = dg(2,0)+fc(ix,iy)*xx(ix-2)*yy(iy)*dble(ix*(ix-1))
       dg(1,1) = dg(1,1)+fc(ix,iy)*xx(ix-1)*yy(iy-1)*dble(iy*ix)
       dg(0,2) = dg(0,2)+fc(ix,iy)*xx(ix)*yy(iy-2)*dble(iy*(iy-1))
       dg(2,1) = dg(2,1)+fc(ix,iy)*xx(ix-2)*yy(iy-1)*dble(ix*(ix-1)*iy)
       dg(1,2) = dg(1,2)+fc(ix,iy)*xx(ix-1)*yy(iy-2)*dble(ix*iy*(iy-1))
       dg(2,2) = dg(2,2)+fc(ix,iy)*xx(ix-2)*yy(iy-2)*&
         dble(ix*(ix-1)*iy*(iy-1))
     end do
   end do

   ! rescaling
   if (dx > 0.d00) then
     dg(1,0) = dg(1,0)/dx
     dg(2,0) = dg(2,0)/dx**2
   end if
   if (dy > 0.d00) then
     dg(0,1) = dg(0,1)/dy
     dg(0,2) = dg(0,2)/dy**2
   end if
   if ((dx > 0.d00).and.(dy > 0.)) then
     dg(1,1) = dg(1,1)/(dx*dy)
     dg(2,1) = dg(2,1)/((dx**2)*dy)
     dg(1,2) = dg(1,2)/(dx*(dy**2))
     dg(2,2) = dg(2,2)/(dx*dy)**2
   end if
 end if

end subroutine make_interp_xy

!***********************************************************************
SUBROUTINE run_terminal(iwr)
! Stefan Typel for the CompOSE core team, version 1.12, 2017/12/13
 use compose_internal
 use m_get_tables, only : init_eos_table_term


implicit none
integer :: iwr,iterm,iinit

write(*,*)
write(*,*) ' *************************************************'
write(*,*) ' *              Welcome to CompOSE               *'
write(*,*) ' * CompStar Online Supernovae Equations of State *'
write(*,*) ' *                 Version 2.17                  *'
write(*,*) ' *                  2018/09/07                   *'
write(*,*) ' *************************************************'
write(*,*)
write(*,*) ' This program helps to generate user-specified EoS tables'
write(*,*) ' from the EoS tables provided by the CompOSE database at'
write(*,*) ' compose.obspm.fr.'
write(*,*)
write(*,*) ' Please select the task number from the following list:'
iterm = 0
do while (((iterm < 1).or.(iterm > 3)).and.(iterm /= 999))
   write(*,*)
   write(*,*) ' Task 1: Selection of Output Quantities'
   write(*,*) '         (Creates files eos.quantities and eos.init, if not existing)'
   write(*,*) ' Task 2: Definition of Tabulation Scheme and Parameter Values'
   write(*,*) '         (Creates files eos.parameters and eos.init, if not existing)'
   write(*,*) ' Task 3: Generation of EoS Table'
   write(*,*) '         (Creates files eos.table, eos.report,'
   write(*,*) '         eos.beta, if possible, and eos.init, if not existing)'
   write(*,*)
   read(5,*) iterm
end do

call init_eos_table_term(iwr,iterm,iinit)


select case (iterm)
case (1)
   call init_quant()
case (2)
   call init_para()
case (3)
   call get_eos_table_term(iwr,iinit)
case (999)
   write(*,*) 'New File eos.init created'
case default
   write(*,*)
   write(*,*) 'Execution of Programm Cancelled'
   write(*,*)
end select

end SUBROUTINE run_terminal
!***********************************************************************

SUBROUTINE init_quant()
! Stefan Typel for the CompOSE core team, version 1.08, 2017/11/16
USE compose_internal
implicit none
integer :: iunit2,iunit3,ierror,inew,iv,ip,iq,im,ierr,ivar(4),ibeta2,&
     idxp(dim_ip),idxq(dim_iq),idxm(dim_im),&
     idxpp(dim_ip),idxqq(dim_iq),idxmm(dim_im),idxerr(nerr_max)
character(11) :: part(0:999)
character(33) :: micr(0:99)
character(27) :: qerr(nerr_max)

do iv=0,999,1
   part(iv) = 'not defined'
end do

part(0)   = 'electron   '
part(1)   = 'muon       '
part(10)  = 'neutron    '
part(11)  = 'proton     '
part(20)  = 'Delta^-    '
part(21)  = 'Delta^0    '
part(22)  = 'Delta^+    '
part(23)  = 'Delta^++   '
part(100) = 'Lambda     '
part(110) = 'Sigma^-    '
part(111) = 'Sigma^0    '
part(112) = 'Sigma^+    '
part(120) = 'Xi^-       '
part(121) = 'Xi^0       '
part(200) = 'omega      '
part(210) = 'sigma      '
part(220) = 'eta        '
part(230) = 'eta^prime  '
part(300) = 'rho^-      '
part(301) = 'rho^0      '
part(302) = 'rho^+      '
part(310) = 'delta^-    '
part(311) = 'delta^0    '
part(312) = 'delta^+    '
part(320) = 'pi^-       '
part(321) = 'pi^0       '
part(322) = 'pi^+       '
part(400) = 'phi        '
part(410) = 'sigma_s    '
part(420) = 'K^-        '
part(421) = 'K^0        '
part(422) = 'Kbar^0     '
part(423) = 'Kbar^+     '
part(424) = 'K^- (thermal)'
part(425) = 'K^- (condensate)'
part(500) = 'u quark    '
part(501) = 'd quark    '
part(502) = 's quark    '
part(600) = 'photon     '
part(700) = 'nn (1S0)   '
part(701) = 'np (1S0)   '
part(702) = 'pp (1S0)   '
part(703) = 'np (3S1)   '

do iv=0,99,1
   micr(iv) = 'quantity defined in data sheet   '
end do

micr(40) = 'scaled Landau mass m^L/m    []   '
micr(41) = 'scaled Dirac mass  m^D/m    []   '
micr(50) = 'single-particle potential U [MeV]'
micr(51) = 'vector self-energy V        [MeV]'
micr(52) = 'scalar self-energy S        [MeV]'
micr(60) = 'gap Delta                   [MeV]'

qerr(1) = 'absolute error in f/n [MeV]'
qerr(2) = 'relative error in f/n []   '
qerr(3) = 'absolute error in e/n [MeV]'
qerr(4) = 'relative error in e/n []   '
qerr(5) = 'absolute error in p/n [MeV]'
qerr(6) = 'relative error in p/n []   '
qerr(7) = 'absolute error in s/n []   '
qerr(8) = 'relative error in s/n []   '

iunit2 = 31

open(unit=iunit2,file='eos.init',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   read(iunit2,*) (ivar(iv),iv=1,4,1),irpl
   do iv=1,4,1
      read(iunit2,*) para_min(iv),para_max(iv)
   end do
   read(iunit2,*) incl_l,nadd_max,nall_max,ibeta2
   read(iunit2,*) np_max,nq_max
   read(iunit2,*) (idxp(iv),iv=1,np_max,1)
   read(iunit2,*) (idxq(iv),iv=1,nq_max,1)
   read(iunit2,*) nm_max
   read(iunit2,*) (idxm(iv),iv=1,nm_max,1)
else
   write(*,*)
   write(*,*) ' There is no file eos.init.'
   write(*,*) ' Please restart the program compose with task 1'
   write(*,*) ' and generate a new file eos.init.'
   write(*,*)
   stop
end if
close(unit=iunit2)


iunit3 = 32

open(unit=iunit3,file='eos.quantities',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   write(*,*)
   write(*,*) ' The file eos.quantities exists already.'
   write(*,*) ' Do you want to generate a new file eos.quantities?'
   write(*,*) ' Please select:'
   write(*,*) ' 1: Yes'
   write(*,*) ' else: No'
   read(5,*) inew
else
   inew = 1
end if
close(unit=iunit3)

if (inew == 1) then
   open(unit=iunit3,file='eos.quantities',&
        status='unknown',action='write',iostat=ierror)

   write(*,*)
   write(*,*) ' Data are available from this(these) group(s):'
   if (nadd_max > 0) then
      write(*,*) ' - regular and additional thermodynamic quantities'
   else
      write(*,*) ' - regular thermodynamic quantities'
   end if
!2017/05/22
   if (irpl == 0) then
      write(*,*) ' - function values and derivatives of the free energy per baryon'
   end if
   if ((np_max > 0).or.(nq_max > 0)) then
      write(*,*) ' - data on the chemical composition'
   end if
   if (nm_max > 0) then
      write(*,*) ' - microscopic quantities'
   end if
   write(*,*)
   write(*,*) ' The following regular thermodynamic quantities are available:'
   write(*,*) ' (See also table 7.1 of the manual.)'
   write(*,*)
   write(*,*) ' quantity:                                          unit:       index:'
   write(*,*) ' pressure p                                         [MeV fm^-3]    1'
   write(*,*) ' entropy per baryon S                               []             2'
   write(*,*) ' shifted baryon chemical potential mu_b-m_n         [MeV]          3'
   write(*,*) ' charge chemical potential mu_q                     [MeV]          4'
   write(*,*) ' lepton chemical potential mu_l                     [MeV]          5'
   write(*,*) ' scaled free energy per baryon F/m_n-1              []             6'
   write(*,*) ' scaled internal energy per baryon E/m_n-1          []             7'
   write(*,*) ' scaled enthalpy energy per baryon H/m_n-1          []             8'
   write(*,*) ' scaled free enthalpy per baryon G/m_n-1            []             9'
   if ((ivar(1) == 1).and.(ivar(2) == 1)) then
      write(*,*) ' derivative dp/dn_b|E                               [MeV]         10'
      write(*,*) ' derivative p/dE|n_b                                [fm^-3]       11'
      write(*,*) ' square of speed of sound (c_s)^2                   []            12'
   end if
   if (ivar(1) == 1)&
        write(*,*) ' specific heat capacity at constant volume c_V      []            13'
   if ((ivar(1) == 1).and.(ivar(2) == 1)) then
      write(*,*) ' specific heat capacity at constant pressure c_p    []            14'
      write(*,*) ' adiabatic index Gamma                              []            15'
      write(*,*) ' expansion coefficient at constant pressure alpha_p [MeV^-1]      16'
      write(*,*) ' tension coefficient at constant volume beta_V      [fm^-3]       17'
   end if
   if (ivar(2) == 1)&
        write(*,*) ' isothermal compressibility kappa_T                 [MeV^-1 fm^3] 18'
   if ((ivar(1) == 1).and.(ivar(2) == 1))&
        write(*,*) ' isentropic compressibility kappa_S                 [MeV^-1 fm^3] 19'
   write(*,*) ' free energy per baryon F                           [MeV]         20'
   write(*,*) ' internal energy per baryon E                       [MeV]         21'
   write(*,*) ' enthalpy per baryon H                              [MeV]         22'
   write(*,*) ' free enthalpy per baryon G                         [MeV]         23'
   write(*,*) ' energy density                       [MeV/fm^3]         24'
   write(*,*)
   n_qty = -1
   do while ((n_qty < 0).or.(n_qty > dim_qtyt))
      write(*,*) ' How many regular thermodynamic quantities do you want to select for the file eos.table?'
      read(*,*) n_qty
      if ((n_qty < 0).or.(n_qty > dim_qtyt)) then
         write(*,*) ' Number out of range!'
      end if
   end do
!2017/05/22
   if (n_qty > 0) then
      write(*,*) ' Please select the indices of the thermodynamic quantities'
      write(*,*) ' in the order of your choice.'
      do iv=1,n_qty,1
         idx_qty(iv) = 0
         do while ((idx_qty(iv) < 1).or.(idx_qty(iv) > dim_qtyt))
            write(*,*) 'Index #',iv,'?'
            read(5,*) idx_qty(iv)
            if ((idx_qty(iv) < 1).or.(idx_qty(iv) > dim_qtyt)) then
               write(*,*) ' Number out of range!'
            end if
         end do
      end do
   end if

   if (nadd_max > 0) then
      write(*,*)
      if (nadd_max == 1) then
         write(*,*) ' There is',nadd_max,' additional thermodynamic quantities available.'
         write(*,*) ' It is defined in the data sheet eos.pdf on the website of the EoS.'
      else
         write(*,*) ' There are',nadd_max,' additional thermodynamic quantities available.'
         write(*,*) ' They are defined in the data sheet eos.pdf on the website of the EoS.'
      end if
      write(*,*)
      if (nadd_max == 1) then
         write(*,*) ' Do you want to include this additional quantity in the file eos.table?'
         write(*,*) ' Please select: 1 for yes, any other number for no'
         read(*,*) n_add
         if (n_add /= 1) n_add = 0
         idx_a(1) = n_add
      else
         n_add = -1
         do while ((n_add < 0).or.(n_add > nadd_max))
            write(*,*) ' How many additional thermodynamic quantities do you want to select for the file eos.table?'
            read(*,*) n_add
            if ((n_add < 0).or.(n_add > nadd_max)) then
               write(*,*) ' Number out of range!'
            end if
         end do
!2017/05/22
         if (n_add > 0) then
            write(*,*) ' Please select the additional thermodynamic quantities'
            write(*,*) ' in the order of your choice.'
            do iv=1,n_add,1
               idx_a(iv) = 0
               do while ((idx_a(iv) < 1).or.(idx_a(iv) > nadd_max))
                  write(*,*) 'Additional quantity #',iv,'?'
                  read(5,*) idx_a(iv)
                  if ((idx_a(iv) < 1).or.(idx_a(iv) > nadd_max)) then
                     write(*,*) ' Number out of range!'
                  end if
               end do
            end do
         end if
      end if
   else
      n_add = 0
   end if

!2017/05/22
   if (irpl == 0) then
      write(*,*)
      write(*,*) ' The following function values and derivatives'
      write(*,*) ' of the free energy per baryon are available:'
      write(*,*)
      write(*,*) ' quantity:         unit:    index:'
      write(*,*) ' F                 [MeV]       1'
      if (ivar(1) == 1) then
         write(*,*) ' dF/dT             []          2'
         write(*,*) ' d^2F/(dT^2)       [MeV^-1]    3'
         if (ivar(2) == 1) then
            write(*,*) ' d^2F/(dT dn_b)    [fm^3]      4'
         end if
         if (ivar(3) == 1) then
            write(*,*) ' d^2F/(dT dY_q)    []          5'
         end if
      endif
      if (ivar(2) == 1) then
         write(*,*) ' dF/dn_b           [MeV fm^3]  6'
         write(*,*) ' d^2F/(dn_b^2)     [MeV fm^6]  7'
         if (ivar(3) == 1) then
            write(*,*) ' d^2F/(dn_b dY_q)  [MeV fm^3]  8'
         end if
      endif
      if (ivar(3) == 1) then
         write(*,*) ' dF/dY_q           [MeV]       9'
         write(*,*) ' d^2F/(dY_q^2)     [MeV]      10'
      end if
   end if
   write(*,*)
   n_df = -1
   do while ((n_df < 0).or.(n_df > dim_df))
      write(*,*) ' How many quantities do you want to select for the file eos.table?'
      read(*,*) n_df
      if ((n_df < 0).or.(n_df > dim_df)) then
         write(*,*) ' Number out of range!'
      end if
   end do
   if (n_df > 0) then
      write(*,*) ' Please select the indices of the thermodynamic quantities'
      write(*,*) ' in the order of your choice.'
      do iv=1,n_df,1
         idx_df(iv) = 0
         do while ((idx_df(iv) < 1).or.(idx_df(iv) > dim_df))
            write(*,*) 'Index #',iv,'?'
            read(5,*) idx_df(iv)
            if ((idx_df(iv) < 1).or.(idx_df(iv) > dim_df)) then
               write(*,*) ' Number out of range!'
            end if
         end do
      end do
   end if

   if (np_max > 0) then
      write(*,*)
      write(*,*) ' There are composition data available on the following particle species:'
      write(*,*)
      write(*,*) ' particle number:  particle name:  particle index:'
      do iv=1,np_max,1
         if (idxp(iv) < 1000) then
            write(*,FMT=1000) iv,part(idxp(iv)),idxp(iv)
         else
            select case (idxp(iv))
            case (2001)
               write(*,FMT=1000) iv,'deuteron 2H',idxp(iv)
            case (3001)
               write(*,FMT=1000) iv,'triton 3H  ',idxp(iv)
            case (3002)
               write(*,FMT=1000) iv,'helion 3He ',idxp(iv)
            case (4002)
               write(*,FMT=1000) iv,'alpha 4He  ',idxp(iv)
            case default
               write(*,FMT=1000) iv,'nucleus    ',idxp(iv)
            end select
         end if
      end do
      write(*,*)
      n_p = -1
      do while ((n_p < 0).or.(n_p > np_max))
         write(*,*) ' How many particles do you want to select for the file eos.table?'
         read(*,*) n_p
         if ((n_p < 0).or.(n_p > np_max)) then
            write(*,*) ' Number out of range!'
         end if
      end do
      if (n_p > 0) then
         write(*,*) ' Please select the particle numbers'
         write(*,*) ' in the order of your choice.'
         do iv=1,n_p,1
            ip = 0
            do while ((ip < 1).or.(ip > np_max))
               write(*,*) ' Particle # ?'
               read(5,*) ip
               if ((ip < 1).or.(ip > np_max)) then
                  write(*,*) ' Number out of range!'
               end if
            end do
            idxpp(iv) = idxp(ip)
         end do
      end if
   end if

   if (nq_max > 0) then
      write(*,*)
      write(*,*) ' There are average mass, charge and neutron numbers'
      write(*,*) ' available for the following sets of particles'
      write(*,*) ' (see data sheet eos.pdf for definition)'
      write(*,*)
      write(*,*) ' set of particles     index'
      do iv=1,nq_max,1
         write(*,FMT=2000) iv,idxq(iv)
      end do
      write(*,*)
      n_q = -1
      do while ((n_q < 0).or.(n_q > nq_max))
         write(*,*) ' How many sets do you want to select for the file eos.table?'
         read(*,*) n_q
         if ((n_q < 0).or.(n_q > nq_max)) then
            write(*,*) ' Number out of range!'
         end if
      end do
      if (nq_max == 1) then
         idxqq(1) = idxq(1)
      else if (n_q > 0) then
         write(*,*) ' Please select the set numbers for the file eos.table'
         write(*,*) ' in the order of your choice.'
         do iv=1,n_q,1
            iq = 0
            do while ((iq < 1).or.(iq > nq_max))
               write(*,*) ' Set # ?'
               read(5,*) iq
               if ((iq < 1).or.(iq > nq_max)) then
                  write(*,*) ' Number out of range!'
               end if
            end do
            idxqq(iv) = idxq(iq)
         end do
      end if
   end if

   if (nm_max > 0) then
      write(*,*)
      write(*,*) ' There are microscopic data available of the following type:'
      write(*,*)
      write(*,*) ' quantity number: quantity type:                      of particle:  quantity index:'
      do iv=1,nm_max,1
         ip = INT(0.001*DBLE(idxm(iv))+0.001)
         iq = idxm(iv)-1000*ip
         write(*,FMT=3000) iv,micr(iq),'of ',part(ip),idxm(iv)
      end do
      write(*,*)
      n_m = -1
      do while ((n_m < 0).or.(n_m > nm_max))
         write(*,*) ' How many quantities do you want to select for the file eos.table?'
         read(*,*) n_m
         if ((n_m < 0).or.(n_m > nm_max)) then
            write(*,*) ' Number out of range!'
         end if
      end do
      if (n_m > 0) then
         write(*,*) ' Please select the quantity numbers for the file eos.table'
         write(*,*) ' in the order of your choice.'
         do iv=1,n_m,1
            im = 0
            do while ((im < 1).or.(im > nm_max))
               write(*,*) ' Quantity # ?'
               read(5,*) im
               if ((im < 1).or.(im > nm_max)) then
                  write(*,*) ' Number out of range!'
               end if
            end do
            idxmm(iv) = idxm(im)
         end do
      end if
   end if

   write(*,*)
   write(*,*) ' There are error estimates available of the following type:'
   write(*,*)
   write(*,*) ' quantity number: quantity type:'
   do iv=1,nerr_max,1
      write(*,FMT=4000) iv,qerr(iv)
   end do
   write(*,*)
   n_err = -1
   do while ((n_err < 0).or.(n_err > nerr_max))
      write(*,*) ' How many quantities do you want to select for the file eos.table?'
      read(*,*) n_err
      if ((n_err < 0).or.(n_err > nerr_max)) then
         write(*,*) ' Number out of range!'
      end if
   end do
   if (n_err > 0) then
      write(*,*) ' Please select the quantity numbers for the file eos.table'
      write(*,*) ' in the order of your choice.'
      do iv=1,n_err,1
         ierr = 0
         do while ((ierr < 1).or.(ierr > nerr_max))
            write(*,*) ' Quantity # ?'
            read(5,*) ierr
            if ((ierr < 1).or.(ierr > nerr_max)) then
               write(*,*) ' Number out of range!'
            end if
         end do
         idxerr(iv) = ierr
      end do
   end if

   iout = 1 ! Default is ASCII format
#if defined have_hdf5
   write(*,*)
   write(*,*) ' Please select the format of the file eos.table from'
   write(6,*) ' 1: ASCII, else: HDF5'
   read(5,*) iout
#endif

1000 format(7x,i4,11x,a11,6x,i6)
2000 format(7x,i4,9x,i6)
3000 format(7x,i4,8x,a33,3x,a3,a11,3x,i6)
4000 format(7x,i4,8x,a27)

!2017/05/22
   write(iunit3,*) '# number of regular, additional and derivative quantities (see table 7.1)'
   write(iunit3,*) n_qty,n_add,n_df
   write(iunit3,*) '# indices of regular, additional and derivative quantities'
   write(iunit3,*) (idx_qty(iv),iv=1,n_qty,1),(idx_a(iv),iv=1,n_add,1),(idx_df(iv),iv=1,n_df,1)
   write(iunit3,*) '# number of pairs and quadruples for composition data (see table 3.2/3.3)'
   write(iunit3,*) n_p,n_q
   write(iunit3,*) '# indices for pairs and quadruples for composition data'
   write(iunit3,*) (idxpp(iv),iv=1,n_p,1),(idxqq(iv),iv=1,n_q,1)
   write(iunit3,*) '# number of microscopic quantities (see section 4.2.4)'
   write(iunit3,*) n_m
   write(iunit3,*) '# indices of microscopic quantities'
   write(iunit3,*) (idxmm(iv),iv=1,n_m,1)
   write(iunit3,*) '# number of error quantities'
   write(iunit3,*) n_err
   write(iunit3,*) '# indices of error quantities'
   write(iunit3,*) (idxerr(iv),iv=1,n_err,1)
   write(iunit3,*) '# format of output file, ASCII (1) or HDF5 (else)'
   write(iunit3,*) iout

   close(unit=iunit3)

   write(*,*)
   write(*,*) ' new file eos.quantities generated'
   write(*,*)

end if

end SUBROUTINE init_quant
!***********************************************************************
SUBROUTINE init_para()
! Stefan Typel for the CompOSE core team, version 1.07, 2017/12/13
 use compose_internal
 use general_var, only: tabulation_schema
implicit none
integer :: iunit2,iunit4,ierror,iv,idim,ibeta,i_entr,inew,iflag,ibeta2,&
     npt(4),ivar(4),ipl(4),ll(4)
double precision :: x,y,z,zero,min(4),max(4)
character(18) :: cdim(3)

cdim(1) = 'one-dimensional.  '
cdim(2) = 'two-dimensional.  '
cdim(3) = 'three-dimensional.'

iunit2 = 31

zero = 0.d00

open(unit=iunit2,file='eos.init',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   read(iunit2,*) (ivar(iv),iv=1,4,1),irpl
   do iv=1,4,1
      read(iunit2,*) para_min(iv),para_max(iv)
   end do
   read(iunit2,*) incl_l,nadd_max,nall_max,ibeta2
else
   write(*,*)
   write(*,*) ' There is no file eos.init.'
   write(*,*) ' Please restart the program compose with task 1'
   write(*,*) ' and generate a new file eos.init.'
   write(*,*)
   stop
end if
close(unit=iunit2)

iunit4 = 33

open(unit=iunit4,file='eos.parameters',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   write(*,*)
   write(*,*) ' The file eos.parameters exists already.'
   write(*,*) ' Do you want to generate a new file eos.parameters?'
   write(*,*) ' Please select:'
   write(*,*) ' 1: Yes'
   write(*,*) ' else: No'
   read(5,*) inew
else
   inew = 1
end if
close(unit=iunit4)

if (inew == 1) then

   idim = 0
   do iv=1,4,1
      idim = idim+ivar(iv)
   end do

   write(*,*)
   write(*,FMT=1000) ' The EoS table is ',cdim(idim)
   write(*,*) ' It depends on'
   if (ivar(1) == 1) write(*,*) ' - the temperature T'
   if (ivar(2) == 1) write(*,*) ' - the baryon density n_b'
   if (ivar(3) == 1) write(*,*) ' - the hadronic charge fraction Y_q'
   if (ivar(4) == 1) write(*,*) ' - the magnetic field strength B'
   write(*,*)
   if (idim == 1) then
      write(*,*) ' The interpolation in this parameter can be'
   else
      write(*,*) ' The interpolation in these parameters can be'
   end if
   write(*,*) ' (1) continuous in the function values'
   write(*,*) ' (2) continuous in the function values and first derivatives'
   write(*,*) ' (3) continuous in the function values, the first and second derivatives'
   do iv=1,4,1
      ipl(iv) = 0
      if (ivar(iv) == 1) then
         write(*,*)
         write(*,*) ' Please select the interpolation order (1, 2, or 3) '
         if (iv == 1) write(*,*) ' for the temperature T'
         if (iv == 2) write(*,*) ' for the baryon density n_b'
         if (iv == 3) write(*,*) ' for the hadronic charge fraction Y_q'
         if (iv == 4) write(*,*) ' for the magnetic field strength B'
         do while ((ipl(iv) < 1).or.(ipl(iv) > 3))
            read(5,*) ipl(iv)
            if ((ipl(iv) < 1).or.(ipl(iv) > 3)) then
               write(*,*) ' Number out of range!'
            end if
         end do
      else
         ipl(iv) = 1
      end if
   end do

   if ((incl_l == 1).and.(ivar(3) == 1)) then
      write(*,*)
      write(*,*) ' Please select if you want to calculate the EoS of matter in beta-equilibrium.'
      write(*,*) ' 1: yes, else: no'
      read(5,*) ibeta
      if (ibeta /= 1) ibeta = 0
   else
      ibeta = 0
   end if

   i_entr = 0
   if (ivar(1) == 1) then
      write(*,*)
      write(*,*) ' Please select if you want the calculate the EoS for given entropy per baryon.'
      write(*,*) ' 1: yes, else: no'
      read(5,*) i_entr
      if (i_entr /= 1) i_entr = 0
   end if

   write(*,*)
   write(*,*) ' Please select the tabulation scheme for the parameters from'
   write(*,*) ' 0: explicit listing of parameter values'
   write(*,*) ' 1: loop form of parameter values'
   tabulation_schema = -1
   do while ((tabulation_schema < 0).or.(tabulation_schema > 1))
      read(5,*) tabulation_schema
      if ((tabulation_schema < 0).or.(tabulation_schema > 1)) then
         write(*,*) ' Number out of range!'
      end if
   end do

   open(unit=iunit4,file='eos.parameters',&
        status='unknown',action='write',iostat=ierror)
   write(iunit4,*) '# order of interpolation in first, second and third index'
   select case (irpl)
   case (0)
      write(iunit4,*) ipl(1),ipl(2),ipl(3)
   case (1)
      write(iunit4,*) ipl(4),ipl(2),ipl(3)
   case (2)
      write(iunit4,*) ipl(1),ipl(4),ipl(3)
   case (3)
      write(iunit4,*) ipl(1),ipl(2),ipl(4)
   end select
   write(iunit4,*) '# calculation of beta-equilibrium (1: yes, else: no) and for given entropy (1: yes, else: no)'
   write(iunit4,*) ibeta,i_entr
   write(iunit4,*) '# tabulation scheme (0 = explicit listing, 1 = loops, see manual)'
   write(iunit4,*) tabulation_schema
   write(iunit4,*) '# parameter values (first, second and third index) depending on tabulation scheme'

   if (tabulation_schema == 0) then
      write(*,*)
      write(*,*) ' How many data points of your EoS table do you want to generate?'
      npt(1) = -1
      do while (npt(1) < 0)
         read(5,*) npt(1)
         if (npt(1) < 0) then
            write(*,*) ' Number out of range!'
         end if
      end do
      write(iunit4,*) npt(1)
      do iv=1,npt(1),1
         if (ibeta == 0) then
            iflag = 0
            select case (irpl)
            case (0)
               do while (iflag == 0)
                  iflag = 1
                  select case (idim)
                  case (3)
                     if (i_entr == 1) then
                        write(*,FMT=2000) ' Please enter values for S, n_b, and Y_q for data point #',iv
                     else
                        write(*,FMT=2000) ' Please enter values for T, n_b, and Y_q for data point #',iv
                     end if
                     read(5,*) x,y,z
                  case (2)
                     if (ivar(1) == 0) then
                        write(*,FMT=2000) ' Please enter values for n_b and Y_q for data point #    ',iv
                        read(5,*) y,z
                     end if
                     if (ivar(2) == 0) then
                        if (i_entr == 1) then
                           write(*,FMT=2000) ' Please enter values for S and Y_q for data point #      ',iv
                        else
                           write(*,FMT=2000) ' Please enter values for T and Y_q for data point #      ',iv
                        endif
                        read(5,*) x,z
                     end if
                     if (ivar(3) == 0) then
                        if (i_entr == 1) then
                           write(*,FMT=2000) ' Please enter values for S and n_b for data point #      ',iv
                        else
                           write(*,FMT=2000) ' Please enter values for T and n_b for data point #      ',iv
                        end if
                        read(5,*) x,y
                     end if
                  case (1)
                     if (ivar(1) == 1) then
                        if (i_entr == 1) then
                           write(*,FMT=2000) ' Please enter values for S for data point #              ',iv
                        else
                           write(*,FMT=2000) ' Please enter values for T for data point #              ',iv
                        end if
                        read(5,*) x
                     end if
                     if (ivar(2) == 1) then
                        write(*,FMT=2000) ' Please enter values for n_b for data point #            ',iv
                        read(5,*) y
                     end if
                     if (ivar(3) == 1) then
                        write(*,FMT=2000) ' Please enter values for Y_q for data point #            ',iv
                        read(5,*) z
                     end if
                  end select

                  if (i_entr /= 1) then
                     if (ivar(1) == 1) then
                        if ((x < para_min(1)).or.(x > para_max(1))) then
                           write(*,*) ' Temperature T out of range'
                           iflag = 0
                        end if
                     else
                        x = para_min(1)
                     end if
                  end if
                  if (ivar(2) == 1) then
                     if ((y < para_min(2)).or.(y > para_max(2))) then
                        write(*,*) ' Baryon density n_b out of range'
                        iflag = 0
                     end if
                  else
                     y = para_min(2)
                  end if
                  if (ivar(3) == 1) then
                     if ((z < para_min(3)).or.(z > para_max(3))) then
                        write(*,*) ' Hadronic charge fraction Y_q out of range'
                        iflag = 0
                     end if
                  else
                     z = para_min(3)
                  end if
               end do
               write(iunit4,*) x,y,z
            case (1)
               do while (iflag == 0)
                  write(*,FMT=2000) ' Please enter values for B, n_b, and Y_q for data point #',iv
                  iflag = 1
                  read(5,*) x,y,z
                  if ((x < para_min(4)).or.(x > para_max(4))) then
                     write(*,*) ' Magnetic field strength B out of range'
                     iflag = 0
                  end if
                  if ((y < para_min(2)).or.(y > para_max(2))) then
                     write(*,*) ' Baryon density n_b out of range'
                     iflag = 0
                  end if
                  if ((z < para_min(3)).or.(z > para_max(3))) then
                     write(*,*) ' Hadronic charge fraction Y_q out of range'
                     iflag = 0
                  end if
               end do
               write(iunit4,*) zero,y,z,x
            case (2)
               do while (iflag == 0)
                  if (i_entr == 1) then
                     write(*,FMT=2000) ' Please enter values for S, B, and Y_q for data point #',iv
                  else
                     write(*,FMT=2000) ' Please enter values for T, B, and Y_q for data point #',iv
                  end if
                  iflag = 1
                  read(5,*) x,y,z
                  if (i_entr /= 1) then
                     if ((x < para_min(1)).or.(x > para_max(1))) then
                        write(*,*) ' Temperature T out of range'
                        iflag = 0
                     end if
                  end if
                  if ((y < para_min(4)).or.(y > para_max(4))) then
                     write(*,*) ' Magnetic field strength B out of range'
                     iflag = 0
                  end if
                  if ((z < para_min(3)).or.(z > para_max(3))) then
                     write(*,*) ' Hadronic charge fraction Y_q out of range'
                     iflag = 0
                  end if
               end do
               write(iunit4,*) x,zero,z,y
            case (3)
               do while (iflag == 0)
                  if (i_entr == 1) then
                     write(*,FMT=2000) ' Please enter values for S, n_b, and B for data point #',iv
                  else
                     write(*,FMT=2000) ' Please enter values for T, n_b, and B for data point #',iv
                  end if
                  iflag = 1
                  read(5,*) x,y,z
                  if (i_entr /= 1) then
                     if ((x < para_min(1)).or.(x > para_max(1))) then
                        write(*,*) ' Temperature T out of range'
                        iflag = 0
                     end if
                  end if
                  if ((y < para_min(2)).or.(y > para_max(2))) then
                     write(*,*) ' Baryon density n_b out of range'
                     iflag = 0
                  end if
                  if ((z < para_min(4)).or.(z > para_max(4))) then
                     write(*,*) ' Magnetic field strength B out of range'
                     iflag = 0
                  end if
               end do
               write(iunit4,*) x,y,zero,z
            end select
         else
            iflag = 0
            select case (irpl)
            case (0)
               do while (iflag == 0)
                  if (i_entr == 1) then
                     write(*,FMT=2000) ' Please enter values for S and n_b for data point #',iv
                  else
                     write(*,FMT=2000) ' Please enter values for T and n_b for data point #',iv
                  end if
                  iflag = 1
                  read(5,*) x,y
                  if (i_entr /= 1) then
                     if ((x < para_min(1)).or.(x > para_max(1))) then
                        write(*,*) ' Temperature T out of range'
                        iflag = 0
                     end if
                  end if
                  if ((y < para_min(2)).or.(y > para_max(2))) then
                     write(*,*) ' Baryon density n_b out of range'
                     iflag = 0
                  end if
               end do
               write(iunit4,*) x,y,zero
            case (1)
               do while (iflag == 0)
                  write(*,FMT=2000) ' Please enter values for B and n_b for data point #',iv
                  iflag = 1
                  read(5,*) x,y,z
                  if ((x < para_min(4)).or.(x > para_max(4))) then
                     write(*,*) ' Magnetic field strength B out of range'
                     iflag = 0
                  end if
                  if ((y < para_min(2)).or.(y > para_max(2))) then
                     write(*,*) ' Baryon density n_b out of range'
                     iflag = 0
                  end if
               end do
               write(iunit4,*) zero,y,zero,x
            case (2)
               do while (iflag == 0)
                  if (i_entr == 1) then
                     write(*,FMT=2000) ' Please enter values for S and B for data point #',iv
                  else
                     write(*,FMT=2000) ' Please enter values for T and B for data point #',iv
                  end if
                  iflag = 1
                  read(5,*) x,y,z
                  if ((x < para_min(1)).or.(x > para_max(1))) then
                     write(*,*) ' Temperature T out of range'
                     iflag = 0
                  end if
                  if ((y < para_min(4)).or.(y > para_max(4))) then
                     write(*,*) ' Magnetic field strength B out of range'
                     iflag = 0
                  end if
               end do
               write(iunit4,*) x,zero,zero,y
            end select
         end if

      end do
   else
      do iv=1,4,1
         if (ivar(iv) == 1) then
            if ((ibeta == 0).or.(iv /= 3)) then
               iflag = 0
               do while (iflag == 0)
                  write(*,*)
                  select case (iv)
                  case (1)
                     if (i_entr == 1) then
                        write(*,*) ' Please enter the minimum and maximum values for the entropy per baryon S'
                     else
                        write(*,*) ' Please enter the minimum and maximum values for the temperature T'
                     end if
                  case (2)
                     write(*,*) ' Please enter the minimum and maximum values for the baryon density n_b'
                  case (3)
                     write(*,*) ' Please enter the minimum and maximum values for the hadronic charge fraction Y_q'
                  case (4)
                     write(*,*) ' Please enter the minimum and maximum values for the magnetic field strength B'
                  end select
                  read(5,*) min(iv),max(iv)
                  iflag = 1
                  if ((i_entr /= 1).or.(iv /= 1)) then
                     if ((min(iv) < para_min(iv)).or.(min(iv) > para_max(iv))) then
                        write(*,*) ' Minimum out of range.'
                        iflag = 0
                     end if
                     if ((max(iv) < min(iv)).or.(max(iv) > para_max(iv))) then
                        write(*,*) ' Maximum out of range'
                        iflag = 0
                     end if
                  end if
               end do
               if (min(iv) < max(iv)) then
                  npt(iv) = 0
                  do while (npt(iv) < 2)
                     write(*,*) ' Please enter the number of grid points'
                     read(5,*) npt(iv)
                     if (npt(iv) < 2) write(*,*) ' Number out of range'
                  end do
                  write(*,*) ' Please select the scaling of the grid points from'
                  write(*,*) ' 0: linear'
                  write(*,*) ' else: logarithmic'
                  read(5,*) ll(iv)
                  if (ll(iv) /= 0) ll(iv) = 1
               else
                  npt(iv) = 1
                  ll(iv) = 0
               end if
            else
               min(3) = 0.5d00
               max(3) = 0.5d00
               npt(3) = 1
               ll(3) = 0
            end if
         else
            min(iv) = para_min(iv)
            max(iv) = para_max(iv)
            npt(iv) = 1
            ll(iv)  = 0
         end if
      end do

      if (irpl == 0) then
         write(iunit4,*) min(1),min(2),min(3)
         write(iunit4,*) max(1),max(2),max(3)
         write(iunit4,*) npt(1),npt(2),npt(3)
         write(iunit4,*) ll(1),ll(2),ll(3)
      else
         write(iunit4,*) min(1),min(2),min(3),min(4)
         write(iunit4,*) max(1),max(2),max(3),max(4)
         write(iunit4,*) npt(1),npt(2),npt(3),npt(4)
         write(iunit4,*) ll(1),ll(2),ll(3),ll(4)
      end if


   end if

1000 format(1x,a18,a18)
2000 format(1x,a57,i5)


   close(unit=iunit4)

   write(*,*)
   write(*,*) ' new file eos.parameters generated'
   write(*,*)

end if

end SUBROUTINE init_para
!***********************************************************************
SUBROUTINE get_eos_table_term(iwr,iinit)
! Stefan Typel for the CompOSE core team, version 1.02, 2017/11/16
 use compose_internal
 use omp_lib
 use m_get_tables
implicit none
integer :: iwr,iinit,nbl,iunit2,iunit3,iunit4,ierror,iv,ibeta,iyq,&
     ivar(4),idxp(dim_ip),idxq(dim_iq),idxm(dim_im),init_flg
double precision :: timei,timef


! maximum number of subtables
nbl = 10

iunit2 = 31
open(unit=iunit2,file='eos.init',&
     status='old',action='read',iostat=ierror)
if (ierror /= 0) then
   write(*,*)
   write(*,*) ' There is no file eos.init.'
   write(*,*) ' Please restart the program compose with task 1'
   write(*,*) ' and generate a new file eos.init.'
   write(*,*)
   stop
else
   read(iunit2,*) (ivar(iv),iv=1,4,1),irpl
   do iv=1,4
      read(iunit2,*) para_min(iv),para_max(iv)
   end do
   read(iunit2,*) incl_l,nadd_max,nall_max,ibeta
   read(iunit2,*) np_max,nq_max
   read(iunit2,*) (idxp(iv),iv=1,np_max,1)
   read(iunit2,*) (idxq(iv),iv=1,nq_max,1)
   read(iunit2,*) nm_max
   read(iunit2,*) (idxm(iv),iv=1,nm_max,1)
end if
close(unit=iunit2)

iunit3 = 32
open(unit=iunit3,file='eos.quantities',&
     status='old',action='read',iostat=ierror)
if (ierror /= 0) then
   write(*,*)
   write(*,*) ' There is no file eos.quantities.'
   write(*,*) ' Please restart the program compose with task 1'
   write(*,*) ' and generate a new file eos.quantities.'
   write(*,*)
   stop
end if
close(unit=iunit3)

iunit4 = 33
open(unit=iunit4,file='eos.parameters',&
     status='old',action='read',iostat=ierror)
if (ierror /= 0) then
   write(*,*)
   write(*,*) ' There is no file eos.parameters.'
   write(*,*) ' Please restart the program compose with task 2'
   write(*,*) ' and generate a new file eos.parameters.'
   write(*,*)
   stop
end if
close(unit=iunit4)

init_flg = iinit

call read_eos_4_tables(iwr,nbl,iyq,ii_thermo=2,ii_tynb=0,unit=0,iinit=init_flg)

call get_diff_rules()

call init_ipl_rule()

call get_eos_report(iwr)

call define_eos_table(iwr)

call get_eos_table(iwr)

end SUBROUTINE get_eos_table_term
!**********************************************************************
