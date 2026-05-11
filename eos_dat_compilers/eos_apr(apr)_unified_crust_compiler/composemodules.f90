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
! modules
!***********************************************************************
!***********************************************************************

!***********************************************************************
module general_var
 ! Marco Mancini for the CompOSE core team,  2018/10/16
 ! contains common variables for compose

 implicit none
 private

 integer,parameter :: dp = kind(1.d0)

 type :: tabulation_var_i
   integer :: t = 0, nb = 0, yq = 0, b = 0, tnyb = 0
 end type tabulation_var_i
 type :: tabulation_var_r
   real(dp) :: t = 0._dp, nb = 0._dp, yq = 0._dp, b = 0._dp
 end type tabulation_var_r

 ! tabulation scheme is red here : 0 explicit , 1 loop
 integer,public :: tabulation_schema  ! replace i_tab
 type(tabulation_var_i),public :: pts
 type(tabulation_var_r),public :: tabMin, tabMax

end module general_var


!***********************************************************************
module eos_tables
 ! Stefan Typel for the CompOSE core team, version 1.11, 2017/05/22

 ! parameters
 ! dim_reg:   number of standard thermodynamic quantities in eos.thermo (7)
 ! dim_qtyt:  regular number of quantities in array eos_thermo (24)
 ! dim_qtye:  maximum number of error quantities (8)

 implicit none

 integer, parameter :: dim_reg=7, dim_qtyt=24, dim_qtye=8, dim_df=10

 integer :: nall_max, nadd_max, np_max, nq_max, nm_max

 integer, allocatable,dimension(:,:,:,:) :: idxp_compo, idxq_compo, idx_mic

 integer, allocatable, dimension(:) :: idx_compo_p, idx_compo_q, idx_micro

 double precision :: arg(3), eos_thermo(dim_qtyt), eos_thermo_err(dim_qtye),&
   &                 eos_df(dim_df)

 double precision, allocatable,dimension(:) ::  eos_thermo_add, eos_compo_p, eos_micro, eos_q

 double precision, allocatable,dimension(:,:) :: tab_para(:,:), eos_compo_q(:,:)


 double precision, allocatable,dimension(:,:,:,:) :: tab_thermo, tabp_compo, tabq_compo, tab_mic

 ! neutron mass in MeV/c^{2}, standard value
 double precision :: m_n = 939.565379d00

 ! proton mass in MeV/c^{2}, standard value
 double precision :: m_p = 938.272046d00

 real:: time0,time1
END MODULE eos_tables



!***********************************************************************
module compose_internal
 ! Stefan Typel for the CompOSE core team, version 1.14, 2017/05/22
 USE eos_tables
 implicit none

 integer, parameter :: dim_err = 100, dim_ip = 25, dim_iq = 10,&
   &                   dim_im = 100, dim_ia = 5, nerr_max = 8

 integer :: eos_dim, idx_ex(0:3), idx_ipl(3), &
   &        dim_a, dim_idx(3), min_idx(4), &
   &        idx_arg2(-4:5, -4:5, 0:2), &
   &        idx_arg1(-4:5, 0:2), imap(3), jmap(4), &
   &        ipl_rule(0:1, 3, -8:8), &
   &        error_msg(0:dim_err), &
   &        n_qty, n_add, idx_qty(dim_qtyt), idx_a(dim_ia), &
   &        n_df, idx_df(dim_df), &
   &        n_p, n_q, n_m, n_err, idx_err(dim_qtye), iout, &
   &        incl_l, inbyq, dim_p(4), irpl

 integer, allocatable,dimension(:) :: idx_add, idx_thermo, idx_p, idx_q, idx_m
 integer, allocatable,dimension(:,:,:,:) :: idx_arg, idx_argx

 double precision :: d1x(-4:5, -4:5, -4:4), &
   &                 d2x(-4:5, -4:5, -4:4), &
   &                 d1y(-4:5, -4:5, -4:4), &
   &                 d2y(-4:5, -4:5, -4:4), &
   &                 dx, dy, dx2, dy2, &
   ! &                 df(0:3, 0:3, -4:5, -4:5), &
   ! &                 fc(0:5, 0:5), &
   &                 val_rpl, &
   &                 para_min(4), &
   &                 para_max(4), &
   &                 arg2(4)

 double precision, allocatable, dimension(:,:) :: v_thermo
 double precision, allocatable, dimension(:,:,:,:) :: r1d, r2d




contains

 !***********************************************************************
 INTEGER FUNCTION get_ipl_rule(ir, idx, ipl)
  ! Stefan Typel for the CompOSE core team, version 0.02, 2016/06/20

  integer, intent(in) :: ir,idx,ipl

  if (idx > 0) then
    get_ipl_rule = ipl_rule(1,ipl,ir)
  else
    get_ipl_rule = ipl_rule(0,ipl,ir)
  end if

  return
 end FUNCTION get_ipl_rule


END MODULE compose_internal
!***********************************************************************
