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


MODULE precision

  integer,parameter :: dp = selected_real_kind(15)

END MODULE precision


MODULE hdfparameters

  USE precision

  real(dp), allocatable :: nb_hdf5(:),t_hdf5(:),y_q_hdf5(:)
  integer :: m_nb,n_temp,o_y_q

  real(dp), allocatable :: thermo_hdf5(:,:,:,:), thermo_hdf5_add(:,:,:,:),deriv_hdf5(:,:,:,:)
  real(dp), allocatable :: yi_hdf5(:,:,:,:),aav_hdf5(:,:,:,:),zav_hdf5(:,:,:,:),yav_hdf5(:,:,:,:)
  real(dp), allocatable :: nav_hdf5(:,:,:,:)
  real(dp), allocatable :: micro_hdf5(:,:,:,:),err_hdf5(:,:,:,:)
  integer,allocatable :: index_thermo(:),index_thermo_add(:),index_err(:),index_deriv(:)
  integer, allocatable :: index_yi(:),index_av(:),index_micro(:)

contains

  subroutine initialise_hdf5(n_nb,n_t,n_yq,i_beta,i_entr)

    use general_var, only : tabulation_schema
    USE compose_internal
    IMPLICIT NONE

    integer n_nb,n_t,n_yq,i_beta,i_entr

    m_nb = n_nb
    n_temp = n_t
    o_y_q = n_yq

! memory allocation for the data arrays, thermo

    IF(n_qty.ne.0) then
       allocate(thermo_hdf5(m_nb,n_temp,o_y_q,n_qty))
       thermo_hdf5 = 0._dp
       allocate(index_thermo(n_qty))
       index_thermo = 0
    else
       write(*,*) 'n_qty equal to zero, no thermodynamical quantities'
    end IF
    IF(n_add.ne.0) then
       allocate(thermo_hdf5_add(m_nb,n_temp,o_y_q,n_add))
       thermo_hdf5_add = 0._dp
       allocate(index_thermo_add(n_add))
       index_thermo_add = 0
    else
       write(*,*) 'n_add equal to zero, no additional thermodynamic quantities'
    end IF

    IF(n_df.ne.0) then
       allocate(deriv_hdf5(m_nb,n_temp,o_y_q,n_df))
       deriv_hdf5 = 0._dp
       allocate(index_deriv(n_df))
       index_deriv = 0
    else
       write(*,*) 'n_df equal to zero, no derivative quantities'
    end IF

! memory allocation for the data arrays, composition

    IF(n_p.ne.0) then
       allocate(yi_hdf5(m_nb,n_temp,o_y_q,n_p))
       allocate(index_yi(n_p))
       yi_hdf5 = 0._dp
       index_yi = 0
    else
       write(*,*) 'n_p equal to zero, no pairs'
    end if
    IF(n_q.ne.0) Then
       allocate(yav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(aav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(zav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(nav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(index_av(n_q))
       yav_hdf5 = 0._dp
       aav_hdf5 = 0._dp
       zav_hdf5 = 0._dp
       nav_hdf5 = 0._dp
       index_av = 0
    else
       write(*,*) 'n_q equal to zero, no quadruples'
    end if

! memory allocation for the data arrays, microscopic
    IF(n_m.ne.0) THEN
       allocate(micro_hdf5(m_nb,n_temp,o_y_q,n_m))
       allocate(index_micro(n_m))
       micro_hdf5 = 0._dp
       index_micro = 0
    else
       write(*,*) 'n_m equal to zero, no microscopic quantities'
    end if

    IF(n_err.ne.0) THEN
       allocate(err_hdf5(m_nb,n_temp,o_y_q,n_err))
       allocate(index_err(n_err))
       err_hdf5 = 0._dp
       index_err = 0
    else
       write(*,*) 'n_err equal to zero, no error quantities'
    end if


    IF(tabulation_schema.eq.0) then
       allocate(nb_hdf5(m_nb))
       allocate(t_hdf5(m_nb))
       allocate(y_q_hdf5(m_nb))
    else
       allocate(nb_hdf5(m_nb))
       if(i_entr == 0) then
          allocate(t_hdf5(n_temp))
       else
          allocate(t_hdf5(m_nb*n_temp))
          ! if we fix entropy per baryon, T is given
       end if
       if(i_beta==0) then
          allocate(y_q_hdf5(o_y_q))
       else
          allocate(y_q_hdf5(m_nb)) ! for beta_equilibrim, Y_e is determined
       end if
    end IF
    nb_hdf5 = 0._dp
    t_hdf5 = 0._dp
    y_q_hdf5 = 0._dp
  end subroutine initialise_hdf5

  subroutine write_hdf5(i_beta,i_entr)

    use general_var, only : tabulation_schema
    USE compose_internal
    use hdf5
    use modhdf5

    IMPLICIT NONE

    integer :: i_beta,i_entr
    
    integer n4d(4)
    integer(hid_t) ::  h5id, h5file_write
    character(LEN=15) :: file_write = 'eoscompose.h5'

    call hdf5_init()
    call hdf5_create_file(file_write,h5file_write)

    call hdf5_create_group(h5file_write,'Parameters', h5id)
    call hdf5_write_attr(h5id,'pointsnb',m_nb)
    call hdf5_write_attr(h5id,'tabulation_scheme',tabulation_schema)
    if(i_entr == 0) then
       call hdf5_write_attr(h5id,'pointst',n_temp)
    else
       call hdf5_write_attr(h5id,'pointst',n_temp*m_nb)
    end if
    if(i_beta == 0 ) then
       call hdf5_write_attr(h5id,'pointsyq',o_y_q)
    else
       call hdf5_write_attr(h5id,'pointsyq',m_nb)
    end if
    call hdf5_write_data(h5id,'nb',m_nb,nb_hdf5)
    if(tabulation_schema.eq.0) then
       call hdf5_write_data(h5id,'t',m_nb,t_hdf5)
       call hdf5_write_data(h5id,'yq',m_nb,y_q_hdf5)
    else
       if(i_entr == 0 ) then
          call hdf5_write_data(h5id,'t',n_temp,t_hdf5)
       else
          call hdf5_write_data(h5id,'t',n_temp*m_nb,t_hdf5)
       end if
       if(i_beta == 0) then
          call hdf5_write_data(h5id,'yq',o_y_q,y_q_hdf5)
       else
          call hdf5_write_data(h5id,'yq',m_nb,y_q_hdf5)
       end if
    end if
    call hdf5_close_group(h5id)
    n4d(1) = m_nb
    n4d(2) = n_temp
    n4d(3) = o_y_q

    write(*,*) 'writing ',n_qty,' thermodynamic quantities into file '

    if(n_qty.ne.0) then ! write thermo quantities
       call hdf5_create_group(h5file_write,'Thermo_qty', h5id)
       call hdf5_write_attr(h5id,'pointsqty',n_qty)
       n4d(4) = n_qty
       call hdf5_write_data(h5id,'thermo',n4d,thermo_hdf5)
       call hdf5_write_data(h5id,'index_thermo',n_qty,index_thermo)
       call hdf5_close_group(h5id)
    end if

    write(*,*) 'writing ',n_add,' additional thermodynamic quantities into file'
    if(n_add.ne.0) then
       call hdf5_create_group(h5file_write,'Thermo_add', h5id)
       call hdf5_write_attr(h5id,'pointsadd',n_add)
       n4d(4) = n_add
       call hdf5_write_data(h5id,'thermo_add',n4d,thermo_hdf5_add)
       call hdf5_write_data(h5id,'index_thermo_add',n_add,index_thermo_add)
       call hdf5_close_group(h5id)
    end if


    write(*,*) 'writing ',n_df,' derivative quantities into file'
    if(n_df.ne.0) then ! write derivative quantities
       call hdf5_create_group(h5file_write,'deriv_qty', h5id)
       call hdf5_write_attr(h5id,'pointsderiv',n_df)
       n4d(4) = n_df
       call hdf5_write_data(h5id,'deriv',n4d,deriv_hdf5)
       call hdf5_write_data(h5id,'index_deriv',n_df,index_deriv)
       call hdf5_close_group(h5id)
    end if


    write(*,*) 'writing ',n_p,' pairs into file'

    if(n_p.ne.0) then
       call hdf5_create_group(h5file_write,'Composition_pairs', h5id)
       call hdf5_write_attr(h5id,'pointspairs',n_p)
       n4d(4) = n_p
       call hdf5_write_data(h5id,'yi',n4d,yi_hdf5)
       call hdf5_write_data(h5id,'index_yi',n_p,index_yi)
       call hdf5_close_group(h5id)
    end if

    write(*,*) 'writing ',n_q,' quadruples into file'

    if(n_q.ne.0) then
       call hdf5_create_group(h5file_write,'Composition_quadrupels', h5id)
       call hdf5_write_attr(h5id,'pointsav',n_q)
       n4d(4) = n_q
       call hdf5_write_data(h5id,'yav',n4d,yav_hdf5)
       call hdf5_write_data(h5id,'aav',n4d,aav_hdf5)
       call hdf5_write_data(h5id,'zav',n4d,zav_hdf5)
       call hdf5_write_data(h5id,'nav',n4d,nav_hdf5)
       call hdf5_write_data(h5id,'index_av',n_q,index_av)
       call hdf5_close_group(h5id)
    end if


    write(*,*) 'writing ',n_m,' microscopic quantities into file'

    if(n_m.ne.0) then
       call hdf5_create_group(h5file_write,'Micro_qty', h5id)
       call hdf5_write_attr(h5id,'pointsmicro',n_m)
       n4d(4) = n_m
       call hdf5_write_data(h5id,'micro',n4d,micro_hdf5)
       call hdf5_write_data(h5id,'index_micro',n_m,index_micro)
       call hdf5_close_group(h5id)
    end if


    write(*,*) 'writing ',n_err,' error quantities into file'

    if(n_err.ne.0) then
       call hdf5_create_group(h5file_write,'Error_qty', h5id)
       call hdf5_write_attr(h5id,'pointserr',n_err)
       n4d(4) = n_err
       call hdf5_write_data(h5id,'error',n4d,err_hdf5)
       call hdf5_write_data(h5id,'index_err',n_err,index_err)
       call hdf5_close_group(h5id)
    end if

    call hdf5_close_file(h5file_write)


  end subroutine write_hdf5


  subroutine close_hdf5

    IMPLICIT NONE

    IF(allocated(thermo_hdf5)) deallocate(thermo_hdf5)
    IF(allocated(thermo_hdf5_add)) deallocate(thermo_hdf5_add)
    IF(allocated(deriv_hdf5)) deallocate(deriv_hdf5)
    IF(allocated(yi_hdf5)) deallocate(yi_hdf5)
    IF(allocated(yav_hdf5)) deallocate(yav_hdf5)
    IF(allocated(zav_hdf5)) deallocate(zav_hdf5)
    IF(allocated(nav_hdf5)) deallocate(nav_hdf5)
    IF(allocated(aav_hdf5)) deallocate(aav_hdf5)
    IF(allocated(micro_hdf5)) deallocate(micro_hdf5)
    IF(allocated(err_hdf5)) deallocate(err_hdf5)
    IF(allocated(index_thermo)) deallocate(index_thermo)
    IF(allocated(index_thermo_add)) deallocate(index_thermo_add)
    IF(allocated(index_deriv)) deallocate(index_deriv)
    IF(allocated(index_err)) deallocate(index_err)
    IF(allocated(index_yi)) deallocate(index_yi)
    IF(allocated(index_av)) deallocate(index_av)
    IF(allocated(index_micro)) deallocate(index_micro)
    IF(allocated(nb_hdf5)) deallocate(nb_hdf5)
    IF(allocated(t_hdf5)) deallocate(t_hdf5)
    IF(allocated(y_q_hdf5)) deallocate(y_q_hdf5)

  end subroutine close_hdf5

end MODULE hdfparameters
