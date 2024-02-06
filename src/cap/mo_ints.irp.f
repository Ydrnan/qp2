 BEGIN_PROVIDER [double precision, mo_wx_cap , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_wy_cap , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_wz_cap , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * W_x MO_j for CAP
 ! array of the integrals of MO_i * W_y MO_j for CAP
 ! array of the integrals of MO_i * W_z MO_j for CAP
 END_DOC
 implicit none

  call ao_to_mo(                                                     &
      ao_wx_cap,                                                   &
      size(ao_wx_cap,1),                                           &
      mo_wx_cap,                                                   &
      size(mo_wx_cap,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_wy_cap,                                                   &
      size(ao_wy_cap,1),                                           &
      mo_wy_cap,                                                   &
      size(mo_wy_cap,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_wz_cap,                                                   &
      size(ao_wz_cap,1),                                           &
      mo_wz_cap,                                                   &
      size(mo_wz_cap,1)                                            &
      )
END_PROVIDER

 BEGIN_PROVIDER [double precision, mo_spread_x_exact, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_spread_y_exact, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_spread_z_exact, (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * x^2 MO_j
 ! array of the integrals of MO_i * y^2 MO_j
 ! array of the integrals of MO_i * z^2 MO_j
 END_DOC
 implicit none
  call ao_to_mo(                                                     &
      ao_spread_x_exact,                                                   &
      size(ao_spread_x_exact,1),                                           &
      mo_spread_x_exact,                                                   &
      size(mo_spread_x_exact,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_spread_y_exact,                                                   &
      size(ao_spread_y_exact,1),                                           &
      mo_spread_y_exact,                                                   &
      size(mo_spread_y_exact,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_spread_z_exact,                                                   &
      size(ao_spread_z_exact,1),                                           &
      mo_spread_z_exact,                                                   &
      size(mo_spread_z_exact,1)                                            &
      )
END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_ints_cap_alpha,(mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_ints_cap_beta,(mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Cap integrals in mo basis.
  END_DOC

  integer :: i,j

  mo_ints_cap_alpha = mo_wx_cap + mo_wy_cap + mo_wz_cap
  mo_ints_cap_beta = mo_ints_cap_alpha

  mo_ints_cap_alpha(1:elec_alpha_num,:) = 0d0
  mo_ints_cap_alpha(:,1:elec_alpha_num) = 0d0
   
  mo_ints_cap_beta(1:elec_beta_num,:) = 0d0
  mo_ints_cap_beta(:,1:elec_beta_num) = 0d0

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_one_e_integrals_cap,(mo_num,mo_num)]
  implicit none
  BEGIN_DOC
 ! Cap integrals in mo basis.
  END_DOC

  mo_one_e_integrals_cap = mo_ints_cap_alpha + mo_ints_cap_beta
  !mo_one_e_integrals_cap = mo_wx_cap + mo_wy_cap + mo_wz_cap
  !mo_one_e_integrals_cap(1,1) = 0d0
  !mo_one_e_integrals_cap(1,:) = 0d0
  !mo_one_e_integrals_cap(:,1) = 0d0

  integer :: i,j
!  print*,'AO'
!  do i = 1, ao_num
!    !write(*,'(100(F12.6))') ,ao_spread_x_exact(i,:)
!    write(*,'(100(F12.6))') ,ao_wx_cap(i,:)
!  enddo
!  print*,''
!  do i = 1, ao_num
!    !write(*,'(100(F12.6))') ,ao_spread_y_exact(i,:)
!    write(*,'(100(F12.6))') ,ao_wy_cap(i,:)
!  enddo
!  print*,''
!  do i = 1, ao_num
!    !write(*,'(100(F12.6))') ,ao_spread_z_exact(i,:)
!    write(*,'(100(F12.6))') ,ao_wz_cap(i,:)
!  enddo
!
!  print*,'MO coef'
!  do i = 1, ao_num
!    write(*,'(100(F12.6))') ,mo_coef(i,:)
!  enddo
!
!  print*,'MO'
!  do i = 1, mo_num
!    write(*,'(100(F12.6))') ,mo_wx_cap(i,:)
!  enddo
!  print*,''
!  do i = 1, mo_num
!    write(*,'(100(F12.6))') ,mo_wy_cap(i,:)
!  enddo
!  print*,''
!  do i = 1, mo_num
!    write(*,'(100(F12.6))') ,mo_wz_cap(i,:)
!  enddo
!  
!  !PROVIDE ao_spread_x_exact ao_spread_y_exact ao_spread_z_exact  
!  double precision :: norm
!  norm = 0d0
!  do j = 1, ao_num
!    do i = 1, ao_num
!      norm += (eta_cap *(ao_wx_cap(i,j) + ao_wy_cap(i,j) + ao_wz_cap(i,j)))**2
!      !norm += (eta_cap *(ao_spread_x_exact(i,j) + ao_spread_y_exact(i,j) + ao_spread_z_exact(i,j)))**2
!      !norm += (eta_cap *(ao_spread_x(i,j) + ao_spread_y(i,j) + ao_spread_z(i,j)))**2
!    enddo
!  enddo
!  print*,'norm',dsqrt(norm)
!  
!
!  print*,'CAP int'
!  do i = 1, mo_num
!    write(*,'(100(F12.6))') mo_one_e_integrals_cap(i,:)
!  enddo
  !print*,nucl_coord
  !print*,'d',mo_spread_z +mo_spread_y+mo_spread_z

END_PROVIDER
