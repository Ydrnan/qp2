program diagonalize_h_cap

  implicit none

  complex*16, allocatable :: psi_cap(:,:)

  if (cap_onset_x >= 0d0 .and. cap_onset_y >= 0d0 .and. cap_onset_z >= 0d0) then
    if (eta_min >= 0d0 .and. eta_max > eta_min .and. eta_step_size > 0d0) then

      allocate(psi_cap(N_det,N_states))
      psi_cap = dcmplx(psi_coef,0d0)

      eta_cap = eta_min
      touch eta_cap
      do while (eta_cap <= eta_max)

        call diagonalize_ci_cap(psi_cap)
        
        eta_cap += eta_step_size
        touch eta_cap
      enddo

      deallocate(psi_cap)

    else if (eta_cap >= 0d0) then

      allocate(psi_cap(N_det,N_states))
      psi_cap = dcmplx(psi_coef,0d0)

      call diagonalize_ci_cap(psi_cap)
        
      deallocate(psi_cap)

    else
      print*,''
      print*,'Parameters for eta not setted correctly, no CAP calculation will be performed.'
      print*,''
      return  
    endif
  endif
end
