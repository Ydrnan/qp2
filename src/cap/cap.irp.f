subroutine cap()

  implicit none

  complex*16, allocatable :: psi_cap(:,:), energy(:), l_energy(:,:)
  double precision :: write_e(2,N_states)
  integer :: i,j,n_eta

  if (cap_onset_x >= 0d0 .and. cap_onset_y >= 0d0 .and. cap_onset_z >= 0d0) then
    if (eta_min >= 0d0 .and. eta_max > eta_min .and. eta_step_size > 0d0) then

      n_eta = 0
      eta_cap = eta_min
      do while (eta_cap <= eta_max)
        eta_cap += eta_step_size
        n_eta += 1
      enddo

      allocate(psi_cap(N_det,N_states),energy(N_states), l_energy(N_states,n_eta))

      psi_cap = dcmplx(psi_coef,0d0)

      eta_cap = eta_min
      touch eta_cap
      do i = 1, n_eta

        call diagonalize_ci_cap(psi_cap,energy)
        l_energy(1:N_states,i) = energy(1:N_states)

        eta_cap += eta_step_size
        touch eta_cap
      enddo

      write(*,*) ''
      write(*,'(A,I8)') ' Number of eta values: ',n_eta
      write(*,*) ''
      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      write(*,'(A15)',advance='no') ' '
      do i = 1, N_states-1
        write(*,'(A17,I4,A12)',advance='no') 'State ',i,' '
      enddo
      write(*,'(A17,I4,A12)') '          State ',N_states,'     '

      !write(*,'(A15)',advance='no') '============== '
      write(*,'(A15)',advance='no') '   Eta      '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      write(*,'(A15)',advance='no') ' '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') '   Energy (Re)     Energy (Im)  '
      enddo
      write(*,'(A33)') '   Energy (Re)     Energy (Im)  '


      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      eta_cap = eta_min
      do i = 1, n_eta
        do j = 1, N_states
          write_e(1,j) = dble(l_energy(j,i))
          write_e(2,j) = dimag(l_energy(j,i))
        enddo
        write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') eta_cap, write_e(1:2,1:N_states)
        eta_cap += eta_step_size
      enddo

      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '
      write(*,*) ''
        
      deallocate(psi_cap,energy,l_energy)

    else if (eta_cap >= 0d0) then

      allocate(psi_cap(N_det,N_states),energy(N_states))
      psi_cap = dcmplx(psi_coef,0d0)

      call diagonalize_ci_cap(psi_cap,energy)
        
      deallocate(psi_cap)

    else
      print*,''
      print*,'Parameters for eta not setted correctly, no CAP calculation will be performed.'
      print*,''
      return  
    endif
  endif
end
