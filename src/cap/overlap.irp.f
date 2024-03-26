subroutine overlap_cap(psi_cap)

  implicit none

  complex*16, intent(in) :: psi_cap(N_det,N_states)
  double precision, allocatable :: overlap(:,:)
  double precision :: tmp
  integer :: states(N_states)
  complex*16 :: res

  integer :: i,j,k

  allocate(overlap(N_states,N_states))

  overlap = 0d0

  do j = 1, N_states
    do i = 1, N_states
      tmp = 0d0
      do k = 1, N_det
        tmp += psi_coef(k,j) * dble(psi_cap(k,i))
      enddo
      overlap(i,j) = dabs(tmp)
    enddo
  enddo 

  do i = 1, N_states
    states(i) = i
  enddo

  write(*,*) ''
  write(*,*) 'Overlap of the real part with the unperturbed wave function:'
  write(*,*) ''
  write(*,'(6X,100(I5,1X))') states(1:N_states)
  do i = 1, N_states
    write(*,'(I5,1X,100(F5.2,1X))') i, overlap(i,1:N_states)
  enddo
  write(*,*) ''

  deallocate(overlap)

end
