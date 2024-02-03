program test_dav_complex
  implicit none
  BEGIN_DOC
  ! To test the complex davidson routine
  END_DOC
  print *, 'Hello world'
  read_wf = .True.
  touch read_wf
  PROVIDE threshold_davidson nthreads_davidson
  call dav_complex
end

subroutine dav_complex
 implicit none
 complex*16, allocatable :: u_in(:,:), H_jj(:), energies(:), s2_out(:)
 integer :: sze,N_st,N_st_diag_in
 logical :: converged
 integer :: i,j,info
 double precision, allocatable :: u_in2(:,:), s2_out2(:), energies2(:)

 N_st = N_states
 N_st_diag_in = N_states_diag
 sze = N_det
 PROVIDE nuclear_repulsion

 !!! MARK THAT u_in mut dimensioned with "N_st_diag_in" as a second dimension 
 allocate(u_in2(sze,N_st_diag_in),energies2(N_st_diag_in),s2_out2(N_st_diag_in))
 allocate(u_in(sze,N_st_diag_in),H_jj(sze),energies(N_st_diag_in),s2_out(N_st_diag_in))

 u_in = 0.d0
 do i = 1, N_st
  u_in(1,i) = (1.d0,0d0)
 enddo
 do i = 1, sze
  H_jj(i) = H_matrix_all_dets_complex(i,i)
 enddo
 !call davidson_general_complex(u_in,H_jj,energies,sze,sze,N_st,N_st_diag_in,converged,H_matrix_all_dets_complex)
 !u_in2 = dble(u_in)
 !call davidson_diag_hs2(psi_det,u_in2,s2_out2,sze,energies2,sze,N_st,N_st_diag_in,N_int,0,converged)

 ! Lapack
 complex*16, allocatable :: eigvalues(:), eigvectors(:,:), h(:,:),h_cp(:,:)
 double precision, allocatable :: h_im(:,:), tmp(:,:)
 allocate(eigvalues(sze),eigvectors(sze,sze),h(sze,sze), h_im(sze,sze),h_cp(sze,sze))

 !do i = 1, sze
 !   write(*,'(100(F12.6))') H_matrix_all_dets_complex(i,:)
 !enddo
 call diag_general_complex(eigvalues,eigvectors,H_matrix_all_dets_complex,sze,sze,info)
 eigvalues += dcmplx(nuclear_repulsion,0d0)
 print*,'Exact energies:'
 do i = 1, min(10,N_det)
    write(*,'(I6,100(F18.10))') i, eigvalues(i)
 enddo

 !h = dcmplx(H_matrix_all_dets,0d0)
 !
 !mo_one_e_integrals = mo_one_e_integrals_cap
 !touch mo_one_e_integrals psi_det
 !PROVIDE H_matrix_all_dets
 !h_im = H_matrix_all_dets
 !mo_one_e_integrals = 0d0
 !touch mo_one_e_integrals psi_det
 !PROVIDE H_matrix_all_dets
 !h_im = h_im - H_matrix_all_dets
 !h = h + dcmplx(0d0, - eta_cap * h_im)
 !!do i = 1, sze
 !!   write(*,'(100(F12.6))') dble(h(i,:))
 !!enddo
 !! print*,'imag'
 !!do i = 1, sze
 !!   write(*,'(100(F12.6))') dimag(h(i,:))
 !!enddo

 !h_cp = h
 !call diag_general_complex(eigvalues,eigvectors,h,sze,sze,info)

 !eigvalues += dcmplx(nuclear_repulsion,0d0)
 !print*,'Exact energies:'
 !do i = 1, min(10,N_det)
 !   write(*,'(I6,100(F18.10))') i, eigvalues(i)
 !enddo

! do i = 1, sze
!   call print_det(psi_det(1,1,i),N_int)
!   write(*,'(100(F12.6))') eigvectors(i,1) / eigvectors(1,1)
!enddo
! call modified_gram_schmidt_c(eigvectors,sze,sze)
! do i = 1, sze
!   call print_det(psi_det(1,1,i),N_int)
!   write(*,'(100(F12.6))') eigvectors(i,1) / eigvectors(1,1)
!enddo
!
 !complex*16, allocatable :: rdm(:,:,:)
 !allocate(rdm(mo_num,mo_num,sze))
 !call modified_gram_schmidt_c(eigvectors,sze,sze)
 !call mo_one_rdm_cap(eigvectors, sze, sze, rdm)

 !complex*16, allocatable :: tmpc(:,:), w(:,:)
 ! complex*16 :: trace
 !allocate(tmpc(mo_num,mo_num),w(mo_num,mo_num))
 !w = dcmplx(0d0,mo_one_e_integrals_cap)
 !do j = 1, 8
 !  tmpc = matmul(rdm(:,:,j),w)
 !  !do i = 1, mo_num
 !  !  write(*,'(100(F12.6))') tmpc(i,:)
 !  !enddo 
 !  trace = (0d0,0d0)
 !  do i = 1, mo_num
 !    trace += tmpc(i,i)
 !  enddo
 !  write(*,'(I8,F12.6,F12.6)') j, eta_cap * trace
 !enddo

 !print*,'c-Energies:'
 !call modified_gram_schmidt_c(eigvectors,sze,sze)
 !call check_c_energy(h_cp,eigvectors,sze,sze)

 u_in = (0d0,0d0)
 do i = 1, N_st
  u_in(1,i) = (1.d0,0d0)
 enddo
 call davidson_diag_hs2_complex(psi_det,u_in,s2_out,sze,energies,sze,N_st,N_st_diag_in,N_int,converged)
 !u_in2 = dble(u_in)
 !call davidson_diag_hs2(psi_det,u_in2,s2_out2,sze,energies2,sze,N_st,N_st_diag_in,N_int,0,converged)

 !u_in2 = 0d0
 !do i = 1, N_st
 !  u_in2(1,i) = 1.d0
 !enddo
 !call davidson_diag_hs2(psi_det,u_in2,s2_out2,sze,energies2,sze,N_st,N_st_diag_in,N_int,0,converged)
 !call davidson_diag_hs2(psi_det,u_in2,s2_out2,sze,energies2,sze,N_st,N_st_diag_in,N_int,0,converged)

end
