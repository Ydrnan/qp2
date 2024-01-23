program test_dav_complex
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
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
 N_st = N_states
 N_st_diag_in = N_states_diag
 sze = N_det
 !!!! MARK THAT u_in mut dimensioned with "N_st_diag_in" as a second dimension 
 allocate(u_in(sze,N_st_diag_in),H_jj(sze),energies(N_st_diag_in),s2_out(N_st_diag_in))
 u_in = 0.d0
 do i = 1, N_st
  u_in(1,i) = (1.d0,0d0)
 enddo
 do i = 1, sze
  H_jj(i) = H_matrix_all_dets_complex(i,i) + nuclear_repulsion
 enddo
 call davidson_general_complex(u_in,H_jj,energies,sze,sze,N_st,N_st_diag_in,converged,H_matrix_all_dets_complex)

 complex*16, allocatable :: eigvalues(:), eigvectors(:,:)
 allocate(eigvalues(sze),eigvectors(sze,sze))
 call diag_general_complex(eigvalues,eigvectors,H_matrix_all_dets_complex,sze,sze,info)

 eigvalues += dcmplx(nuclear_repulsion,0d0)
 print*,'Exact energies:'
 do i = 1, min(10,N_det)
    write(*,'(I6,100(F18.10))') i, eigvalues(i)
 enddo

 u_in = 0.d0
 do i = 1, N_st
  u_in(1,i) = (1.d0,0d0)
 enddo
 call davidson_diag_hs2_complex(psi_det,u_in,s2_out,sze,energies,sze,N_st,N_st_diag_in,N_int,converged)

end
