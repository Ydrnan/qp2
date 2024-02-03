
subroutine diagonalize_ci_cap(u_in)

   implicit none
   complex*16, intent(inout) :: u_in(N_det,N_states)

   complex*16, allocatable :: ci_eigenvectors_cap(:,:)
   complex*16, allocatable :: ci_s2_cap(:)
   complex*16, allocatable :: ci_electronic_energy_cap(:)
   complex*16, allocatable :: h_tmp(:,:)
   
   integer                        :: i_good_state
   integer, allocatable           :: index_good_state_array(:)
   logical, allocatable           :: good_state_array(:)
   complex*16, allocatable  :: s2_values_tmp(:)
   integer                        :: i_other_state
   complex*16, allocatable  :: eigenvectors(:,:), eigenvalues(:), H_prime(:,:)
   integer                        :: i_state
   complex*16               :: e_0
   integer                        :: i,j,k, info
   complex*16, allocatable  :: s2_eigvalues(:)
   !complex*16, allocatable  :: e_array(:)
   integer, allocatable           :: iorder(:)
   logical                        :: converged
   logical                        :: do_csf

   PROVIDE threshold_davidson nthreads_davidson distributed_davidson

   allocate(ci_eigenvectors_cap(N_det,N_states_diag))
   allocate(ci_s2_cap(N_states_diag))
   allocate(ci_electronic_energy_cap(N_states_diag))

   ! Guess values for the "N_states" states of the |ci| eigenvectors
   do j=1,min(N_states,N_det)
     do i=1,N_det
       ci_eigenvectors_cap(i,j) = u_in(i,j)
     enddo
   enddo

   do j=min(N_states,N_det)+1,N_states_diag
     do i=1,N_det
       ci_eigenvectors_cap(i,j) = dcmplx(0.d0,0d0)
     enddo
   enddo

   if (diag_algorithm == "Davidson") then

     call davidson_diag_HS2_complex(psi_det,ci_eigenvectors_cap, ci_s2_cap, &
       size(ci_eigenvectors_cap,1),ci_electronic_energy_cap,               &
       N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,converged)

     integer :: N_states_diag_save
     N_states_diag_save = N_states_diag
     do while (.not.converged)
       complex*16, allocatable :: ci_electronic_energy_cap_tmp (:)
       complex*16, allocatable :: ci_eigenvectors_cap_tmp (:,:)
       complex*16, allocatable :: ci_s2_cap_tmp (:)

       N_states_diag  *= 2
       !TOUCH N_states_diag

       allocate (ci_electronic_energy_cap_tmp (N_states_diag) )
       allocate (ci_eigenvectors_cap_tmp (N_det,N_states_diag) )
       allocate (ci_s2_cap_tmp (N_states_diag) )

       ci_electronic_energy_cap_tmp(1:N_states_diag_save) = ci_electronic_energy_cap(1:N_states_diag_save)
       ci_eigenvectors_cap_tmp(1:N_det,1:N_states_diag_save) = ci_eigenvectors_cap(1:N_det,1:N_states_diag_save)
       ci_s2_cap_tmp(1:N_states_diag_save) = ci_s2_cap(1:N_states_diag_save)

       call davidson_diag_HS2_complex(psi_det,ci_eigenvectors_cap_tmp, ci_s2_cap_tmp, &
         size(ci_eigenvectors_cap_tmp,1),ci_electronic_energy_cap_tmp,               &
         N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,0,converged)

       ci_electronic_energy_cap(1:N_states_diag_save) = ci_electronic_energy_cap_tmp(1:N_states_diag_save)
       ci_eigenvectors_cap(1:N_det,1:N_states_diag_save) = ci_eigenvectors_cap_tmp(1:N_det,1:N_states_diag_save)
       ci_s2_cap(1:N_states_diag_save) = ci_s2_cap_tmp(1:N_states_diag_save)

       deallocate (ci_electronic_energy_cap_tmp)
       deallocate (ci_eigenvectors_cap_tmp)
       deallocate (ci_s2_cap_tmp)
     enddo

     if (N_states_diag > N_states_diag_save) then
       N_states_diag = N_states_diag_save
       deallocate(ci_eigenvectors_cap, ci_electronic_energy_cap, ci_s2_cap)
       allocate(ci_eigenvectors_cap(N_det,N_states_diag))
       allocate(ci_s2_cap(N_states_diag))
       allocate(ci_electronic_energy_cap(N_states_diag))
       
       !TOUCH N_states_diag
     endif

   else if (diag_algorithm == "Lapack") then

     print *,  'Diagonalization of H using Lapack'
     allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
     allocate (eigenvalues(N_det))
     if (s2_eig) then
       complex*16, parameter :: alpha = dcmplx(0.1d0,0d0)
       allocate (H_prime(N_det,N_det) )
       H_prime(1:N_det,1:N_det) = H_matrix_all_dets_complex(1:N_det,1:N_det) +  &
         alpha * dcmplx(S2_matrix_all_dets(1:N_det,1:N_det),0d0)
       do j=1,N_det
         H_prime(j,j) = H_prime(j,j) - alpha * dcmplx(expected_s2,0d0)
       enddo
       call diag_general_complex(eigenvalues,eigenvectors,H_prime,size(H_prime,1),N_det,info)
       call nullify_small_elements_complex(N_det,N_det,eigenvectors,size(eigenvectors,1),1.d-12)
       ci_electronic_energy_cap(:) = dcmplx(0.d0,0d0)
       i_state = 0
       allocate (s2_eigvalues(N_det))
       allocate(index_good_state_array(N_det),good_state_array(N_det))
       good_state_array = .False.
       call u_0_S2_u_0_complex(s2_eigvalues,eigenvectors,N_det,psi_det,N_int,&
         N_det,size(eigenvectors,1))
       if (only_expected_s2) then
         do j=1,N_det
           ! Select at least n_states states with S^2 values closed to "expected_s2"
           if(dabs(dble(s2_eigvalues(j))-expected_s2).le.0.5d0)then
             i_state +=1
             index_good_state_array(i_state) = j
             good_state_array(j) = .True.
           endif
           if(i_state.eq.N_states) then
             exit
           endif
         enddo
       else
         do j=1,N_det
           index_good_state_array(j) = j
           good_state_array(j) = .True.
         enddo
       endif
       if(i_state .ne.0)then
         ! Fill the first "i_state" states that have a correct S^2 value
         do j = 1, i_state
           do i=1,N_det
             ci_eigenvectors_cap(i,j) = eigenvectors(i,index_good_state_array(j))
           enddo
           ci_electronic_energy_cap(j) = eigenvalues(index_good_state_array(j))
           ci_s2_cap(j) = s2_eigvalues(index_good_state_array(j))
         enddo
         i_other_state = 0
         do j = 1, N_det
           if(good_state_array(j))cycle
           i_other_state +=1
           if(i_state+i_other_state.gt.n_states_diag)then
             exit
           endif
           do i=1,N_det
             ci_eigenvectors_cap(i,i_state+i_other_state) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_cap(i_state+i_other_state) = eigenvalues(j)
           ci_s2_cap(i_state+i_other_state) = s2_eigvalues(i_state+i_other_state)
         enddo

       else
         print*,''
         print*,'!!!!!!!!   WARNING  !!!!!!!!!'
         print*,'  Within the ',N_det,'determinants selected'
         print*,'  and the ',N_states_diag,'states requested'
         print*,'  We did not find only states with S^2 values close to ',expected_s2
         print*,'  We will then set the first N_states eigenvectors of the H matrix'
         print*,'  as the ci_eigenvectors_cap'
         print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the ci space'
         print*,''
         do j=1,min(N_states_diag,N_det)
           do i=1,N_det
             ci_eigenvectors_cap(i,j) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_cap(j) = eigenvalues(j)
           ci_s2_cap(j) = s2_eigvalues(j)
         enddo
       endif
       deallocate(index_good_state_array,good_state_array)
       deallocate(s2_eigvalues)
     else
       allocate(h_tmp(N_det,N_det))
       h_tmp = H_matrix_all_dets_complex
       call diag_general_complex(eigenvalues,eigenvectors,                    &
           h_tmp,size(h_tmp,1),N_det,info)
       deallocate(h_tmp)
       ci_electronic_energy_cap(:) = dcmplx(0.d0,0d0)
       call u_0_S2_u_0_complex(ci_s2_cap,eigenvectors,N_det,psi_det,N_int,       &
           min(N_det,N_states_diag),size(eigenvectors,1))
       ! Select the "N_states_diag" states of lowest energy
       do j=1,min(N_det,N_states_diag)
         do i=1,N_det
           ci_eigenvectors_cap(i,j) = eigenvectors(i,j)
         enddo
         ci_electronic_energy_cap(j) = eigenvalues(j)
     write(*,'(I8,2(F22.10))') j, ci_electronic_energy_cap(j) + dcmplx(nuclear_repulsion,0d0)
       enddo
     endif

     do k=1,N_states_diag
       ci_electronic_energy_cap(k) = dcmplx(0.d0,0d0)
       do j=1,N_det
         do i=1,N_det
           ci_electronic_energy_cap(k) +=                                &
               CONJG(ci_eigenvectors_cap(i,k)) * ci_eigenvectors_cap(j,k) *         &
               H_matrix_all_dets_complex(i,j)
         enddo
       enddo
     enddo
     deallocate(eigenvectors,eigenvalues)
   endif

   do k = 1, N_states
     do i = 1, N_det
       u_in(i,k) = ci_eigenvectors_cap(i,k)
     enddo
   enddo
   
   print*,'Energy of the states for eta =', eta_cap
   do k = 1, N_states
     write(*,'(I8,4(F22.10))') k, ci_electronic_energy_cap(k) + dcmplx(nuclear_repulsion,0d0), ci_s2_cap(k)
   enddo
   print*,''

end
