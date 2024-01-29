subroutine davidson_diag_hs2_complex(dets_in,u_in,s2_out,dim_in,energies,sze,N_st,N_st_diag,Nint,converged)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization.
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  complex*16, intent(inout) :: u_in(dim_in,N_st_diag)
  complex*16, intent(out)  :: energies(N_st_diag), s2_out(N_st_diag)
  logical, intent(out)           :: converged
  complex*16, allocatable  :: H_jj(:)

  double precision, external     :: diag_H_mat_elem, diag_S_mat_elem
  integer                        :: i,k,l
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_two_e_integrals_in_map
  allocate(H_jj(sze))

  H_jj(1) = dcmplx(diag_h_mat_elem(dets_in(1,1,1),Nint), 0d0)
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(sze,H_jj, dets_in,Nint)                    &
      !$OMP  PRIVATE(i)
  !$OMP DO SCHEDULE(static)
  do i=2,sze
    H_jj(i)  = dcmplx(diag_H_mat_elem(dets_in(1,1,i),Nint), 0d0)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call davidson_diag_hjj_sjj_complex(dets_in,u_in,H_jj,S2_out,energies,dim_in,sze,N_st,N_st_diag,Nint,converged)
  deallocate (H_jj)
end


subroutine davidson_diag_hjj_sjj_complex(dets_in,u_in,H_jj,s2_out,energies,dim_in,sze,N_st,N_st_diag_in,Nint,converged)
  use bitmasks
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! S2_out : Output : s^2
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! N_st_diag_in : Number of states in which H is diagonalized. Assumed > sze
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)           :: dim_in, sze, N_st, N_st_diag_in, Nint
  integer(bit_kind), intent(in) :: dets_in(Nint,2,sze)
  complex*16,  intent(in)       :: H_jj(sze)
  complex*16,  intent(inout)    :: s2_out(N_st_diag_in)
  complex*16, intent(inout)     :: u_in(dim_in,N_st_diag_in)
  complex*16, intent(out)       :: energies(N_st_diag_in)

  integer                       :: iter, N_st_diag
  integer                       :: i,j,k,l,m
  logical, intent(inout)        :: converged

  integer                       :: k_pairs, kl

  integer                       :: iter2, itertot
  complex*16, allocatable       :: y(:,:), h(:,:), h_p(:,:), lambda(:), s2(:), h_cp(:,:)
  complex*8, allocatable        :: y_s(:,:)
  complex*16, allocatable       :: s_(:,:), s_tmp(:,:), y_left(:,:), s_cp(:,:)
  double precision              :: diag_h_mat_elem
  complex*16, allocatable       :: residual_norm(:)
  character*(16384)             :: write_buffer
  complex*16                    :: to_print(3,N_st)
  double precision              :: cpu, wall
  integer                       :: shift, shift2, itermax, istate
  double precision              :: r1, r2, alpha, t1, t2
  logical                       :: state_ok(N_st_diag_in*davidson_sze_max)
  integer                       :: nproc_target, info
  integer                       :: order(N_st_diag_in)
  double precision              :: cmax, val
  complex*16, allocatable       :: U(:,:), overlap(:,:), S_d(:,:)
  double precision, allocatable :: max_U(:)
  integer, allocatable          :: pos_U(:)
  complex*16, pointer           :: W(:,:)
  complex*8, pointer            :: S(:,:)
  logical                       :: disk_based
  complex*16                    :: energy_shift(N_st_diag_in*davidson_sze_max), res

  include 'constants.include.F'

  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, S, y, y_s, S_d, h, lambda
  if (N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag*3
    stop -1
  endif

  itermax = max(2,min(davidson_sze_max, sze/N_st_diag))+1
  itertot = 0

  if (state_following) then
    allocate(overlap(N_st_diag*itermax, N_st_diag*itermax))
  else
    allocate(overlap(1,1))  ! avoid 'if' for deallocate
  endif
  overlap = (0.d0,0d0)

  PROVIDE nuclear_repulsion expected_s2 psi_bilinear_matrix_order psi_bilinear_matrix_order_reverse threshold_davidson_pt2 threshold_davidson_from_pt2

  call write_time(6)
  write(6,'(A)') ''
  write(6,'(A)') 'Davidson Diagonalization'
  write(6,'(A)') '------------------------'
  write(6,'(A)') ''

  ! Find max number of cores to fit in memory
  ! -----------------------------------------

  nproc_target = nproc
  double precision :: rss
  integer :: maxab
  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 2d0 * 8.d0 *                                   &! complex bytes
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.5d0*dble(sze*m)*(N_st_diag*itermax)  &! W,S
         + 1.d0*dble(sze)*(N_st_diag)             &! S_d
         + 4.5d0*(N_st_diag*itermax)**2           &! h,y,y_s,s_,s_tmp
         + 2.d0*(N_st_diag*itermax)               &! s2,lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
         + 3d0*(N_st_diag*itermax)**2             &! h_cp,y_left,s_cp
         + 2d0*(N_st_diag*itermax)                &! max_U,pos_U
                                                   ! In H_S2_u_0_nstates_zmq
         + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on collector
         + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on slave
         + 0.5d0*maxab                            &! idx0 in H_S2_u_0_nstates_openmp_work_*
         + nproc_target *                         &! In OMP section
           ( 1.d0*(N_int*maxab)                   &! buffer
           + 3.5d0*(maxab) )                      &! singles_a, singles_b, doubles, idx
         ) / 1024.d0**3


    if (nproc_target == 0) then
      call check_mem(r1,irp_here)
      nproc_target = 1
      exit
    endif

    if (r1+rss < qp_max_mem) then
      exit
    endif

    if (itermax > 4) then
      itermax = itermax - 1
    else if (m==1.and.disk_based_davidson) then
      m=0
      disk_based = .True.
      itermax = 6
    else
      nproc_target = nproc_target - 1
    endif

  enddo
  nthreads_davidson = nproc_target
  TOUCH nthreads_davidson
  call write_int(6,N_st,'Number of states')
  call write_int(6,N_st_diag,'Number of states in diagonalization')
  call write_int(6,sze,'Number of determinants')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================ =========== =========== =========== ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = 'Iter'
  do i=1,N_st
    if (i==1) then
    write_buffer = trim(write_buffer)//'       Energy                           S^2                   Residual             '
    else
    write_buffer = trim(write_buffer)//'                   Energy                           S^2                   Residual             '
    endif
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================ =========== =========== =========== ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)


  !if (disk_based) then
  !  ! Create memory-mapped files for W and S
  !  type(c_ptr) :: ptr_w, ptr_s
  !  integer :: fd_s, fd_w
  !  call mmap(trim(ezfio_work_dir)//'davidson_w', (/int(sze,8),int(N_st_diag*itermax,8)/),&
  !      8, fd_w, .False., ptr_w)
  !  call mmap(trim(ezfio_work_dir)//'davidson_s', (/int(sze,8),int(N_st_diag*itermax,8)/),&
  !      4, fd_s, .False., ptr_s)
  !  call c_f_pointer(ptr_w, w, (/sze,N_st_diag*itermax/))
  !  call c_f_pointer(ptr_s, s, (/sze,N_st_diag*itermax/))
  !else
    allocate(W(sze,N_st_diag*itermax), S(sze,N_st_diag*itermax))
  !endif

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax), S_d(sze,N_st_diag),                  &

      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
!      h_p(N_st_diag*itermax,N_st_diag*itermax),                      &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      s_(N_st_diag*itermax,N_st_diag*itermax),                       &
      s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
      residual_norm(N_st_diag),                                      &
      s2(N_st_diag*itermax),                                         &
      y_s(N_st_diag*itermax,N_st_diag*itermax),                      &
      lambda(N_st_diag*itermax),                                     &
      h_cp(N_st_diag*itermax,N_st_diag*itermax),                     &
      y_left(N_st_diag*itermax,N_st_diag*itermax),                   &
      s_cp(N_st_diag*itermax,N_st_diag*itermax),                     &
      max_U(N_st_diag*itermax), pos_U(N_st_diag*itermax))

  h = (0.d0,0d0)
  U = (0.d0,0d0)
  y = (0.d0,0d0)
  s_ = (0.d0,0d0)
  s_tmp = (0.d0,0d0)
  pos_U = 0
  max_U = 0d0

  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)

  ! Davidson iterations
  ! ===================

  converged = .False.

  ! Guess
  do k=N_st+1,N_st_diag
    do i=1,sze
        call random_number(r1)
        call random_number(r2)
        r1 = dsqrt(-2.d0*dlog(r1))
        r2 = dtwo_pi*r2
        u_in(i,k) = dcmplx(r1*dcos(r2),0d0) * u_in(i,k-N_st)
    enddo
    u_in(k,k) = u_in(k,k) + (10.d0,0d0)
  enddo
  do k=1,N_st_diag
    call normalize_complex(u_in(1,k),sze)
  enddo
  
  do k=1,N_st_diag
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
  enddo

  do while (.not.converged)
    itertot = itertot+1
    if (itertot == 8) then
      exit
    endif

    iter = 0
    do while (iter < itermax-1)
      iter += 1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter

        ! Orthogonalization of the guess vectors
        call ortho_qr_complex(U,size(U,1),sze,shift2)
        call ortho_qr_complex(U,size(U,1),sze,shift2)

        ! Change the sign of the guess vectors to match with the ones of the previous 
        ! iterations
        if (iter > 1) then
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
          do j = 1, shift
            if (sign(1d0,max_U(j)) * sign(1d0,dble(U(pos_U(j),j))) < -1d-30) then
              do i = 1, sze
                U(i,j) = - U(i,j)
              enddo
              !print*,'sign change!'
            endif
          enddo
          !$OMP END PARALLEL DO
        endif

        ! Sign of each vector
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,val)
        do j = shift+1, shift2
          val = 0d0
          do i = 1, sze
            if (dabs(dble(U(i,j))) > dabs(val)) then
              val = dble(U(i,j))
              max_U(j) = dble(U(i,j))
              pos_U(j) = i
            endif
          enddo
        enddo
        !$OMP END PARALLEL DO

        ! Computes HU
        call H_S2_u_0_nstates_openmp_complex(W(1,shift+1),S_d,U(1,shift+1),N_st_diag,sze)

        S(1:sze,shift+1:shift+N_st_diag) = cmplx(real(dble(S_d(1:sze,1:N_st_diag))), real(dimag(S_d(1:sze,1:N_st_diag))))

      ! Compute s_kl = <u_k | S_l> = <u_k| S2 |u_l>
      ! -------------------------------------------

       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) COLLAPSE(2)
       do j=1,shift2
         do i=1,shift2
           s_(i,j) = 0.d0
           do k=1,sze
             s_(i,j) = s_(i,j) + U(k,i) * dcmplx(dble(S(k,j)),dble(aimag(S(k,j))))
           enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call zgemm('C','N', shift2, shift2, sze,                       &
          (1.d0,0d0), U, size(U,1), W, size(W,1),                          &
          (0.d0,0d0), h, size(h,1))

      call zgemm('C','N', shift2, shift2, sze,                       &
          (1.d0,0d0), U, size(U,1), U, size(U,1),                          &
          (0.d0,0d0), s_tmp, size(s_tmp,1))

      ! Diagonalize h_p
      ! ---------------

       h_cp = h
       s_cp = s_tmp
       call lapack_zggev(h_cp,s_tmp,size(h,1),shift2,lambda,y_left,y,info)

       ! Normalization to have Y^* s_tmp Y = Id
       ! => y = 1/sqrt(y^* s_tmp y)
       ! --------------------------------------

       call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), s_cp, size(h,1), y, size(y,1),                          &
          (0d0,0.d0), s_tmp, size(s_tmp,1))

       call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), s_tmp, size(s_tmp,1),                  &
          (0.d0,0d0), s_cp, size(h,1))

       do i = 1, shift2
         y(:,i) = y(:,i) / cdsqrt(s_cp(i,i))
       enddo 

       if (info > 0) then
         ! Numerical errors propagate. We need to reduce the number of iterations
         itermax = iter-1
         exit
       endif

      ! Compute Energy for each eigenvector
      ! -----------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), h, size(h,1), y, size(y,1),                          &
          (0d0,0.d0), s_tmp, size(s_tmp,1))

      call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), s_tmp, size(s_tmp,1),                  &
          (0.d0,0d0), h, size(h,1))

      do k=1,shift2
        lambda(k) = h(k,k)
      enddo

      ! Compute S2 for each eigenvector
      ! -------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), s_, size(s_,1), y, size(y,1),                        &
          (0.d0,0d0), s_tmp, size(s_tmp,1))

      call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), s_tmp, size(s_tmp,1),                  &
          (0.d0,0d0), s_, size(s_,1))

      do k=1,shift2
        s2(k) = s_(k,k)
      enddo

      if (only_expected_s2) then
          do k=1,shift2
            state_ok(k) = (dabs(dble(s2(k))-expected_s2) < 0.6d0)
          enddo
      else
        do k=1,size(state_ok)
          state_ok(k) = .True.
        enddo
      endif

      do k=1,shift2
        if (.not. state_ok(k)) then
          do l=k+1,shift2
            if (state_ok(l)) then
              call zswap(shift2, y(1,k), 1, y(1,l), 1)
              call zswap(1, s2(k), 1, s2(l), 1)
              call zswap(1, lambda(k), 1, lambda(l), 1)
              state_ok(k) = .True.
              state_ok(l) = .False.
              exit
            endif
          enddo
        endif
      enddo

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------
      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), U, size(U,1), y, size(y,1), (0.d0,0d0), U(1,shift2+1), size(U,1))

      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), W, size(W,1), y, size(y,1), (0.d0,0d0), W(1,shift2+1), size(W,1))

      y_s(:,:) = cmplx(real(dble(y(:,:))), real(dimag(y(:,:))))

      call cgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.,0.), S, size(S,1), y_s, size(y_s,1), (0.,0.), S(1,shift2+1), size(S,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag
        do i=1,sze
           if (dble(H_jj(i) - lambda (k)) >= 1d-2) then 
             U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k) ) /(H_jj(i) - lambda (k))
           else
             U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
                / dcmplx(1d-2, dimag(H_jj(i) - lambda (k)))
           endif
        enddo

        if (k <= N_st) then
          call inner_product_complex(U(1,shift2+k),U(1,shift2+k),sze,residual_norm(k))
          to_print(1,k) = lambda(k) + dcmplx(nuclear_repulsion,0d0)
          to_print(2,k) = s2(k)
          to_print(3,k) = residual_norm(k)
        endif
      enddo
      !$OMP END PARALLEL DO

      if ((itertot>1).and.(iter == 1)) then
        !don't print
        continue
      else
        !write(*,'(1X,I3,1X,100(1X,F16.10,1X,F16.10,1X,F11.6,1X,F11.6,1X,ES11.3,1X,ES11.3))') iter-1, to_print(1:3,1:N_st)
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,F11.6,1X,ES11.3))') iter-1, dble(to_print(1:3,1:N_st))
      endif

      ! Check convergence
      if (iter > 1) then
        if (threshold_davidson_from_pt2) then
          converged = dabs(maxval(dble(residual_norm(1:N_st)))) < threshold_davidson_pt2
        else                                               
          converged = dabs(maxval(dble(residual_norm(1:N_st)))) < threshold_davidson
        endif
      endif

      do k=1,N_st
        if (dble(residual_norm(k)) > 1.d8) then
          print *, 'Davidson failed'
          stop -1
        endif
      enddo
      if (converged) then
        exit
      endif

      logical, external :: qp_stop
      if (qp_stop()) then
        converged = .True.
        exit
      endif

    enddo

    ! Re-contract U and update S and W
    ! --------------------------------

    call cgemm('N','N', sze, N_st_diag, shift2, (1.,0.),      &
        S, size(S,1), y_s, size(y_s,1), (0.,0.), S(1,shift2+1), size(S,1))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        S(i,k) = S(i,shift2+k)
      enddo
    enddo
    !$OMP END PARALLEL DO

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,1d0),      &
        W, size(W,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo
    !$OMP END PARALLEL DO

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        U, size(U,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    !$OMP END PARALLEL DO

  enddo

  ! Sign of each vector
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,val)
  do j = 1, N_st_diag 
    val = 0d0
    do i = 1, sze
      if (dabs(dble(U(i,j))) > dabs(val)) then
        val = dble(U(i,j))
        max_U(j) = dble(U(i,j))
        pos_U(j) = i
      endif
    enddo
  enddo
  !$OMP END PARALLEL DO

  call nullify_small_elements_complex(sze,N_st_diag,U,size(U,1),threshold_davidson_pt2)
  call ortho_qr_complex(U,size(U,1),sze,N_st_diag)
  call ortho_qr_complex(U,size(U,1),sze,N_st_diag)

  ! Change the sign of the guess vectors to match with the ones before the QR
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
  do j = 1, N_st_diag
    if (sign(1d0,max_U(j)) * sign(1d0,dble(U(pos_U(j),j))) < -1d-30) then
      do i = 1, sze
        u_in(i,j) = - U(i,j)
      enddo
    else
      do i = 1, sze
        u_in(i,j) = U(i,j)
      enddo
    endif
  enddo
  !$OMP END PARALLEL DO

  do k=1,N_st_diag
    energies(k) = lambda(k)
    s2_out(k) = s2(k)
  enddo
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================ =========== =========== =========== ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  deallocate(W,S)

  deallocate (                                                       &
      residual_norm,                                                 &
      U, overlap,                                                    &
      h, y_s, S_d,                                                   &
      y, s_, s_tmp,                                                  &
      lambda                                                         &
      )
  FREE nthreads_davidson
end







