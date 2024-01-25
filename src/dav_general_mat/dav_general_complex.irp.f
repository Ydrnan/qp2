
subroutine davidson_general_complex(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag_in,converged,h_mat)
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! u_in : guess coefficients on the various states. Overwritten on exit
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
  integer, intent(in)             :: dim_in, sze, N_st, N_st_diag_in
  complex*16,  intent(in)   :: H_jj(sze),h_mat(sze,sze)
  complex*16, intent(inout) :: u_in(dim_in,N_st_diag_in)
  complex*16, intent(out)   :: energies(N_st)

  integer                        :: iter, N_st_diag
  integer                        :: i,j,k,l,m
  logical, intent(inout)         :: converged

  double precision, external     :: u_dot_v, u_dot_u

  integer                        :: k_pairs, kl

  integer                        :: iter2, itertot
  complex*16, allocatable        :: y(:,:), h(:,:), lambda(:)
  complex*16                     :: diag_h_mat_elem
  complex*16, allocatable        :: residual_norm(:)
  character*(16384)              :: write_buffer
  complex*16                     :: to_print(2,N_st)
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax, istate
  double precision               :: r1, r2, alpha
  integer                        :: nproc_target, info
  integer                        :: order(N_st_diag_in)
  double precision               :: cmax
  complex*16      , allocatable  :: U(:,:), overlap(:,:)!, S_d(:,:)
  complex*16      , pointer      :: W(:,:)
  logical                        :: disk_based
  complex*16                     :: energy_shift(N_st_diag_in*davidson_sze_max)

  include 'constants.include.F'

  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, y, h, lambda
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
  overlap = 0.d0

  provide threshold_davidson !nthreads_davidson
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
  maxab = sze 

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 2d0 * 8.d0 *                                   &!complex  bytes
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.d0*dble(sze*m)*(N_st_diag*itermax)  &! W
         + 2.0d0*(N_st_diag*itermax)**2           &! h,y
         + 2.d0*(N_st_diag*itermax)               &! s2,lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
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
  call write_int(6,sze,'Number of basis functions')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================  ===========  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = 'Iter'
  do i=1,N_st
    if (i==1) then
    write_buffer = trim(write_buffer)//'       Energy            Im            Residual       Im     '
    else 
    write_buffer = trim(write_buffer)//'            Energy            Im            Residual       Im     '
    endif
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================  ===========  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)


!  if (disk_based) then
!    ! Create memory-mapped files for W and S
!    type(c_ptr) :: ptr_w, ptr_s
!    integer :: fd_s, fd_w
!    call mmap(trim(ezfio_work_dir)//'davidson_w', (/int(sze,8),int(N_st_diag*itermax,8)/),&
!        8, fd_w, .False., ptr_w)
!    call mmap(trim(ezfio_work_dir)//'davidson_s', (/int(sze,8),int(N_st_diag*itermax,8)/),&
!        4, fd_s, .False., ptr_s)
!    call c_f_pointer(ptr_w, w, (/sze,N_st_diag*itermax/))
!    call c_f_pointer(ptr_s, s, (/sze,N_st_diag*itermax/))
!  else
    allocate(W(sze,N_st_diag*itermax))
!  endif

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax),                                      &
      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      residual_norm(N_st_diag),                                      &
      lambda(N_st_diag*itermax))

  !print*,'H'
  !do i = 1, sze
  !   do j = 1, sze
  !      print*, i-1,j-1,REALPART(h_mat(i,j))
  !      !write(*,'(100(F10.4))') REALPART(h_mat(i,:))
  !   enddo
  !enddo
  !call diag_general_complex(lambda,y,h_mat,size(h_mat,1),sze)
  !print*,''
  !print*,lambda(1:sze)
  !call abort()

  h = (0.d0, 0d0)
  U = (0.d0, 0d0)
  y = (0.d0, 0d0)


  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)

  ! Davidson iterations
  ! ===================

  converged = .False.

  ! Initialize from N_st to N_st_diat with gaussian random numbers
  ! to be sure to have overlap with any eigenvectors
  do k=N_st+1,N_st_diag
    u_in(k,k) = 10.d0
    do i=1,sze
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      u_in(i,k) = dcmplx(r1*dcos(r2), r1)
    enddo
  enddo
  ! Normalize all states 
  complex*16 :: norm
  do k=1,N_st_diag
    !call normalize(u_in(:,k),sze)
    call normalize_complex(u_in(:,k),sze)
    !call inner_product_complex(u_in(:,k),u_in(:,k),sze,norm)
    !print*,k,norm
  enddo

  ! Copy from the guess input "u_in" to the working vectors "U"
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

    do iter=1,itermax-1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter
   
      !print*,iter,itertot
      !print*,'cond',(iter > 1).or.(itertot == 1), (iter > 1), (itertot == 1)
      !if ((iter > 1).or.(itertot == 1)) then
        ! Compute |W_k> = \sum_i |i><i|H|u_k>
        ! -----------------------------------

         ! Gram-Smitt to orthogonalize all new guess with the previous vectors 
          call ortho_qr_complex(U,size(U,1),sze,shift2)
          call ortho_qr_complex(U,size(U,1),sze,shift2)

          !print*,'check ortho'
          !do k = 1, N_st_diag
          !  do i = 1, N_st_diag
          !  call inner_product_complex(U(:,k),U(:,i),sze,norm)
          !  print*,k,i,norm
          !  enddo
          !enddo
          !print*,'W',U(:,shift+1)
          !print*,'U',U(:,shift+1)

          !print*,'hpsi'
         call hpsi_complex(W(:,shift+1),U(:,shift+1),N_st_diag,sze,h_mat)
          !print*,'done'
      !else
         ! Already computed in update below
      !   print*,'Already computed in update below'
      !   continue
      !endif

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      !print*, 'U* W', shift2
      call zgemm('C','N', shift2, shift2, sze,                       &
          (1.d0,0d0), U, size(U,1), W, size(W,1),                          &
          (0.d0,0d0), h, size(h,1))

      ! Diagonalize h y = lambda y
      ! ---------------

      !print*,'diag'
      !print*,'h',h
      !call lapack_diag(lambda,y,h,size(h,1),shift2)
      call diag_general_complex(lambda,y,h,size(h,1),shift2,info)
      !print*,'diag'
      !write(*,'(100(F11.6))') lambda(1:shift2)

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------

      !print*,'Uy'
      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), U, size(U,1), y, size(y,1), (0d0,0.d0), U(:,shift2+1), size(U,1))
      !print*,'Wy'
      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), W, size(W,1), y, size(y,1), (0.d0,0d0), W(:,shift2+1), size(W,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------
      !print*,'residual'
      do k=1,N_st_diag
        do i=1,sze
          U(i,shift2+k) =  &
            (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
              /(H_jj(i) - lambda (k))
            ! Removed max here
        enddo
        !print*,'new U',U(:,shift2+k)

        if (k <= N_st) then
          !residual_norm(k) = u_dot_u(U(:,shift2+k),sze)
          call inner_product_complex(U(:,shift2+k),U(:,shift2+k),sze, residual_norm(k))
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = residual_norm(k)
        endif
      enddo

      if ((itertot>1).and.(iter == 1)) then
        !don't print 
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,F16.10,1X,ES12.3,1X,ES12.3,1X))') iter-1, to_print(1:2,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
          converged = dabs(maxval(dble(residual_norm(1:N_st)))) < threshold_davidson
      endif   
      

      do k=1,N_st
        if (dble(residual_norm(k)) > 1.e8) then
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

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        W, size(W,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        U, size(U,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    !call ortho_qr_complex(U,size(U,1),sze,N_st_diag)
    !call ortho_qr_complex(U,size(U,1),sze,N_st_diag)
    !do j=1,N_st_diag
    !  k=1
    !  do while ((k<sze).and.(U(k,j) == (0.d0,0d0)))
    !    k = k+1
    !  enddo
    !  if (dble(U(k,j) * u_in(k,j)) < 0.d0) then
    !    do i=1,sze
    !      W(i,j) = -W(i,j)
    !    enddo
    !  endif
    !enddo
  enddo

  do k=1,N_st
    energies(k) = lambda(k)
  enddo
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================  ===========  ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''

  call write_time(6)
  print*,'Energies:'
  call ortho_qr_complex(u_in,size(u_in,1),sze,N_st_diag)
  call check_energy(h_mat,u_in,N_st,sze)
  print*,'c-Energies:'
  call modified_gram_schmidt_c(u_in,N_st_diag,sze)
  call check_c_energy(h_mat,u_in,N_st,sze)
  write(6,'(A)') ''

    deallocate(W)

  deallocate (                                                       &
      residual_norm,                                                 &
      U, h,                                                 &
      y,                                                             &
      lambda                                                         &
      )
  deallocate(overlap)
  FREE nthreads_davidson
end

subroutine hpsi_complex(v,u,N_st,sze,h_mat)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v = H | u \rangle$ and 
  END_DOC
  integer, intent(in)              :: N_st,sze
  complex*16, intent(in)     :: u(sze,N_st),h_mat(sze,sze)
  complex*16, intent(inout)  :: v(sze,N_st)
  integer :: i,j,istate
  v = (0.d0, 0d0)
  !call zgemv('N', sze, sze, (1d0,0d0), h_mat, sze, u, 1, (1d0,0d0), v, 1)
  do istate = 1, N_st
   do i = 1, sze
    do j = 1, sze
      v(i,istate) += h_mat(j,i) * u(j,istate)
    enddo
   enddo
  enddo
end

subroutine check_energy(h,psi,N_st,sze)
  implicit none

  integer, intent(in) :: N_st, sze
  complex*16, intent(in) :: h(sze,sze), psi(sze,N_st)

  complex*16 :: res
  complex*16, allocatable :: h_psi(:,:), energy(:)
  integer :: i,j

  allocate(h_psi(sze,N_st),energy(N_st))

  call hpsi_complex(h_psi,psi,N_st,sze,h)

  do j = 1, N_st
    energy(j) = (0d0,0d0)
    do i = 1, sze
      energy(j) += CONJG(psi(i,j)) * h_psi(i,j)
    enddo
    call inner_product_complex(psi(1,j),psi(1,j),sze,res)
    energy(j) = energy(j) / res + dcmplx(nuclear_repulsion,0d0)
    write(*,'(I6,2(F18.10))') j, dble(energy(j)), dimag(energy(j))
  enddo
  
end

subroutine check_c_energy(h,psi,N_st,sze)
  implicit none

  integer, intent(in) :: N_st, sze
  complex*16, intent(in) :: h(sze,sze), psi(sze,N_st)

  complex*16 ::norm
  complex*16, allocatable :: h_psi(:,:), energy(:)
  integer :: i,j

  allocate(h_psi(sze,N_st),energy(N_st))

  call hpsi_complex(h_psi,psi,N_st,sze,h)

  do j = 1, N_st
    energy(j) = (0d0,0d0)
    do i = 1, sze
      energy(j) += psi(i,j) * h_psi(i,j)
    enddo
    call inner_prod_c(psi(1,j),psi(1,j),sze,norm)
    !call c_norm(psi(1,j),sze,norm)
    energy(j) = energy(j) / norm + dcmplx(nuclear_repulsion,0d0)
    write(*,'(I6,2(F18.10))') j, dble(energy(j)), dimag(energy(j))
  enddo

end

subroutine modified_gram_schmidt_c(v,N_st,sze)
  implicit none

  integer, intent(in) :: N_st, sze
  complex*16, intent(inout) :: v(sze,N_st)

  complex*16, allocatable :: u(:,:), p(:)
  integer :: i,j
  complex*16 :: res

  allocate(u(sze,N_st),p(sze))

  u = (0d0,0d0)

  do i = 1, N_st
    u(:,i) = v(:,i)
    do j = 1, i-1
      if (j==1) then
        call proj_c(u(1,j),v(1,i),sze,p)
      else
        call proj_c(u(1,j),u(1,i),sze,p)
      endif
      u(:,i) = u(:,i) - p(:)
    enddo
  enddo
  do i = 1, N_st
    call normalize_c(u(1,i),sze)
  enddo
  
  v = u

  ! Check orthogonality
  do i = 1, N_st
    call inner_prod_c(u(1,i),u(1,i),sze,res)
    if (dabs(cdabs(res) - 1d0) > 1d-12) then
      print*,'Gram-Schmidt orthogonalization failed.'
      call abort()
    endif
  enddo

  do i = 2, N_st
    do j = 1, i-1
    call inner_prod_c(u(1,i),u(1,j),sze,res)
    if (cdabs(res) > 1d-12) then
      print*,'Gram-Schmidt orthogonalization failed.'
      call abort()
    endif
    enddo
  enddo

end

subroutine proj_c(u,v,sze,p)
  implicit none

  integer, intent(in) :: sze
  complex*16, intent(in) :: u(sze), v(sze)
  complex*16, intent(out) :: p(sze)
  complex*16 :: vu,uu, f
  integer :: i

  call inner_prod_c(v,u,sze,vu)
  call inner_prod_c(u,u,sze,uu)

  f = vu/uu

  do i = 1, sze
    p(i) = f * u(i)
  enddo
  
end

subroutine inner_prod_c(u,v,sze,res)
  implicit none
  
  integer, intent(in) :: sze
  complex*16, intent(in) :: u(sze), v(sze)
  complex*16, intent(out) :: res

  integer :: i
  
  res = (0d0,0d0)
  do i = 1, sze
    res = res + u(i) * v(i)
  enddo
  
end

subroutine normalize_c(u,sze)
  implicit none

  integer, intent(in) :: sze
  complex*16, intent(inout) :: u(sze)

  complex*16 :: n, z
  !double precision :: phi
  integer :: i

  call inner_prod_c(u,u,sze,z)

  !phi = datan2(IMAGPART(z), REALPART(Z))
  !n = (1d0,0d0) / cdsqrt(dcmplx(cdabs(z),0d0)) * cdexp(dcmplx(0d0,-0.5d0*phi))
  n = (1d0,0d0) / cdsqrt(z)
  do i = 1, sze
    u(i) = u(i) * n
  enddo

end

subroutine inner_product_complex(u,v,sze,res)
  implicit none
  BEGIN_DOC
  ! Computes the c-norm (u|u) of a complex wave function u.
  END_DOC

  integer, intent(in) :: sze
  complex*16, intent(in) :: u(sze), v(sze)
  complex*16, intent(out) :: res
  integer :: i

  res = (0d0,0d0)

  do i = 1, sze
    !print*,u(i),v(i),DCONJG(v(i))
    res = res + u(i) * dconjg(v(i))
  enddo

end

subroutine normalize_complex(u,sze)
  implicit none
  BEGIN_DOC
  ! Normalizes u s.t. <u|u> = 1.0.
  END_DOC

  integer, intent(in) :: sze
  complex*16, intent(inout) :: u(sze)
  complex*16 :: norm
  integer :: i

  call inner_product_complex(u,u,sze,norm)

  do i = 1, sze
    u(i) = u(i) * (1d0,0d0)/cdsqrt(norm)
  enddo

end

subroutine sort_complex(v,iorder,sze)
  implicit none
  integer, intent(in) :: sze
  complex*16, intent(inout) :: v(sze)
  integer, intent(out) :: iorder(sze)
  complex*16, allocatable :: tmp(:)
  double precision, allocatable :: re(:), im(:)

  integer :: i,j,k

  allocate(tmp(sze),re(sze),im(sze))

  do i = 1, sze
   iorder(i) = i
  enddo

  do i = 1, sze
    tmp(i) = v(i)
  enddo

  do i = 1, sze
    re(i) = dble(v(i))
  enddo

  call dsort(re, iorder, sze)

  j = 0
  do i = 2, sze
    if (re(i) == re(i-1)) then
      j = j + 1
      im(j) = dimag(v(iorder(i-1)))
    else if (j /= 0) then
      j = j + 1
      im(j) = dimag(v(iorder(i-1)))
      call dsort(im, iorder(i-j:i-1),j)
      j = 0
    endif
  enddo

  do i = 1, sze
    v(i) = tmp(iorder(i))
  enddo
end

subroutine diag_general_complex(eigvalues,eigvectors,H,nmax,n,info)
  implicit none
  BEGIN_DOC
  ! Diagonalization of a complex matrix
  END_DOC

  integer, intent(in) :: nmax,n
  complex*16, intent(in) :: H(nmax,n)
  complex*16, intent(out) :: eigvalues(n), eigvectors(nmax,n)
  integer, intent(out) ::info
  complex*16, allocatable :: w(:), vl(:,:), vr(:,:), work(:)
  double precision, allocatable :: rwork(:)
  integer :: lwork,i,j
  integer, allocatable :: iorder(:)

  lwork = max(1,2*n)

  allocate(w(n),vl(nmax,n),vr(nmax,n),work(lwork), rwork(2*n),iorder(n))
  !print*,nmax,n
  !print*,size(H,1),size(H,2)
  !print*,'H',H
  call zgeev('N', 'V', n, H, nmax, w, vl, nmax, vr, nmax, work, lwork, rwork, info)

  if (info < 0) then
    print*,' the ', abs(info), '-th argument had an illegal value.'
  else if (info > 0) then
    print*, 'the QR algorithm failed to compute all the eigenvalues.'
  endif

  call sort_complex(w, iorder, n)
  eigvalues = w

  do j = 1, n
    do i = 1, n
      eigvectors(i,j) = vr(i,iorder(j))
    enddo
  enddo

end

subroutine lapack_zggev(A,B,nmax,n,eigvalues,l_eigvectors,r_eigvectors,info)
  implicit none

  integer, intent(in)    :: nmax, n
  complex*16, intent(in) :: A(nmax,n), B(nmax,n)
  complex*16, intent(out) :: eigvalues(n), l_eigvectors(nmax,n), r_eigvectors(nmax,n)
  integer, intent(out) ::info
  complex*16, allocatable :: alpha(:), beta(:), vl(:,:), vr(:,:), work(:), e(:), tmp(:)
  double precision, allocatable :: rwork(:)
  integer :: lwork,i,j
  integer, allocatable :: iorder(:)
    
  lwork = max(1,2*n)

  allocate(alpha(n), beta(n), e(n), vl(nmax,n), vr(nmax,n), work(lwork), rwork(max(1,8*n)), iorder(n), tmp(n))

  call zggev('V','V', n, a, size(a,1), b, size(b,1), alpha, beta, vl, size(vl,1), vr, size(vr,1), work, lwork, rwork, info)

  if (info < 0) then
    print*,' the ', abs(info), '-th argument had an illegal value.'
  else if (info > 0) then
    print*, 'the QZ iteration failed. No eigenvectors have been calculated.'
  endif

  do i = 1, n
    if (cdabs(beta(i)) < 1d-10) then
      info = 1
      e(i) = dcmplx(0d0,0d0)
    else
      e(i) = alpha(i)/beta(i)
    endif
  enddo

  tmp = e
  call sort_complex(e, iorder, n)

  do i = 1, n
    eigvalues(i) = tmp(iorder(i))
  enddo

  do j = 1, n
    do i = 1, n
      r_eigvectors(i,j) = vr(i,iorder(j))
    enddo
  enddo

  do j = 1, n
    do i = 1, n
      l_eigvectors(i,j) = vl(i,iorder(j))
    enddo
  enddo

end

BEGIN_PROVIDER [ complex*16, H_matrix_all_dets_complex,(N_det,N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |H| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i,j,k,h1,p1,h2,p2,s1,s2
 double precision :: hij, wij, phase
 double precision, external :: diag_H_mat_elem_cap
 integer :: degree, exc(0:2,2,2)
 call  i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 print*,'Providing the H_matrix_all_dets_complex ...'
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,wij,hij,degree,k,&
 !$OMP phase, h1,p1,h2,p2,s1,s2,exc) &
 !$OMP SHARED (N_det, psi_det, N_int,H_matrix_all_dets_complex)
 do i =1,N_det
   do j = i, N_det
    call  i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)

    if (degree == 0) then
      wij = diag_H_mat_elem_cap(psi_det(1,1,i),N_int)
    else if (degree == 1) then
      call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
      call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
      call i_W_j_single_spin_cap(psi_det(1,1,i),psi_det(1,1,j),N_int,s1,wij)   
    else
      wij = 0d0
    endif

    H_matrix_all_dets_complex(i,j) = dcmplx(hij,wij)
    H_matrix_all_dets_complex(j,i) = dcmplx(hij,wij)
  enddo
 enddo
 !$OMP END PARALLEL DO
 print*,'H_matrix_all_dets_complex done '
END_PROVIDER

subroutine test_ortho_qr_complex(A, lda, m, n)
  implicit none

  integer, intent(in) :: lda, m, n
  complex*16, intent(inout) :: A(lda,n)
  
  integer :: lwork, info
  complex*16, allocatable :: work(:), tau(:)

  lwork = -1
  allocate(work(1),tau(min(m, n)))

  if (m < n) then
    print*,'Error in parameters for qr'
    call abort()
  endif
  
  call zgeqrf(m, n, A, lda, tau, work, lwork, info)

  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))

  call zgeqrf(m, n, A, lda, tau, work, lwork, info)

  if (info < 0) then
     print*,'the', info,'-th parameter had an illegal value.'
  endif
  
  lwork = -1
  deallocate(work)
  allocate(work(1))
  call zungqr(m, n, n, A, lda, tau, work, lwork, info) 

  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))

  call zungqr(m, n, n, A, lda, tau, work, lwork, info)
  
  if (info < 0) then
     print*,'the', info,'-th parameter had an illegal value.'
  endif

  deallocate(tau,work)
end

subroutine test_ortho_qr_complex2(A, lda, m, n)
  implicit none

  integer, intent(in) :: lda, m, n
  complex*16, intent(inout) :: A(lda,n)

  integer :: lwork, info
  complex*16, allocatable :: work(:), tau(:)
  integer, allocatable :: jpvt(:)
  double precision, allocatable :: rwork(:)

  lwork = -1
  allocate(work(1),tau(min(m, n)),jpvt(n),rwork(2*n))

  call zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
  
  jpvt = 1
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))

  call zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)

  if (info < 0) then
     print*,'the', info,'-th parameter had an illegal value.'
  endif

  lwork = -1
  deallocate(work)
  allocate(work(1))
  call zungqr(m, n, n, A, lda, tau, work, lwork, info)

  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))

  call zungqr(m, n, n, A, lda, tau, work, lwork, info)

  if (info < 0) then
     print*,'the', info,'-th parameter had an illegal value.'
  endif

  deallocate(tau,work)

end
