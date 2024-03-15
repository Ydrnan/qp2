subroutine qr_decomposition_c(A,lda,m,n)
  implicit none

  integer, intent(in) ::  lda,m,n 
  complex*16, intent(inout) :: A(lda,n)

  complex*16, allocatable :: R(:,:)
  integer :: i,j,k

  allocate(R(n,n))

  !do i = 1, lda
  !  write(*,'(100(ES12.4))') A(i,:)
  !enddo

  do j = 1, n
    do i = 1, n
       R(i,j) = (0d0,0d0)
    enddo
  enddo

  do j = 1, n
    do i = 1, j-1
      do k = 1, m
        R(i,j) = R(i,j) + A(k,i) * A(k,j)
        !R(i,j) = R(i,j) + dconjg(A(k,i)) * A(k,j)
      enddo
      do k = 1, m
        !A(:,j) = A(:,j) - R(i,j) * A(:,i)
        A(k,j) = A(k,j) - R(i,j) * A(k,i)
      enddo
    enddo

    do k = 1, m
      R(j,j) = R(j,j) + A(k,j) * A(k,j)
      !R(j,j) = R(j,j) + dconjg(A(k,j)) * A(k,j)
    enddo
    R(j,j) = cdsqrt(R(j,j))
    do k = 1, m
      !A(:,j) = A(:,j) / R(j,j)
      A(k,j) = A(k,j) / R(j,j)
    enddo
  enddo

  deallocate(R)

end

subroutine normalize_xt_a_x(A,lda,X,ldx,n)

  implicit none

  integer, intent(in) :: lda, ldx, n
  complex*16, intent(in) :: A(lda,n)
  complex*16, intent(inout) :: X(ldx,n)

  complex*16, allocatable :: A_cp(:,:), e(:), tmp(:,:), X_tmp(:,:), eigvec(:,:), diag(:,:)
  complex*16 :: res
  integer :: info, i,j

  allocate(A_cp(n,n),e(n),tmp(n,n),X_tmp(n,n),eigvec(n,n),diag(n,n))
    
  A_cp(1:n,1:n) = A(1:n,1:n)
  print*,'A'
  do i = 1,n
    write(*,'(100(ES12.3))') A_cp(i,:)
  enddo
  print*,''
  do i = 1, n
    call normalize_c(A_cp(:,i),n)
  enddo
  print*,'A'
  do i = 1,n
    write(*,'(100(ES12.3))') A_cp(i,:)
  enddo
  print*,''
  call diag_general_complex(e,eigvec,A_cp,n,n,info)
  call qr_decomposition_c(eigvec,size(eigvec,1),n,n)
  !print*,'e',e
  !do i = 1,n
  !  write(*,'(100(ES12.3))') eigvec(i,:)
  !enddo
  !print*,''

  diag = (0d0,0d0)
  do i = 1, n
    if (cdabs(e(i)) > 1d-12) then
      diag(i,i) = e(i)
    endif
  enddo

  call zgemm('N','N',n,n,n, (1d0, 0d0), eigvec, n, diag, n, (0d0,0d0), tmp, n)
  call zgemm('N','T',n,n,n, (1d0, 0d0), tmp, n, eigvec, n, (0d0,0d0), X_tmp, n)
  print*,'bt'
  do i = 1, n
    write(*,'(100(ES12.3))') X_tmp(i,:) 
  enddo
  print*,''

 do i = 1, n
    if (cdabs(e(i)) > 1d-12) then
      diag(i,i) = (1d0,0d0) / cdsqrt(e(i))
    endif
  enddo

  call zgemm('N','N',n,n,n, (1d0, 0d0), eigvec, n, diag, n, (0d0,0d0), tmp, n)
  call zgemm('N','T',n,n,n, (1d0, 0d0), tmp, n, eigvec, n, (0d0,0d0), X_tmp, n)
  X(1:n,1:n) = X_tmp(1:n,1:n)

  call zgemm('T','N',n,n,n, (1d0, 0d0), X, ldx, A, lda, (0d0,0d0), tmp, n)
  call zgemm('N','N',n,n,n, (1d0, 0d0), tmp, n, X, ldx, (0d0,0d0), X_tmp, n)
!
!  print*,'error'
  do i = 1, n
    X_tmp(i,i) = X_tmp(i,i) - (1d0,0d0)
!    write(*,'(100(ES12.3))') X_tmp(i,:) 
  enddo

  print*,'check', maxval(cdabs(X_tmp))

 double precision, allocatable :: dA(:,:), dX(:,:)
 allocate(dA(n,n), dX(n,n))
 do j = 1, n
   do i = 1, n
     dA(i,j) = dble(A(i,j))
   enddo
   da(j,j) = da(j,j) + 1d0
 enddo
 do i = 1, n
   write(*,'(100(ES12.3))') da(i,:)
 enddo
 dX = 0d0
 call dnormalize_xt_a_x(dA,dX,n)

 deallocate(tmp,x_tmp,diag,e,eigvec)

end

subroutine dnormalize_xt_a_x(A, X, n)
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: A(n, n)
  double precision, intent(inout) :: X(n, n)

  double precision, allocatable :: A_cp(:,:), e(:), tmp(:,:), X_tmp(:,:), eigvec(:,:), diag(:,:)
  double precision :: res
  integer :: i

  allocate(A_cp(n, n), e(n), tmp(n, n), X_tmp(n, n), eigvec(n, n), diag(n, n))

  A_cp(1:n, 1:n) = A(1:n, 1:n)

  print*,1
  do i = 1, n
    write(*, '(100(ES12.3))') A_cp(i, :)
  enddo
  print*, ''

  call lapack_diag(e, eigvec, A_cp, n, n)

  diag = 0d0
  do i = 1, n
      if (dabs(e(i)) > 1d-12) then
        diag(i, i) = e(i)
      endif
  enddo

  call dgemm('N', 'N', n, n, n, 1d0, eigvec, n, diag, n, 0d0, tmp, n)
  call dgemm('N', 'T', n, n, n, 1d0, tmp, n, eigvec, n, 0d0, X_tmp, n)

  print*,2
  do i = 1, n
    write(*, '(100(ES12.3))') X_tmp(i, :)
  enddo
  print*, ''

  write(*, '(A,100(ES12.3))') 'e',e(:)
  do i = 1, n
    if (e(i) > 1d-12) then
      diag(i,i) = 1d0 / dsqrt(e(i))
    endif
  enddo

  call dgemm('N','N',n,n,n, 1d0, eigvec, n, diag, n, 0d0, tmp, n)
  call dgemm('N','T',n,n,n, 1d0, tmp, n, eigvec, n, 0d0, X_tmp, n)
  X(1:n,1:n) = X_tmp(1:n,1:n)

  call dgemm('T','N',n,n,n, 1d0, X, n, A, n, 0d0, tmp, n)
  call dgemm('N','N',n,n,n, 1d0, tmp, n, X, n, 0d0, X_tmp, n)

  print*,'error'
  do i = 1, n
    X_tmp(i,i) = X_tmp(i,i) - 1d0
    write(*,'(100(ES12.3))') X_tmp(i,:)
  enddo

  print*,'check', maxval(dabs(X_tmp))


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
  complex*16 :: res
  integer :: i

  call inner_product_complex(u,u,sze,res)

  do i = 1, sze
    u(i) = u(i) * (1d0,0d0)/cdsqrt(res)
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
  complex*16, allocatable :: w(:), vl(:,:), vr(:,:), work(:), tmp(:)
  double precision, allocatable :: rwork(:)
  integer :: lwork,i,j
  integer, allocatable :: iorder(:)

  lwork = -1

  allocate(w(n),vl(nmax,n),vr(nmax,n),work(1), rwork(2*n),iorder(n), tmp(n))

  call zgeev('N', 'V', n, H, nmax, w, vl, nmax, vr, nmax, work, lwork, rwork, info)

  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))

  call zgeev('N', 'V', n, H, nmax, w, vl, nmax, vr, nmax, work, lwork, rwork, info)

  if (info < 0) then
    print*,' the ', abs(info), '-th argument had an illegal value.'
    !print*,H
  else if (info > 0) then
    print*, 'the QR algorithm failed to compute all the eigenvalues.'
  endif

  tmp = w
  call sort_complex(w, iorder, n)

  do i = 1, n
    eigvalues(i) = tmp(iorder(i))
  enddo

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

  !print*,'a'
  !do i = 1, n
  !    write(*,'(100(ES12.3))') dble(a(i,:))
  !enddo
  !print*,'b'
  !do i = 1, n
  !    write(*,'(100(ES12.3))') dble(b(i,:))
  !enddo
  !print*,''
    
  lwork = -1
  allocate(alpha(n), beta(n), e(n), vl(nmax,n), vr(nmax,n), work(1), rwork(max(1,8*n)), iorder(n), tmp(n))

  call zggev('V','V', n, a, size(a,1), b, size(b,1), alpha, beta, vl, size(vl,1), vr, size(vr,1), work, lwork, rwork, info)

  lwork=int(work(1))
  deallocate(work)
  allocate(work(lwork))

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

  !print*,'e',dble(e)
  tmp = e
  call sort_complex(e, iorder, n)
  !print*,'e',dble(e)

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

complex*16 function trace_complex(a,n)
  implicit none

  integer, intent(in) :: n
  complex*16, intent(in) :: a(n,n)
    
  integer :: i

  trace_complex = (0d0,0d0)
  do i = 1, n
    trace_complex += a(i,i)
  enddo

end
