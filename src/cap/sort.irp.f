! Recursive sort to sort the determinants based on the integers that compose them
! The sort is done on the first integer, then the second, ...
recursive subroutine recursive_int_sort(l_det,key,sze,Nint,idx)
  use bitmasks
  implicit none

  integer, intent(in) :: sze, Nint, idx
  integer(bit_kind), intent(inout) :: l_det(2*Nint, sze)
  integer, intent(inout) :: key(sze)

  integer :: i,j,k,l,nb_u
  integer, allocatable :: iorder(:),nu(:),pu(:)

  if (sze == 0) return

  if (idx < 2*Nint) then

    ! Sort
    call multiple_int_sort(l_det,key,sze,Nint,idx)

    allocate(pu(sze),nu(sze))
    ! Unique, nb and position
    call search_unique_int(l_det,sze,Nint,idx,nb_u,nu,pu)

    do i = 1, nb_u
      call recursive_int_sort(l_det(1,pu(i)),key(pu(i)),nu(i),Nint,idx+1)
    enddo
    deallocate(pu,nu)

  else

    ! Sort
    call multiple_int_sort(l_det,key,sze,Nint,idx)

  endif

end

subroutine multiple_int_sort(l_det,key,sze,Nint,idx)
  use bitmasks
  implicit none

  integer, intent(in) :: sze,Nint,idx
  integer(bit_kind), intent(inout) :: l_det(2*Nint,sze)
  integer, intent(inout) :: key(sze)

  integer :: i,j,k,l,val
  integer(bit_kind), allocatable :: tmp(:),tmp_int(:,:)
  integer, allocatable :: iorder(:), tmp_k(:)

  ! Sort
  allocate(tmp(sze),tmp_int(2*Nint,sze),tmp_k(sze),iorder(sze))

  do i = 1, sze
    tmp(i) = l_det(idx,i)
    tmp_int(:,i) = l_det(:,i)
    tmp_k(i) = key(i)
    iorder(i) = i
  enddo

  call i8sort(tmp,iorder,sze)

  do i = 1, sze
    l_det(:,i) = tmp_int(:,iorder(i))
    key(i) = tmp_k(iorder(i))
  enddo

  deallocate(tmp,tmp_int,tmp_k,iorder)
end

subroutine search_unique_int(l_det,sze,Nint,idx,nb_u,nu,pu)
  use bitmasks
  implicit none

  integer, intent(in) :: sze,Nint,idx
  integer(bit_kind), intent(in) :: l_det(2*Nint,sze)
  integer(bit_kind) :: val

  integer, intent(out) :: nb_u, nu(sze), pu(sze)

  integer :: i,j,k,l

  ! Unique, nb and position
  k = 1
  pu = 0 ! starting position
  nu = 0 ! nb
  pu(1) = 1
  nu(1) = 1
  val = l_det(idx,1)
  do i = 2, sze
    if (val /= l_det(idx,i)) then
      k = k + 1
      pu(k) = i
      nu(k) = nu(k) + 1
      val = l_det(idx,i)
    else
      nu(k) = nu(k) + 1
    endif
  enddo

  nb_u = k

end

