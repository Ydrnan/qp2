subroutine write_wf()

    implicit none

    integer :: i
    character(len=256) :: fmt
    character*(2048)                :: output(2)

    write(fmt,*) '(I12,X,',3*N_states,'(X,ES24.12))'

    open(unit=11, file=trim(ezfio_filename)//'/work/wf_c.txt')   
    do i = 1, N_det
      write(11,fmt) i, psi_coef(i,:), psi_cap_coef(i,:)
    enddo
    close(11)
    open(unit=11, file=trim(ezfio_filename)//'/work/wf_d.txt')   
    do i = 1, N_det
      write(11,*) i
      call bitstring_to_str( output(1), psi_det(1,1,i), N_int )
      call bitstring_to_str( output(2), psi_det(1,2,i), N_int )
      write(11,*) trim(output(1))
      write(11,*) trim(output(2))
    enddo
    close(11)

end
