program testFRB
    use FRB
    implicit none
    integer, parameter :: nsamples=10000
    double precision, dimension(3,nsamples) :: FRB_pos
    integer :: u, i

    do i=1,nsamples
        FRB_pos(:,i) = get_FRB_position(0.1d0, 3d0, 1d10, 5d10, 10d0)
    enddo

    open (newunit=u, file='testFRB_position.dat', status='new',form='unformatted')
    write (u) FRB_pos
    close (u)


end program testFRB
