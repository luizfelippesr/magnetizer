program testFRB
    use FRB
    implicit none
    integer, parameter :: nsamples=10000
    double precision, dimension(3,nsamples) :: FRB_pos
    integer :: u, i

    ! The following will generate nsamples positions (x,y,z) of FRBs to check
    ! whether this is working, one can open testFRB_position.dat and plot:
    ! * the distribution of z showing that this is proportional to exp(-|z|/0.1)
    ! * the distribution of sqrt(x^2+y^2) and check whether it qualitatively
    !   agrees with the expected molecular gas density distribution
    ! For the (unphysical) choice below with large gas mass and no stars, the
    ! approximately all the gas will in molecular phase, therefore the surface
    ! density should be that of an exponential disc
    do i=1,nsamples
        FRB_pos(:,i) = get_FRB_position(h_FRB=0.1d0, r_disk=3d0, Mgas_disk=1d13, &
                                        Mstars_disk=0d0, rmax=10d0)
    enddo
    open (newunit=u, file='testFRB_position.dat', status='new',form='unformatted')
    write (u) FRB_pos
    close (u)

end program testFRB
