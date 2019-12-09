program testFRB
    use FRB
    use tools, only: linspace
    implicit none
    double precision, parameter :: pi=     3.14156295358d0
    integer, parameter :: nsamples=10000
    integer, parameter :: npoints=1001
    double precision, dimension(3,nsamples) :: FRB_pos
    double precision, dimension(npoints) :: x, y, z
    double precision :: theta
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

    ! The follwing will FRB positions along a particular LoS. The distribution
    ! of positions should again follow the weighting by the molecular density.
    z = linspace(-2.5d0,2.5d0,npoints)
    y = 2.0
    theta = 30.d0  * pi/180d0
    x = z / tan(theta)

    do i=1,nsamples
        FRB_pos(:,i) = get_FRB_LoS_position(x,y(1),z,h_FRB=0.033d0, r_disk=3d0, Mgas_disk=1d13, &
                                        Mstars_disk=0d0, rmax=10d0)
    enddo
    open (newunit=u, file='testFRB_LoS_position.dat', status='new',form='unformatted')
    write (u) FRB_pos
    close (u)

end program testFRB
