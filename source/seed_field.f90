module seed_field  !Set initial conditions
  use global_input_parameters
  use grid
  use random
  use profiles, only: l, h, Beq, delta_r_floor
  use floor_field
  implicit none

  contains
    subroutine init_seed_field(f)
      integer :: i,var
      double precision, dimension(nx) :: Bseed
      double precision, dimension(:,:), intent(inout) :: f
      double precision :: r1
      ! Initializes the seed field to a fraction of the equipartition field
      Bseed=frac_seed*Beq
      ! Initializes f
      f(:,:)=0.d0

      select case (trim(p_seed_choice))
        case('fraction')
          ! The field is a fixed fraction of the equipartition field
          f(:,1)=-Bseed
          f(:,2)= Bseed
        case('decaying')
          ! The field
          r1 = p_r_seed_decay/r_max_kpc
          f(:,1)=-Bseed*r*(1.d0-r)**p_nn_seed*dexp(-r/r1)
          f(:,2)= Bseed*r*(1.d0-r)**p_nn_seed*dexp(-r/r1)
        case('random')
          ! The field at each point (and component) is a Gaussian drawing from
          ! a distribution with variance Bseed
          do var=1,2 !Seed r and phi components of B
            do i=1,nx
              f(i,var)=Bseed(i)*random_normal()
            enddo
          enddo
        case('floor')
          ! Uses floor target field as seed (leave Br=0)
          Bseed = frac_seed*compute_floor_target_field(r,l,h,Beq,delta_r_floor)
          f(:,1)= 0d0
          f(:,2)= Bfloor_sign*Bseed
        case default
          print *, 'init_seed: Invalid option ',p_seed_choice
          stop
      end select
    end subroutine init_seed_field
end module seed_field