!# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
!#
!# This file is part of Magnetizer.
!#
!# Magnetizer is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# Magnetizer is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
!#
module deriv  !Finite differencing routine for spatial differentiation (same as PENCIL Code)
  use grid
!
  contains
    function xder(f)
!
      implicit none
!
      integer :: ix
      double precision, dimension(nx), intent(in) :: f
      double precision :: fac
      double precision, dimension(nx) :: xder
!
!     assume uniform mesh
!
      fac=1.d0/(60.*(x(2)-x(1)))
!
      do ix=4,nx-3
        xder(ix)=fac*(+45.*(f(ix+1)-f(ix-1)) &
                      - 9.*(f(ix+2)-f(ix-2)) &
                      +    (f(ix+3)-f(ix-3)) &
                     )
      enddo
!
      xder(1)=fac*(-147.d0*f(1)+360.d0*f(2)-450.d0*f(3)+400.d0*f(4)-225.d0*f(5)+72.d0*f(6)-10.d0*f(7))
      xder(2)=fac*( -10.d0*f(1)- 77.d0*f(2)+150.d0*f(3)-100.d0*f(4)+ 50.d0*f(5)-15.d0*f(6)+ 2.d0*f(7))
      xder(3)=fac*(   2.d0*f(1)- 24.d0*f(2)- 35.d0*f(3)+ 80.d0*f(4)- 30.d0*f(5)+ 8.d0*f(6)-      f(7))
!
!     outer points
!
      xder(nx  )=fac*(147.d0*f(nx)-360.d0*f(nx-1)+450.d0*f(nx-2)-400.d0*f(nx-3)+225.d0*f(nx-4)-72.d0*f(nx-5) &
                      +10.d0*f(nx-6))
      xder(nx-1)=fac*( 10.d0*f(nx)+ 77.d0*f(nx-1)-150.d0*f(nx-2)+100.d0*f(nx-3)- 50.d0*f(nx-4)+15.d0*f(nx-5) &
                      - 2.d0*f(nx-6))
      xder(nx-2)=fac*( -2.d0*f(nx)+ 24.d0*f(nx-1)+ 35.d0*f(nx-2)- 80.d0*f(nx-3)+ 30.d0*f(nx-4)- 8.d0*f(nx-5) &
                      +      f(nx-6))
    end function xder
!
    function xder2(f)
!
      implicit none
!
      integer :: ix
      double precision, dimension(nx), intent(in) :: f
      double precision :: fac
      double precision, dimension(nx) :: xder2
!
!     assume uniform mesh
!
      fac=1.d0/(180.d0*(x(2)-x(1))**2)
!
      do ix=4,nx-3
        xder2(ix)=fac*(-490.d0*f(ix)             &
                       +270.d0*(f(ix+1)+f(ix-1)) &
                       - 27.d0*(f(ix+2)+f(ix-2)) &
                       +  2.d0*(f(ix+3)+f(ix-3)) &
                      )
      enddo
!
      xder2(1)=fac*(812.d0*f(1)-3132.d0*f(2)+5265.d0*f(3)-5080.d0*f(4)+2970.d0*f(5)-972.d0*f(6)+137.d0*f(7))
      xder2(2)=fac*(137.d0*f(1)- 147.d0*f(2)- 255.d0*f(3)+ 470.d0*f(4)- 285.d0*f(5)+ 93.d0*f(6)- 13.d0*f(7))
      xder2(3)=fac*(-13.d0*f(1)+ 228.d0*f(2)- 420.d0*f(3)+ 200.d0*f(4)+  15.d0*f(5)- 12.d0*f(6)+  2.d0*f(7))
!
!     outer points
!
      xder2(nx  )=fac*( 812.d0*f(nx)-3132.d0*f(nx-1)+5265.d0*f(nx-2)-5080.d0*f(nx-3)+2970.d0*f(nx-4) &
                       -972.d0*f(nx-5) +137.d0*f(nx-6))
      xder2(nx-1)=fac*( 137.d0*f(nx)- 147.d0*f(nx-1)- 255.d0*f(nx-2)+ 470.d0*f(nx-3)- 285.d0*f(nx-4) &
                       + 93.d0*f(nx-5) - 13.d0*f(nx-6))
      xder2(nx-2)=fac*( -13.d0*f(nx)+ 228.d0*f(nx-1)- 420.d0*f(nx-2)+ 200.d0*f(nx-3)+  15.d0*f(nx-4) &
                        - 12.d0*f(nx-5) + 2.d0*f(nx-6))
    end function xder2
end module deriv
!*****************************************************
