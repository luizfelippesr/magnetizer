!# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
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
!! (Slightly) adapted from Galacticus (https://bitbucket.org/abensonca/galacticus/)
!!
!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which performs numerical integration.

module Integration
  !% Implements numerical integration.
  use FGSL
  implicit none
  private
  public :: Integrate_Done, Integrate

  ! Module scope error status.
  integer :: errorStatusGlobal

  ! Integrand interface.
  abstract interface
     double precision function integrandTemplate(x)
       double precision, intent(in   ) :: x
     end function integrandTemplate
  end interface

  ! Integrand function.
  procedure(integrandTemplate), pointer :: currentIntegrand
  !$omp threadprivate(currentIntegrand)
  
contains

  recursive double precision function Integrate(lowerLimit,upperLimit,integrand,integrandFunction, &
       integrationWorkspace,maxIntervals,toleranceAbsolute,toleranceRelative,hasSingularities, &
       fromInfinity, toInfinity, integrationRule,reset,errorStatus)
    !% Integrates the supplied {\normalfont \ttfamily integrand} function.
    use, intrinsic :: ISO_C_Binding
    implicit none
    double precision, intent(in) :: lowerLimit, upperLimit
    procedure(integrandTemplate) :: integrand
    type(fgsl_function), intent(inout) :: integrandFunction
    type(fgsl_integration_workspace), intent(inout) :: integrationWorkspace
    type(fgsl_error_handler_t)  :: integrationErrorHandler, standardGslErrorHandler
    integer, intent(in   ), optional :: integrationRule, maxIntervals
    double precision, intent(in   ), optional :: toleranceAbsolute, toleranceRelative
    logical, intent(in   ), optional :: hasSingularities, toInfinity, fromInfinity
    logical, intent(inout), optional :: reset
    integer, intent(  out), optional :: errorStatus
    integer, parameter :: maxIntervalsDefault = 1000
    double precision, parameter :: toleranceAbsoluteDefault = 1.0d-10, toleranceRelativeDefault = 1.0d-10
    type(c_ptr                     )                                     :: parameterPointer
    integer :: integrationRuleActual, status
    integer(kind=c_size_t)  :: maxIntervalsActual
    double precision :: integrationError, integrationValue, toleranceAbsoluteActual, toleranceRelativeActual
    logical :: hasSingularitiesActual, resetActual, fromInfinityActual, toInfinityActual
    procedure(integrandTemplate), pointer :: previousIntegrand
    
    ! Store the current integrand function so that we can restore it on exit. This allows the integration function to be called recursively.
    previousIntegrand => currentIntegrand
    currentIntegrand  => integrand
    ! Set optional parameters to specified or default values.
    if (present(maxIntervals)) then
       maxIntervalsActual=maxIntervals
    else
       maxIntervalsActual=maxIntervalsDefault
    end if
    if (present(toleranceAbsolute)) then
       toleranceAbsoluteActual=toleranceAbsolute
    else
       toleranceAbsoluteActual=toleranceAbsoluteDefault
    end if
    if (present(toleranceRelative)) then
       toleranceRelativeActual=toleranceRelative
    else
       toleranceRelativeActual=toleranceRelativeDefault
    end if
    if (present(hasSingularities)) then
       hasSingularitiesActual=hasSingularities
    else
       hasSingularitiesActual=.false.
    end if
    if (present(toInfinity)) then
       toInfinityActual=toInfinity
    else
       toInfinityActual=.false.
    end if
    if (present(fromInfinity)) then
       fromInfinityActual=fromInfinity
    else
       fromInfinityActual=.false.
    end if
    if (present(integrationRule)) then
       integrationRuleActual=integrationRule
    else
       integrationRuleActual=FGSL_Integ_Gauss61
    end if

    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.true.
    end if

    ! Initialize the integration variables if necessary.
    if (resetActual) then
       integrationWorkspace=FGSL_Integration_Workspace_Alloc(maxIntervalsActual)
       integrandFunction   =FGSL_Function_Init(integrandWrapper,parameterPointer)
    end if

    ! Set error handler if necessary.
    if (present(errorStatus)) then
       integrationErrorHandler=FGSL_Error_Handler_Init(Integration_GSL_Error_Handler)
       standardGslErrorHandler=FGSL_Set_Error_Handler (integrationErrorHandler      )
       errorStatusGlobal=FGSL_Success
    end if

    ! Do the integration
    if ((.not.fromInfinityActual) .and. (.not.toInfinityActual)) then ! A -> B
      if (hasSingularitiesActual) then
        status=FGSL_Integration_QAGS(integrandFunction,lowerLimit,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
              &,maxIntervalsActual,integrationWorkspace,integrationValue,integrationError)
      else
        status=FGSL_Integration_QAG(integrandFunction,lowerLimit,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
              &,maxIntervalsActual,integrationRuleActual,integrationWorkspace,integrationValue,integrationError)
      end if
    else if (fromInfinityActual .and. toInfinityActual) then ! -\infty -> \infty
        status=FGSL_Integration_QAGI(integrandFunction,toleranceAbsoluteActual,toleranceRelativeActual &
              &,maxIntervalsActual,integrationWorkspace,integrationValue,integrationError)
    else if (fromInfinityActual .and. .not.toInfinityActual) then ! -\infty -> B
        status=FGSL_Integration_QAGIL(integrandFunction,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
              &,maxIntervalsActual,integrationWorkspace,integrationValue,integrationError)
    else if (.not.fromInfinityActual .and. toInfinityActual) then ! A -> \infty
        status=FGSL_Integration_QAGIU(integrandFunction,lowerLimit,toleranceAbsoluteActual,toleranceRelativeActual &
              &,maxIntervalsActual,integrationWorkspace,integrationValue,integrationError)
    else
      stop 'Something weird happened'
    endif
    Integrate=integrationValue

    ! Reset error handler.
    if (present(errorStatus)) then
       errorStatus            =errorStatusGlobal
       standardGslErrorHandler=FGSL_Set_Error_Handler (standardGslErrorHandler)
    end if
    ! Restore the previous integrand.
    currentIntegrand => previousIntegrand
    return
  end function Integrate

  function integrandWrapper(x,parameterPointer) bind(c)
    !% Wrapper function used for GSL integration functions.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)        :: integrandWrapper
    real(c_double), value :: x
    type(c_ptr   ), value :: parameterPointer

    integrandWrapper=currentIntegrand(x)
    return
  end function integrandWrapper
  
  subroutine Integrate_Done(integrandFunction,integrationWorkspace)
    !% Frees up integration objects that are no longer required.
    implicit none
    type(fgsl_function             ), intent(inout) :: integrandFunction
    type(fgsl_integration_workspace), intent(inout) :: integrationWorkspace

    call FGSL_Function_Free(integrandFunction)
    call FGSL_Integration_Workspace_Free(integrationWorkspace)
    return
  end subroutine Integrate_Done

  subroutine Integration_GSL_Error_Handler(reason,file,line,errorNumber) bind(c)
    !% Handle errors from the GSL library during integration.
    use, intrinsic :: ISO_C_Binding
    type   (c_ptr     ), value :: file       , reason
    integer(kind=c_int), value :: errorNumber, line

    errorStatusGlobal=errorNumber
    return
  end subroutine Integration_GSL_Error_Handler

end module Integration
