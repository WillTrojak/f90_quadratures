!**************************************************************************
module precision
   use, intrinsic :: iso_fortran_env
   implicit none

   private

   integer, parameter :: int1  = int32
   integer, parameter :: int2  = int64
   integer, parameter :: real1 = real32 
   integer, parameter :: real2 = real64
   integer, parameter :: comp1 = 2*real1
   integer, parameter :: comp2 = 2*real2

   public :: int1,int2,real1,real2,comp1,comp2

contains
   !**********************************************************************

   !**********************************************************************
end module precision
!*************************************************************************
