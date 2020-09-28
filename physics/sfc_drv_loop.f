!>  \file sfc_drv_loop.f
!!  

!> JP: attempt to bring in surface loop, and sfc_diff into land itself 

      module noah_loop

      use set_soilveg_mod,  only: set_soilveg

      implicit none

      private

      public :: noah_loop_init, noah_loop_run, noah_loop_finalize
      contains

!>\ingroup Noah_LSM
!! This subroutine contains the CCPP-compliant lsm_noah_init to initialize soil vegetation.
!! \section arg_table_lsm_noah_init Argument Table
!! \htmlinclude lsm_noah_init.html
!!
      subroutine noah_loop_init(me, isot, ivegsrc, nlunit,
     & errmsg, errflg)

      implicit none

      integer,              intent(in)  :: me, isot, ivegsrc, nlunit
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!--- initialize soil vegetation
      call set_soilveg(me, isot, ivegsrc, nlunit)

      end subroutine noah_loop_init

      subroutine  noah_loop_finalize
      end subroutine  noah_loop_finalize
      
      subroutine noah_loop_run

      integer :: iter

      do iter=1,2

      enddo


      end subroutine noah_loop_run

      end module noah_loop
