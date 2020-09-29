!>  \file sfc_drv_loop.f
!!  

!> JP: attempt to bring in surface loop, and sfc_diff into land itself 

      module noah_loop

      use set_soilveg_mod,  only: set_soilveg
      use lsm_noah, only: lsm_noah_run
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
      
      subroutine noah_loop_run(                                         &
     &       im, km, grav, cp, hvap, rd, eps, epsm1, rvrdm1, ps,        & !  ---  inputs:
     &       t1, q1, soiltyp, vegtype, sigmaf,                          &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, prslki, zf, land, wind, slopetyp,                   &
     &       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
     &       lheatstrg, isot, ivegsrc,                                  &
     &       bexppert, xlaipert, vegfpert,pertvegf,                     &  ! sfc perts, mgehne
!  ---  in/outs:
     &       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
     &       canopy, trans, tsurf, zorl,                                &
!  ---  outputs:
     &       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
     &       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
     &       smcwlt2, smcref2, wet1, errmsg, errflg                     &
     &     )
 

      use machine,           only: kind_phys
! for lsm_noah
      use funcphys, only : fpvs
      use surface_perturbation, only : ppfbet
      
      implicit none
      
      integer :: iter
      ! Interface variables  
      !integer, intent(in)                                :: im
      !real(kind=kind_phys), dimension(im), intent(in)    :: wind
      logical,              dimension(im), intent(inout) :: flag_guess
      logical,              dimension(im), intent(inout) :: flag_iter

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables                                                                                                               
      integer :: i

! for lsm_noah
!  ---  input:
      integer, intent(in) :: im, km, isot, ivegsrc
      real (kind=kind_phys), intent(in) :: grav, cp, hvap, rd, eps,     &
     &       epsm1, rvrdm1
      real (kind=kind_phys), intent(in) :: pertvegf

      integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &       t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
     &       ch, prsl1, prslki, wind, shdmin, shdmax,                   &
     &       snoalb, sfalb, zf,                                         &
     &       bexppert, xlaipert, vegfpert

      real (kind=kind_phys),  intent(in) :: delt

      logical, dimension(im), intent(in) :: land

      logical, intent(in) :: lheatstrg

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, zorl

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
     &       smc, stc, slc

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
     &       qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,     &
     &       evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2,    &
     &       wet1

      
! Initialize CCPP error handling variables                                                                                      
      errmsg = ''
      errflg = 0

      
      do iter=1,2
         
         ! GFS_surface_loop_control_part1
         do i=1,im
            if (iter == 1 .and. wind(i) < 2.0d0) then
               flag_guess(i) = .true.
            endif
         enddo

         ! GFS_surface_loop_control_part2
         do i = 1, im
            flag_iter(i)  = .false.
            flag_guess(i) = .false.

            if (iter == 1 .and. wind(i) < 2.0d0) then
                  flag_iter(i) = .true.
            endif

         enddo

! sfc_diff

!     lsm_noah
      call lsm_noah_run                                                 &
     &     ( im, km, grav, cp, hvap, rd, eps, epsm1, rvrdm1, ps,        & !  ---  inputs:
     &       t1, q1, soiltyp, vegtype, sigmaf,                          &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, prslki, zf, land, wind, slopetyp,                   &
     &       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
     &       lheatstrg, isot, ivegsrc,                                  &
     &       bexppert, xlaipert, vegfpert,pertvegf,                     &  ! sfc perts, mgehne
!  ---  in/outs:
     &       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
     &       canopy, trans, tsurf, zorl,                                &
!  ---  outputs:
     &       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
     &       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
     &       smcwlt2, smcref2, wet1, errmsg, errflg                     &
     &     )
         
      enddo ! iter


      end subroutine noah_loop_run

      end module noah_loop
