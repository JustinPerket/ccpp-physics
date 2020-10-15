!>  \file sfc_drv_loop.f
!!  

!> JP: attempt to bring in surface loop, and sfc_diff into land itself 

module noah_loop

  implicit none

  private

  public :: noah_loop_init, noah_loop_run, noah_loop_finalize
contains

!>\ingroup Noah_Loop
!! This subroutine contains the CCPP-compliant noah_loop_init to initialize soil vegetation.
!! \section arg_table_noah_loop_init Argument Table
!! \htmlinclude noah_loop_init.html
!!
  subroutine noah_loop_init(me, isot, ivegsrc, nlunit,errmsg, errflg)

    use set_soilveg_mod,  only: set_soilveg
    implicit none

    integer,              intent(in)  :: me, isot, ivegsrc, nlunit
    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

    !     Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !---  initialize soil vegetation
    call set_soilveg(me, isot, ivegsrc, nlunit)

  end subroutine noah_loop_init

!! \section arg_table_noah_loop_finalize Argument Table
!! \htmlinclude noah_loop_finalize.html
!!
  subroutine  noah_loop_finalize
  end subroutine  noah_loop_finalize


!>\defgroup Noah_Loop GFS Loop Driver

!! \section arg_table_noah_loop_run Argument Table
!! \htmlinclude noah_loop_run.html
!!
!> \section general_noah_loop GFS sfc_drv_loop General Algorithm
!>  @{
  subroutine noah_loop_run(                                         &
                                !! ARGS FROM NOAH      
                                !  ---  inputs:
       im, km, grav, cp, hvap, rd, eps, epsm1, rvrdm1, ps,          & 
       t1, q1, soiltyp, vegtype, sigmaf,                            &
       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,            &
       prsl1, prslki, zf, land, wind, slopetyp,                     &
       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,        &
       lheatstrg, isot, ivegsrc,                                    &
       bexppert, xlaipert, vegfpert,pertvegf,                       & ! sfc perts, mgehne
                                !     ---  in/outs:
       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,          &
       canopy, trans, tsurf, z0rl,                                  &
                                !     ---  outputs:
       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,        &
       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,              &
       smcwlt2, smcref2, wet1,                                      &
                                !! ARGS FROM  stab_prep_lnd (minus those from noah
                                !  ---  inputs:
       dry,prsik1,z0pert,ztpert,ustar_lnd,                           &
                                !  ---  outputs:
                                !! ARGS FROM stability (minus those above)
                                !  ---  inputs:
                                !  ---  outputs:
       rb_lnd, fm_lnd, fh_lnd, fm10_lnd, fh2_lnd,                   &
       stress_lnd,                                                  &           
                                !!
       errmsg, errflg                                               &
       )

    use lsm_noah, only: lsm_noah_run
    use sfc_diff, only: stab_prep_lnd, stability

    use machine,           only: kind_phys
    !     for lsm_noah
    use funcphys, only : fpvs
    use surface_perturbation, only : ppfbet

    implicit none

    integer :: iter
    ! Interface variables  
    !integer, intent(in)                                :: im
    !real(kind=kind_phys), dimension(im), intent(in)    :: wind
    logical,              dimension(im), intent(inout) :: flag_guess
    logical,              dimension(im), intent(inout) :: flag_iter
    logical,              dimension(im), intent(in) :: dry
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables                                                                                                               
    integer :: i

!!!     for lsm_noah
    !     ---  input:
    integer, intent(in) :: im, km, isot, ivegsrc
    real (kind=kind_phys), intent(in) :: grav, cp, hvap, rd, eps, epsm1, rvrdm1
    real (kind=kind_phys), intent(in) :: pertvegf

    integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

    real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
         t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,      &
         ch, prsl1, prslki, wind, shdmin, shdmax,                     &
         snoalb, sfalb, zf,                                           &
         bexppert, xlaipert, vegfpert

    real (kind=kind_phys),  intent(in) :: delt

    logical, dimension(im), intent(in) :: land

    logical, intent(in) :: lheatstrg

    !     ---  in/out:
    real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
         snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, z0rl

    real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
         smc, stc, slc

    !     ---  output:
    real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
         qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,       &
         evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2,      &
         wet1

!!!     for stab_prep_lnd ----------------------------------------------
    integer, parameter  :: kp = kind_phys
    !  ---  inputs:
    real(kind=kind_phys),  dimension(im), intent(in) :: prsik1,    &
         z0pert,ztpert,ustar_lnd
    !  ---  outputs:

    !  ---  locals:
    real(kind=kind_phys) :: tem1,tem2,czilc,thv1
    real(kind=kind_phys), parameter :: one=1.0_kp, zero=0.0_kp, half=0.5_kp, &
         zmin=1.0e-6_kp, log01=log(0.01_kp), log05=log(0.05_kp), log07=log(0.07_kp)
    real(kind=kind_phys) :: tvs, z0max, ztmax, virtfac

!!!     for stability      ----------------------------------------------
    real(kind=kind_phys), dimension(im), intent(inout) :: rb_lnd, fm_lnd,    &
         fh_lnd, fm10_lnd, fh2_lnd,stress_lnd


    real(kind=kind_phys), parameter :: qmin=1.0e-8_kp

    !     Initialize CCPP error handling variables                                                                                      
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

       !     sfc_diff
       do i=1,im
          if(flag_iter(i)) then
             virtfac = one + rvrdm1 * max(q1(i),qmin)
             thv1    = t1(i) * prslki(i) * virtfac

             if (dry(i)) then ! Some land

                call stab_prep_lnd                                     &
                                !     ---  inputs:
                     (zf(i),prsik1(i),sigmaf(i),vegtype(i),shdmax(i),  &
                     ivegsrc,z0pert(i),ztpert(i),                      &
                     tskin(i),tsurf(i),z0rl(i),            &
                     ustar_lnd(i),virtfac,                             &
                                !     ---  outputs:
                     z0max,ztmax,tvs   ) 

                call stability                                        &
                                !     ---  inputs:
                     (zf(i), snwdph(i), thv1, wind(i),                &
                     z0max, ztmax, tvs, grav,                         &
                                !     ---  outputs:
                     rb_lnd(i), fm_lnd(i), fh_lnd(i), fm10_lnd(i),    &
                     fh2_lnd(i),cm_lnd(i), ch_lnd(i), stress_lnd(i),  &
                     ustar_lnd(i)  )

                !     lsm_noah
                call lsm_noah_run                                               &
                     ( im, km, grav, cp, hvap, rd, eps, epsm1, rvrdm1, ps,      & !  ---  inputs:
                     t1, q1, soiltyp, vegtype, sigmaf,                          &
                     sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
                     prsl1, prslki, zf, land, wind, slopetyp,                   &
                     shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
                     lheatstrg, isot, ivegsrc,                                  &
                     bexppert, xlaipert, vegfpert,pertvegf,                     & ! sfc perts, mgehne
                                !     ---  in/outs:
                     weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
                     canopy, trans, tsurf, z0rl,                                &
                                !     ---  outputs:
                     sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
                     cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
                     smcwlt2, smcref2, wet1, errmsg, errflg                     &
                     )
             endif            ! Dry points   
          endif ! flag_iter
       enddo ! im
    enddo ! iter


  end subroutine noah_loop_run
!-----------------------------------
!> @}
end module noah_loop
