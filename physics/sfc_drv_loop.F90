
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
    !write(*,*) 'JP in noah_loop_init'
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
       prsik1,z0pert,ztpert,ustar,                           &
                                !  ---  outputs:
                                !! ARGS FROM stability (minus those above)
                                !  ---  inputs:
                                !  ---  outputs:
       rb_lnd, fm_lnd, fh_lnd, fm10_lnd, fh2_lnd,                   &
       stress,                                                  &           
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
         t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3,          &
         prsl1, prslki, wind, shdmin, shdmax,                     &
         snoalb, sfalb, zf,                                           &
         bexppert, xlaipert, vegfpert

    real (kind=kind_phys),  intent(in) :: delt

    logical, dimension(im), intent(in) :: land

    logical, intent(in) :: lheatstrg

    !     ---  in/out:
    real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
         snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, z0rl,    &
         cm, ch
    
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
         z0pert,ztpert
    !  ---  outputs:

    !  ---  locals:
    real(kind=kind_phys) :: tem1,tem2,czilc,thv1
    real(kind=kind_phys), parameter :: one=1.0_kp, zero=0.0_kp, half=0.5_kp, &
         zmin=1.0e-6_kp, log01=log(0.01_kp), log05=log(0.05_kp), log07=log(0.07_kp)
    real(kind=kind_phys) :: tvs, z0max, ztmax, virtfac

!!!     for stability      ----------------------------------------------
    real(kind=kind_phys), dimension(im), intent(inout) :: rb_lnd, fm_lnd,    &
         fh_lnd, fm10_lnd, fh2_lnd,stress,ustar


    real(kind=kind_phys), parameter :: qmin=1.0e-8_kp

    !write(*,*) 'JP in noah_loop_run 1'

    !     Initialize CCPP error handling variables                                                                                      
    errmsg = ''
    errflg = 0

    ! revert
    do i=1,im
       flag_guess(i)   = .false.
       flag_iter(i)    = .true.
    end do
       
    write(6,'("sfc_drv_loop: zf   - min/max/avg",3g16.6)') minval(zf),   maxval(zf),   sum(zf)/size(zf)
    write(6,'("sfc_drv_loop: prsik1   - min/max/avg",3g16.6)') minval(prsik1),   maxval(prsik1),   sum(prsik1)/size(prsik1)
    write(6,'("sfc_drv_loop: sigmaf   - min/max/avg",3g16.6)') minval(sigmaf),   maxval(sigmaf),   sum(sigmaf)/size(sigmaf)
    write(6,'("sfc_drv_loop: vegtype   - min/max/avg",3g16.6)') minval(vegtype),   maxval(vegtype),   sum(vegtype)/size(vegtype)
    write(6,'("sfc_drv_loop: shdmax   - min/max/avg",3g16.6)') minval(shdmax),   maxval(shdmax),   sum(shdmax)/size(shdmax)
    write(6,'("sfc_drv_loop: z0pert   - min/max/avg",3g16.6)') minval(z0pert),   maxval(z0pert),   sum(z0pert)/size(z0pert)
    write(6,'("sfc_drv_loop: tskin   - min/max/avg",3g16.6)') minval(tskin),   maxval(tskin),   sum(tskin)/size(tskin)
    write(6,'("sfc_drv_loop: tsurf   - min/max/avg",3g16.6)') minval(tsurf),   maxval(tsurf),   sum(tsurf)/size(tsurf)
    write(6,'("sfc_drv_loop: z0rl   - min/max/avg",3g16.6)') minval(z0rl),   maxval(z0rl),   sum(z0rl)/size(z0rl)
    write(6,'("sfc_drv_loop: ustar   - min/max/avg",3g16.6)') minval(ustar),   maxval(ustar),   sum(ustar)/size(ustar)
    write(6,'("sfc_drv_loop: snwdph   - min/max/avg",3g16.6)') minval(snwdph),   maxval(snwdph),   sum(snwdph)/size(snwdph)
    write(6,'("sfc_drv_loop: wind   - min/max/avg",3g16.6)') minval(wind),   maxval(wind),   sum(wind)/size(wind)

    do iter=1,2

       
       !!!     sfc_diff
       do i=1,im
          if(flag_iter(i)) then
             virtfac = one + rvrdm1 * max(q1(i),qmin)
             thv1    = t1(i) * prslki(i) * virtfac

             if (land(i)) then ! Some land
             !! don't think land is exporting from atm after move to sfc generic pre. test workaround:
             !if (soiltyp(i) < 14.0) then ! Some land

                call stab_prep_lnd                                     &
                                !     ---  inputs:
                     (zf(i),prsik1(i),sigmaf(i),vegtype(i),shdmax(i),  &
                     ivegsrc,z0pert(i),ztpert(i),                      &
                     tskin(i),tsurf(i),z0rl(i),                        &
                     ustar(i),virtfac,                                 &
                                !     ---  outputs:
                     z0max,ztmax,tvs   ) 
                
                call stability                                        &
                                !     ---  inputs:
                     (zf(i), snwdph(i), thv1, wind(i),                &
                     z0max, ztmax, tvs, grav,                         &
                                !     ---  outputs:
                     rb_lnd(i), fm_lnd(i), fh_lnd(i), fm10_lnd(i),    &
                     fh2_lnd(i),cm(i), ch(i), stress(i),              &
                     ustar(i)  )
             end if ! land
          end if ! flag_iter
       end do ! im


       write(6,'("sfc_drv_loop 2: rb_lnd   - min/max/avg",3g16.6)') minval(rb_lnd),   maxval(rb_lnd),   sum(rb_lnd)/size(rb_lnd)
       write(6,'("sfc_drv_loop 2: fm_lnd   - min/max/avg",3g16.6)') minval(fm_lnd),   maxval(fm_lnd),   sum(fm_lnd)/size(fm_lnd)
       write(6,'("sfc_drv_loop 2: fh_lnd   - min/max/avg",3g16.6)') minval(fh_lnd),   maxval(fh_lnd),   sum(fh_lnd)/size(fh_lnd)
       write(6,'("sfc_drv_loop 2: fm10_lnd   - min/max/avg",3g16.6)') minval(fm10_lnd),   maxval(fm10_lnd),   sum(fm10_lnd)/size(fm10_lnd)
       write(6,'("sfc_drv_loop 2: fh2_lnd   - min/max/avg",3g16.6)') minval(fh2_lnd),   maxval(fh2_lnd),   sum(fh2_lnd)/size(fh2_lnd)
       write(6,'("sfc_drv_loop 2: cm   - min/max/avg",3g16.6)') minval(cm),   maxval(cm),   sum(cm)/size(cm)
       write(6,'("sfc_drv_loop 2: ch   - min/max/avg",3g16.6)') minval(ch),   maxval(ch),   sum(ch)/size(ch)
       write(6,'("sfc_drv_loop 2: stress   - min/max/avg",3g16.6)') minval(stress),   maxval(stress),   sum(stress)/size(stress)
       write(6,'("sfc_drv_loop 2: ustar   - min/max/avg",3g16.6)') minval(ustar),   maxval(ustar),   sum(ustar)/size(ustar)

       ! JP end
       
 !!! CCPP  GFS_surface_loop_control_part1
       do i=1,im
          if (iter == 1 .and. wind(i) < 2.0d0) then
             flag_guess(i) = .true.
          endif
       enddo

!!! Land Call
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
       
!!! GFS_surface_loop_control_part2
       do i = 1, im
          flag_iter(i)  = .false.
          flag_guess(i) = .false.

          if (iter == 1 .and. wind(i) < 2.0d0) then
             if (land(i)) then
                flag_iter(i) = .true.
             end if
          endif
       enddo

    enddo ! iter

    ! tmp debug
    write(6,*) "lsm_run: im", im
    write(6,*) "lsm_run: km", km
    write(6,*) "lsm_run: grav", grav
    write(6,*) "lsm_run: cp", cp
    write(6,*) "lsm_run: hvap", hvap
    write(6,*) "lsm_run: rd", rd
    write(6,*) "lsm_run: eps", eps
    write(6,*) "lsm_run: epsm1", epsm1
    write(6,*) "lsm_run: rvrdm1", rvrdm1
    write(6,*) "lsm_run: delt", delt
    write(6,*) "lsm_run: rd", rd
    write(6,*) "lsm_run: rd", rd
    
    write(6,'("lsm_run: ps   - min/max/avg",3g16.6)') minval(ps),   maxval(ps),   sum(ps)/size(ps)
    write(6,'("lsm_run: t1   - min/max/avg",3g16.6)') minval(t1),   maxval(t1),   sum(t1)/size(t1)
    write(6,'("lsm_run: q1   - min/max/avg",3g16.6)') minval(q1),   maxval(q1),   sum(q1)/size(q1)
    write(6,'("lsm_run: soiltyp   - min/max/avg",3g16.6)') minval(soiltyp),   maxval(soiltyp),   sum(soiltyp)/size(soiltyp)
    write(6,'("lsm_run: vegtype   - min/max/avg",3g16.6)') minval(vegtype),   maxval(vegtype),   sum(vegtype)/size(vegtype)
    write(6,'("lsm_run: sigmaf   - min/max/avg",3g16.6)') minval(sigmaf),   maxval(sigmaf),   sum(sigmaf)/size(sigmaf)
    write(6,'("lsm_run: sfcemis   - min/max/avg",3g16.6)') minval(sfcemis),   maxval(sfcemis),   sum(sfcemis)/size(sfcemis)
    write(6,'("lsm_run: dlwflx   - min/max/avg",3g16.6)') minval(dlwflx),   maxval(dlwflx),   sum(dlwflx)/size(dlwflx)
    write(6,'("lsm_run: dswsfc   - min/max/avg",3g16.6)') minval(dswsfc),   maxval(dswsfc),   sum(dswsfc)/size(dswsfc)
    write(6,'("lsm_run: snet   - min/max/avg",3g16.6)') minval(snet),   maxval(snet),   sum(snet)/size(snet)
    write(6,'("lsm_run: tg3   - min/max/avg",3g16.6)') minval(tg3),   maxval(tg3),   sum(tg3)/size(tg3)
    write(6,'("lsm_run: cm   - min/max/avg",3g16.6)') minval(cm),   maxval(cm),   sum(cm)/size(cm)
    write(6,'("lsm_run: ch   - min/max/avg",3g16.6)') minval(ch),   maxval(ch),   sum(ch)/size(ch)
    write(6,'("lsm_run: prsl1   - min/max/avg",3g16.6)') minval(prsl1),   maxval(prsl1),   sum(prsl1)/size(prsl1)
    write(6,'("lsm_run: prslki   - min/max/avg",3g16.6)') minval(prslki),   maxval(prslki),   sum(prslki)/size(prslki)
    write(6,'("lsm_run: zf   - min/max/avg",3g16.6)') minval(zf),   maxval(zf),   sum(zf)/size(zf)
    !write(6,'("lsm_run: land   - min/max/avg",3g16.6)') minval(land),   maxval(land),   sum(land)/size(land)
    write(6,'("lsm_run: wind   - min/max/avg",3g16.6)') minval(wind),   maxval(wind),   sum(wind)/size(wind)
    write(6,'("lsm_run: slopetyp   - min/max/avg",3g16.6)') minval(slopetyp),   maxval(slopetyp),   sum(slopetyp)/size(slopetyp)
    write(6,'("lsm_run: shdmin   - min/max/avg",3g16.6)') minval(shdmin),   maxval(shdmin),   sum(shdmin)/size(shdmin)
    write(6,'("lsm_run: shdmax   - min/max/avg",3g16.6)') minval(shdmax),   maxval(shdmax),   sum(shdmax)/size(shdmax)
    write(6,'("lsm_run: snoalb   - min/max/avg",3g16.6)') minval(snoalb),   maxval(snoalb),   sum(snoalb)/size(snoalb)
    write(6,'("lsm_run: sfalb   - min/max/avg",3g16.6)') minval(sfalb),   maxval(sfalb),   sum(sfalb)/size(sfalb)
    !write(6,'("lsm_run: lheatstrg   - min/max/avg",3g16.6)') minval(lheatstrg),   maxval(lheatstrg),   sum(lheatstrg)/size(lheatstrg)
    !write(6,'("lsm_run: isot   - min/max/avg",3g16.6)') minval(isot),   maxval(isot),   sum(isot)/size(isot)
    !write(6,'("lsm_run: ivegsrc   - min/max/avg",3g16.6)') minval(ivegsrc),   maxval(ivegsrc),   sum(ivegsrc)/size(ivegsrc)
    !write(6,'("lsm_run: bexppert   - min/max/avg",3g16.6)') minval(bexppert),   maxval(bexppert),   sum(bexppert)/size(bexppert)
    !write(6,'("lsm_run: xlaipert   - min/max/avg",3g16.6)') minval(xlaipert),   maxval(xlaipert),   sum(xlaipert)/size(xlaipert)
    !write(6,'("lsm_run: vegfpert   - min/max/avg",3g16.6)') minval(vegfpert),   maxval(vegfpert),   sum(vegfpert)/size(vegfpert)
    !write(6,'("lsm_run: pertvegf   - min/max/avg",3g16.6)') minval(pertvegf),   maxval(pertvegf),   sum(pertvegf)/size(pertvegf)
    write(6,'("lsm_run 2: weasd   - min/max/avg",3g16.6)') minval(weasd),   maxval(weasd),   sum(weasd)/size(weasd)
    write(6,'("lsm_run 2: snwdph   - min/max/avg",3g16.6)') minval(snwdph),   maxval(snwdph),   sum(snwdph)/size(snwdph)
    write(6,'("lsm_run 2: tskin   - min/max/avg",3g16.6)') minval(tskin),   maxval(tskin),   sum(tskin)/size(tskin)
    write(6,'("lsm_run 2: tprcp   - min/max/avg",3g16.6)') minval(tprcp),   maxval(tprcp),   sum(tprcp)/size(tprcp)
    !write(6,'("lsm_run 2: srflag   - min/max/avg",3g16.6)') minval(srflag),   maxval(srflag),   sum(srflag)/size(srflag)
    write(6,'("lsm_run 2: smc   - min/max/avg",3g16.6)') minval(smc),   maxval(smc),   sum(smc)/size(smc)
    write(6,'("lsm_run 2: stc   - min/max/avg",3g16.6)') minval(stc),   maxval(stc),   sum(stc)/size(stc)
    write(6,'("lsm_run 2: slc   - min/max/avg",3g16.6)') minval(slc),   maxval(slc),   sum(slc)/size(slc)
    write(6,'("lsm_run 2: canopy   - min/max/avg",3g16.6)') minval(canopy),   maxval(canopy),   sum(canopy)/size(canopy)
    write(6,'("lsm_run 2: trans   - min/max/avg",3g16.6)') minval(trans),   maxval(trans),   sum(trans)/size(trans)
    write(6,'("lsm_run 2: tsurf   - min/max/avg",3g16.6)') minval(tsurf),   maxval(tsurf),   sum(tsurf)/size(tsurf)
    write(6,'("lsm_run 2: z0rl   - min/max/avg",3g16.6)') minval(z0rl),   maxval(z0rl),   sum(z0rl)/size(z0rl)
    write(6,'("lsm_run 3: sncovr1   - min/max/avg",3g16.6)') minval(sncovr1),   maxval(sncovr1),   sum(sncovr1)/size(sncovr1)
    write(6,'("lsm_run 3: qsurf   - min/max/avg",3g16.6)') minval(qsurf),   maxval(qsurf),   sum(qsurf)/size(qsurf)
    write(6,'("lsm_run 3: gflux   - min/max/avg",3g16.6)') minval(gflux),   maxval(gflux),   sum(gflux)/size(gflux)
    write(6,'("lsm_run 3: drain   - min/max/avg",3g16.6)') minval(drain),   maxval(drain),   sum(drain)/size(drain)
    write(6,'("lsm_run 3: evap   - min/max/avg",3g16.6)') minval(evap),   maxval(evap),   sum(evap)/size(evap)
    write(6,'("lsm_run 3: hflx   - min/max/avg",3g16.6)') minval(hflx),   maxval(hflx),   sum(hflx)/size(hflx)
    write(6,'("lsm_run 3: ep   - min/max/avg",3g16.6)') minval(ep),   maxval(ep),   sum(ep)/size(ep)
    write(6,'("lsm_run 3: runoff   - min/max/avg",3g16.6)') minval(runoff),   maxval(runoff),   sum(runoff)/size(runoff)
    write(6,'("lsm_run 3: cmm   - min/max/avg",3g16.6)') minval(cmm),   maxval(cmm),   sum(cmm)/size(cmm)
    write(6,'("lsm_run 3: chh   - min/max/avg",3g16.6)') minval(chh),   maxval(chh),   sum(chh)/size(chh)
    write(6,'("lsm_run 3: evbs   - min/max/avg",3g16.6)') minval(evbs),   maxval(evbs),   sum(evbs)/size(evbs)
    write(6,'("lsm_run 3: evcw   - min/max/avg",3g16.6)') minval(evcw),   maxval(evcw),   sum(evcw)/size(evcw)
    write(6,'("lsm_run 3: sbsno   - min/max/avg",3g16.6)') minval(sbsno),   maxval(sbsno),   sum(sbsno)/size(sbsno)
    write(6,'("lsm_run 3: snowc   - min/max/avg",3g16.6)') minval(snowc),   maxval(snowc),   sum(snowc)/size(snowc)
    write(6,'("lsm_run 3: stm   - min/max/avg",3g16.6)') minval(stm),   maxval(stm),   sum(stm)/size(stm)
    write(6,'("lsm_run 3: snohf   - min/max/avg",3g16.6)') minval(snohf),   maxval(snohf),   sum(snohf)/size(snohf)
    write(6,'("lsm_run 3: smcwlt2   - min/max/avg",3g16.6)') minval(smcwlt2),   maxval(smcwlt2),   sum(smcwlt2)/size(smcwlt2)
    write(6,'("lsm_run 3: smcref2   - min/max/avg",3g16.6)') minval(smcref2),   maxval(smcref2),   sum(smcref2)/size(smcref2)
    write(6,'("lsm_run 3: wet1   - min/max/avg",3g16.6)') minval(wet1),   maxval(wet1),   sum(wet1)/size(wet1)

  end subroutine noah_loop_run
!-----------------------------------
!> @}
end module noah_loop
