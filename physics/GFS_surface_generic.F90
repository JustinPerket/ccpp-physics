!> \file GFS_surface_generic.F90
!!  Contains code related to all GFS surface schemes.

!>\defgroup mod_GFS_surface_generic_pre GFS Surface Generic Pre module
      module GFS_surface_generic_pre

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_pre_init, GFS_surface_generic_pre_finalize, GFS_surface_generic_pre_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

      subroutine GFS_surface_generic_pre_init ()
      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! \htmlinclude GFS_surface_generic_pre_run.html
!!
      subroutine GFS_surface_generic_pre_run (im, levs, vfrac, islmsk, isot, ivegsrc, stype, vtype, slope, &
                          prsik_1, prslk_1, tsfc, phil, con_g,                                             &
                          sigmaf, soiltyp, vegtype, slopetyp, work3, zlvl,                                 &
                          drain_cpl, dsnow_cpl, rain_cpl, snow_cpl, lndp_type, n_var_lndp, sfc_wts,        &
                          lndp_var_list, lndp_prt_list,                                                    &
                          sfcemis , dlwflx  , snet    , tg3     , cm      , ch      ,   & ! JP add
                          prsl1   , land    , shdmin  , shdmax  , snoalb  , sfalb   ,   & ! JP add
                          weasd , snwdph, tskin , tprcp , srflag, smc   , stc   , slc   , canopy, trans   ,& ! JP add
                          z0rl  , ustar,                                                                   & ! JP add
                          pgr   , tgrs_1, qgrs_1,                                                          & ! JP add                           
                          z01d, zt1d, bexp1d, xlai1d, vegf1d, lndp_vgf,                                    &
                          cplflx, flag_cice, islmsk_cice, slimskin_cpl,                                    &
                          wind, u1, v1, cnvwind, smcwlt2, smcref2,                                         &
                          soiltyp_cpl , vegtype_cpl , sigmaf_cpl  , sfcemis_cpl , dlwflx_cpl  , snet_cpl    , tg3_cpl     , cm_cpl      , ch_cpl      , & ! JP add
                          prsl1_cpl   , prslki_cpl  , zf_cpl      , land_cpl    , slopetyp_cpl, shdmin_cpl  , shdmax_cpl  , snoalb_cpl  , sfalb_cpl   , & ! JP add
                          bexppert_cpl, xlaipert_cpl, vegfpert_cpl,                                                                                     & ! JP add
                          prsik1_cpl, weasd_cpl , snwdph_cpl, tskin_cpl , tprcp_cpl , srflag_cpl, smc_cpl   , stc_cpl   , slc_cpl   ,                   & ! JP add
                          canopy_cpl, trans_cpl , tsurf_cpl , z0rl_cpl  , z0pert_cpl, ztpert_cpl, ustar_cpl , wind_cpl  ,                               & ! JP add
                          ps_cpl    , t1_cpl    , q1_cpl    ,                                                                                           & ! JP add
                          errmsg, errflg)

        use surface_perturbation,  only: cdfnor

        implicit none

        ! Interface variables
        integer, intent(in) :: im, levs, isot, ivegsrc
        integer, dimension(:), intent(in) :: islmsk
        integer, dimension(:), intent(inout) :: soiltyp, vegtype, slopetyp

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(:), intent(in) :: vfrac, stype, vtype, slope, prsik_1, prslk_1

        real(kind=kind_phys), dimension(:), intent(inout) :: tsfc
        real(kind=kind_phys), dimension(:,:), intent(in) :: phil

        real(kind=kind_phys), dimension(:), intent(inout) :: sigmaf, work3, zlvl

        ! Stochastic physics / surface perturbations
        real(kind=kind_phys), dimension(:),   intent(out) :: drain_cpl
        real(kind=kind_phys), dimension(:),   intent(out) :: dsnow_cpl
        real(kind=kind_phys), dimension(:),   intent(in)  :: rain_cpl
        real(kind=kind_phys), dimension(:),   intent(in)  :: snow_cpl
        integer,                              intent(in)  :: lndp_type, n_var_lndp
        character(len=3),     dimension(:),   intent(in)  :: lndp_var_list
        real(kind=kind_phys), dimension(:),   intent(in)  :: lndp_prt_list
        real(kind=kind_phys), dimension(:,:), intent(in)  :: sfc_wts
        real(kind=kind_phys), dimension(:),   intent(out) :: z01d
        real(kind=kind_phys), dimension(:),   intent(out) :: zt1d
        real(kind=kind_phys), dimension(:),   intent(out) :: bexp1d
        real(kind=kind_phys), dimension(:),   intent(out) :: xlai1d
        real(kind=kind_phys), dimension(:),   intent(out) :: vegf1d
        real(kind=kind_phys),                 intent(out) :: lndp_vgf

        logical,                              intent(in)    :: cplflx
        real(kind=kind_phys), dimension(:),   intent(in)    :: slimskin_cpl
        logical,              dimension(:),   intent(inout) :: flag_cice
        integer,              dimension(:),   intent(out)   :: islmsk_cice

        real(kind=kind_phys), dimension(:),   intent(out) :: wind
        real(kind=kind_phys), dimension(:),   intent(in ) :: u1, v1
        ! surface wind enhancement due to convection
        real(kind=kind_phys), dimension(:),   intent(inout ) :: cnvwind
        !
        real(kind=kind_phys), dimension(:),   intent(out)    :: smcwlt2, smcref2

        ! JP add for export to land
        !integer             , dimension(im),  intent(in)  :: soiltyp
        !integer             , dimension(im),  intent(in)  :: vegtype
        !real(kind=kind_phys), dimension(im),  intent(in)  :: sigmaf
        real(kind=kind_phys), dimension(im),  intent(in)  :: sfcemis
        real(kind=kind_phys), dimension(im),  intent(in)  :: dlwflx
        real(kind=kind_phys), dimension(im),  intent(in)  :: snet
        real(kind=kind_phys), dimension(im),  intent(in)  :: tg3
        real(kind=kind_phys), dimension(im),  intent(in)  :: cm
        real(kind=kind_phys), dimension(im),  intent(in)  :: ch
        real(kind=kind_phys), dimension(im),  intent(in)  :: prsl1
        !real(kind=kind_phys), dimension(im),  intent(in)  :: prslki
        !real(kind=kind_phys), dimension(im),  intent(in)  :: zf
        logical             , dimension(im),  intent(in)  :: land
        !integer             , dimension(im),  intent(in)  :: slopetyp
        real(kind=kind_phys), dimension(im),  intent(in)  :: shdmin
        real(kind=kind_phys), dimension(im),  intent(in)  :: shdmax
        real(kind=kind_phys), dimension(im),  intent(in)  :: snoalb
        real(kind=kind_phys), dimension(im),  intent(in)  :: sfalb
        !real(kind=kind_phys), dimension(im),  intent(in)  :: bexppert
        !real(kind=kind_phys), dimension(im),  intent(in)  :: xlaipert
        !real(kind=kind_phys), dimension(im),  intent(in)  :: vegfpert
        !real(kind=kind_phys), dimension(im),  intent(in)  :: prsik1
        real(kind=kind_phys), dimension(im),  intent(in)  :: weasd
        real(kind=kind_phys), dimension(im),  intent(in)  :: snwdph
        real(kind=kind_phys), dimension(im),  intent(in)  :: tskin
        real(kind=kind_phys), dimension(im),  intent(in)  :: tprcp
        real(kind=kind_phys), dimension(im),  intent(in)  :: srflag
        real(kind=kind_phys), dimension(im),  intent(in)  :: smc
        real(kind=kind_phys), dimension(im),  intent(in)  :: stc
        real(kind=kind_phys), dimension(im),  intent(in)  :: slc
        real(kind=kind_phys), dimension(im),  intent(in)  :: canopy
        real(kind=kind_phys), dimension(im),  intent(in)  :: trans
        !real(kind=kind_phys), dimension(im),  intent(in)  :: tsurf
        real(kind=kind_phys), dimension(im),  intent(in)  :: z0rl
        !real(kind=kind_phys), dimension(im),  intent(in)  :: z0pert
        !real(kind=kind_phys), dimension(im),  intent(in)  :: ztpert
        real(kind=kind_phys), dimension(im),  intent(in)  :: ustar
        real(kind=kind_phys), dimension(im),  intent(in)  :: pgr
        real(kind=kind_phys), dimension(im),  intent(in)  :: tgrs_1
        real(kind=kind_phys), dimension(im),  intent(in)  :: qgrs_1
        
        integer             , dimension(im),  intent(out)  :: soiltyp_cpl
        integer             , dimension(im),  intent(out)  :: vegtype_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: sigmaf_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: sfcemis_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: dlwflx_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: snet_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: tg3_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: cm_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: ch_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: prsl1_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: prslki_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: zf_cpl
        logical             , dimension(im),  intent(out)  :: land_cpl
        integer             , dimension(im),  intent(out)  :: slopetyp_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: shdmin_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: shdmax_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: snoalb_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: sfalb_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: bexppert_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: xlaipert_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: vegfpert_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: prsik1_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: weasd_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: snwdph_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: tskin_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: tprcp_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: srflag_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: smc_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: stc_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: slc_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: canopy_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: trans_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: tsurf_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: z0rl_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: z0pert_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: ztpert_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: ustar_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: wind_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: ps_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: t1_cpl
        real(kind=kind_phys), dimension(im),  intent(out)  :: q1_cpl
        ! JP end
        
        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer              :: i, k
        real(kind=kind_phys) :: onebg, cdfz

        ! Set constants
        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Scale random patterns for surface perturbations with perturbation size
        ! Turn vegetation fraction pattern into percentile pattern
        lndp_vgf=-999.

        if (lndp_type==1) then
          do k =1,n_var_lndp
            select case(lndp_var_list(k))
            case ('rz0')
                z01d(:) = lndp_prt_list(k)* sfc_wts(:,k)
            case ('rzt')
                 zt1d(:) = lndp_prt_list(k)* sfc_wts(:,k)
            case ('shc')
                 bexp1d(:) = lndp_prt_list(k) * sfc_wts(:,k)
            case ('lai')
                xlai1d(:) = lndp_prt_list(k)* sfc_wts(:,k)
            case ('vgf')
        ! note that the pertrubed vegfrac is being used in sfc_drv, but not sfc_diff
              do i=1,im
                call cdfnor(sfc_wts(i,k),cdfz)
                vegf1d(i) = cdfz
              enddo
              lndp_vgf = lndp_prt_list(k)
            end select
          enddo
        endif

        ! End of stochastic physics / surface perturbation

        do i=1,im
          sigmaf(i) = max(vfrac(i), 0.01_kind_phys)
          islmsk_cice(i) = islmsk(i)
          if (islmsk(i) == 2) then
            if (isot == 1) then
              soiltyp(i) = 16
            else
              soiltyp(i) = 9
            endif
            if (ivegsrc == 0 .or. ivegsrc == 4) then
              vegtype(i) = 24
            elseif (ivegsrc == 1) then
              vegtype(i) = 15
            elseif (ivegsrc == 2) then
              vegtype(i) = 13
            elseif (ivegsrc == 3 .or. ivegsrc == 5) then
              vegtype(i) = 15
            endif
            slopetyp(i)  = 9
          else
            soiltyp(i)  = int( stype(i)+0.5_kind_phys )
            vegtype(i)  = int( vtype(i)+0.5_kind_phys )
            slopetyp(i) = int( slope(i)+0.5_kind_phys )    !! clu: slope -> slopetyp
            if (vegtype(i)  < 1) vegtype(i)  = 17
            if (slopetyp(i) < 1) slopetyp(i) = 1
          endif

          work3(i)   = prsik_1(i) / prslk_1(i)

          zlvl(i)    = phil(i,1) * onebg
          smcwlt2(i) = zero
          smcref2(i) = zero

          wind(i)  = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))   &
                         + max(zero, min(cnvwind(i), 30.0_kind_phys)), one)
         !wind(i)  = max(sqrt(Statein%ugrs(i,1)*Statein%ugrs(i,1) + &
         !                         Statein%vgrs(i,1)*Statein%vgrs(i,1))  &
         !              + max(zero, min(Tbd%phy_f2d(i,Model%num_p2d), 30.0)), one)
          cnvwind(i) = zero

        enddo

      if (cplflx) then
        do i=1,im
          islmsk_cice(i) = nint(slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)

          ! ! JP add, for export to land comp
          soiltyp_cpl   (i) = soiltyp(i)
          vegtype_cpl   (i) = vegtype(i)
          sigmaf_cpl    (i) = sigmaf(i)
          !sfcemis_cpl   (i) = sfcemis(i) ! move to composites pre
          !dlwflx_cpl    (i) = dlwflx(i)  ! move to composites inter run
          !dswsfc_cpl    (i) = dswsfc(i)
          snet_cpl      (i) = snet(i)
          tg3_cpl       (i) = tg3(i)
          cm_cpl        (i) = cm(i)
          ch_cpl        (i) = ch(i)
          prsl1_cpl     (i) = prsl1(i)
          prslki_cpl    (i) = work3(i)
          zf_cpl        (i) = zlvl(i)
          !land_cpl      (i) = land(i) ! move to composites pre
          slopetyp_cpl  (i) = slopetyp(i)
          shdmin_cpl    (i) = shdmin(i)
          shdmax_cpl    (i) = shdmax(i)
          snoalb_cpl    (i) = snoalb(i)
          sfalb_cpl     (i) = sfalb(i)
          bexppert_cpl  (i) = bexp1d(i)
          xlaipert_cpl  (i) = xlai1d(i)
          vegfpert_cpl  (i) = vegf1d(i)
          prsik1_cpl    (i) = prsik_1(i)
          !weasd_cpl     (i) = weasd(i)  ! move to composites pre
          !snwdph_cpl    (i) = snwdph(i) ! move to composites pre     
          tskin_cpl     (i) = tskin(i)
          !tprcp_cpl     (i) = tprcp(i) ! move to composites pre
          srflag_cpl    (i) = srflag(i)
          smc_cpl       (i) = smc(i)
          stc_cpl       (i) = stc(i)
          slc_cpl       (i) = slc(i)
          canopy_cpl    (i) = canopy(i)
          trans_cpl     (i) = trans(i)
          !tsurf_cpl     (i) = tsurf(i) ! move to composites pre
          !z0rl_cpl      (i) = z0rl(i)  ! move to composites pre
          z0pert_cpl    (i) = z01d(i)
          ztpert_cpl    (i) = zt1d(i)
          ustar_cpl     (i) = ustar(i)
          wind_cpl      (i) = wind(i)
          ps_cpl        (i) = pgr(i)
          t1_cpl        (i) = tgrs_1(i)
          q1_cpl        (i) = qgrs_1(i)


       enddo
       ! JP tmp
       !write(6,'("sfc_gen_pre: zlvl   - min/max/avg",3g16.6)') minval(zlvl),   maxval(zlvl),   sum(zlvl)/size(zlvl)
      endif

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre


      module GFS_surface_generic_post

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_post_init, GFS_surface_generic_post_finalize, GFS_surface_generic_post_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

      subroutine GFS_surface_generic_post_init ()
      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! \htmlinclude GFS_surface_generic_post_run.html
!!
      subroutine GFS_surface_generic_post_run (im, cplflx, cplchm, cplwav, lssav, dry, icy, wet,                                    &
        dtf, ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1,                                                                            &
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu,           &
        adjvisbmu, adjvisdfu,t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf,              &
        epi, gfluxi, t1, q1, u1, v1, &
        ! soiltyp , vegtype , sigmaf  , sfcemis , dlwflx  , snet    , tg3     , cm      , ch      ,   & ! JP add
        ! prsl1   , prslki  , zf      , land    , slopetyp, shdmin  , shdmax  , snoalb  , sfalb   ,   & ! JP add
        ! bexppert, xlaipert, vegfpert,                                                               & ! JP add
        ! prsik1, weasd , snwdph, tskin , tprcp , srflag, smc   , stc   , slc   , canopy,             & ! JP add
        ! tsurf , z0rl  , z0pert, ztpert, ustar,                                               & ! JP add
        dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       & 
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, ep,                                &
        runoff, srunoff, runof, drain, lheatstrg, h0facu, h0facs, zvfun, hflx, evap, hflxq, hffac, errmsg, errflg)

        implicit none

        integer,                                intent(in) :: im
        logical,                                intent(in) :: cplflx, cplchm, cplwav, lssav
        logical, dimension(:),                  intent(in) :: dry, icy, wet
        real(kind=kind_phys),                   intent(in) :: dtf

        real(kind=kind_phys), dimension(:),  intent(in)  :: ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw,  &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
          t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf
        
        ! ! JP add for export to land
        ! integer             , dimension(im),  intent(in)  :: soiltyp
        ! integer             , dimension(im),  intent(in)  :: vegtype
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: sigmaf
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: sfcemis
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: dlwflx
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: snet
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: tg3
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: cm
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: ch
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: prsl1
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: prslki
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: zf
        ! logical             , dimension(im),  intent(in)  :: land
        ! integer             , dimension(im),  intent(in)  :: slopetyp
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: shdmin
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: shdmax
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: snoalb
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: sfalb
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: bexppert
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: xlaipert
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: vegfpert
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: prsik1
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: weasd
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: snwdph
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: tskin
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: tprcp
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: srflag
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: smc
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: stc
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: slc
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: canopy
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: tsurf
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: z0rl
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: z0pert
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: ztpert
        ! real(kind=kind_phys), dimension(im),  intent(in)  :: ustar
        
        ! integer             , dimension(im),  intent(out)  :: soiltyp_cpl
        ! integer             , dimension(im),  intent(out)  :: vegtype_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: sigmaf_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: sfcemis_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: dlwflx_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: snet_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: tg3_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: cm_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: ch_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: prsl1_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: prslki_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: zf_cpl
        ! logical             , dimension(im),  intent(out)  :: land_cpl
        ! integer             , dimension(im),  intent(out)  :: slopetyp_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: shdmin_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: shdmax_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: snoalb_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: sfalb_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: bexppert_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: xlaipert_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: vegfpert_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: prsik1_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: weasd_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: snwdph_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: tskin_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: tprcp_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: srflag_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: smc_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: stc_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: slc_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: canopy_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: trans_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: tsurf_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: z0rl_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: z0pert_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: ztpert_cpl
        ! real(kind=kind_phys), dimension(im),  intent(out)  :: ustar_cpl
        
        ! ! JP end

        real(kind=kind_phys), dimension(:),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl,        &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl,    &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa,       &
          evcwa, transa, sbsnoa, snowca, snohfa, ep

        real(kind=kind_phys), dimension(:), intent(inout) :: runoff, srunoff
        real(kind=kind_phys), dimension(:), intent(in)    :: drain, runof

        ! For canopy heat storage
        logical, intent(in) :: lheatstrg
        real(kind=kind_phys), intent(in) :: h0facu, h0facs
        real(kind=kind_phys), dimension(:), intent(in)  :: zvfun
        real(kind=kind_phys), dimension(:), intent(in)  :: hflx,  evap
        real(kind=kind_phys), dimension(:), intent(out) :: hflxq
        real(kind=kind_phys), dimension(:), intent(out) :: hffac

        ! CCPP error handling variables
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        real(kind=kind_phys), parameter :: albdf = 0.06_kind_phys

        integer :: i
        real(kind=kind_phys) :: xcosz_loc, ocalnirdf_cpl, ocalnirbm_cpl, ocalvisdf_cpl, ocalvisbm_cpl

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          epi(i)    = ep1d(i)
          gfluxi(i) = gflx(i)
          t1(i)     = tgrs_1(i)
          q1(i)     = qgrs_1(i)
          u1(i)     = ugrs_1(i)
          v1(i)     = vgrs_1(i)
        enddo

        if (cplflx .or. cplchm .or. cplwav) then
          do i=1,im
            u10mi_cpl(i) = u10m(i)
            v10mi_cpl(i) = v10m(i)
          enddo
        endif

        if (cplflx .or. cplchm) then
          do i=1,im
            tsfci_cpl(i) = tsfc(i)
          enddo
        endif

        if (cplflx) then
          do i=1,im
            dlwsfci_cpl (i) = adjsfcdlw(i)
            dswsfci_cpl (i) = adjsfcdsw(i)
            dlwsfc_cpl  (i) = dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
            dswsfc_cpl  (i) = dswsfc_cpl(i) + adjsfcdsw(i)*dtf
            dnirbmi_cpl (i) = adjnirbmd(i)
            dnirdfi_cpl (i) = adjnirdfd(i)
            dvisbmi_cpl (i) = adjvisbmd(i)
            dvisdfi_cpl (i) = adjvisdfd(i)
            dnirbm_cpl  (i) = dnirbm_cpl(i) + adjnirbmd(i)*dtf
            dnirdf_cpl  (i) = dnirdf_cpl(i) + adjnirdfd(i)*dtf
            dvisbm_cpl  (i) = dvisbm_cpl(i) + adjvisbmd(i)*dtf
            dvisdf_cpl  (i) = dvisdf_cpl(i) + adjvisdfd(i)*dtf
            nlwsfci_cpl (i) = adjsfcdlw(i)  - adjsfculw(i)
            if (wet(i)) then
              nlwsfci_cpl(i) = adjsfcdlw(i) - adjsfculw_wat(i)
            endif
            nlwsfc_cpl  (i) = nlwsfc_cpl(i) + nlwsfci_cpl(i)*dtf
            t2mi_cpl    (i) = t2m(i)
            q2mi_cpl    (i) = q2m(i)
            psurfi_cpl  (i) = pgr(i)

            ! ! ! JP add, for export to land comp
            ! soiltyp_cpl   (i) = soiltyp(i)
            ! vegtype_cpl   (i) = vegtype(i)
            ! sigmaf_cpl    (i) = sigmaf(i)
            ! sfcemis_cpl   (i) = sfcemis(i)
            ! dlwflx_cpl    (i) = dlwflx(i)
            ! !dswsfc_cpl    (i) = dswsfc(i)
            ! snet_cpl      (i) = snet(i)
            ! tg3_cpl       (i) = tg3(i)
            ! cm_cpl        (i) = cm(i)
            ! ch_cpl        (i) = ch(i)
            ! prsl1_cpl     (i) = prsl1(i)
            ! prslki_cpl    (i) = prslki(i)
            ! zf_cpl        (i) = zf(i)
            ! land_cpl      (i) = land(i)
            ! slopetyp_cpl  (i) = slopetyp(i)
            ! shdmin_cpl    (i) = shdmin(i)
            ! shdmax_cpl    (i) = shdmax(i)
            ! snoalb_cpl    (i) = snoalb(i)
            ! sfalb_cpl     (i) = sfalb(i)
            ! bexppert_cpl  (i) = bexppert(i)
            ! xlaipert_cpl  (i) = xlaipert(i)
            ! vegfpert_cpl  (i) = vegfpert(i)
            ! prsik1_cpl    (i) = prsik1(i)
            ! weasd_cpl     (i) = weasd(i)
            ! snwdph_cpl    (i) = snwdph(i)
            ! tskin_cpl     (i) = tskin(i)
            ! tprcp_cpl     (i) = tprcp(i)
            ! srflag_cpl    (i) = srflag(i)
            ! smc_cpl       (i) = smc(i)
            ! stc_cpl       (i) = stc(i)
            ! slc_cpl       (i) = slc(i)
            ! canopy_cpl    (i) = canopy(i)
            ! trans_cpl     (i) = trans(i)
            ! tsurf_cpl     (i) = tsurf(i)
            ! z0rl_cpl      (i) = z0rl(i)
            ! z0pert_cpl    (i) = z0pert(i)
            ! ztpert_cpl    (i) = ztpert(i)
            ! ustar_cpl     (i) = ustar(i)
          enddo

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

          do i=1,im
!           if (Sfcprop%landfrac(i) < one) then ! Not 100% land
            if (wet(i)) then                    ! some open water
!  ---  compute open water albedo
              xcosz_loc = max( zero, min( one, xcosz(i) ))
              ocalnirdf_cpl = 0.06_kind_phys
              ocalnirbm_cpl = max(albdf, 0.026_kind_phys/(xcosz_loc**1.7_kind_phys+0.065_kind_phys)     &
       &                       + 0.15_kind_phys * (xcosz_loc-0.1_kind_phys) * (xcosz_loc-0.5_kind_phys) &
       &                       * (xcosz_loc-one))
              ocalvisdf_cpl = 0.06_kind_phys
              ocalvisbm_cpl = ocalnirbm_cpl

              nnirbmi_cpl(i) = adjnirbmd(i) * (one-ocalnirbm_cpl)
              nnirdfi_cpl(i) = adjnirdfd(i) * (one-ocalnirdf_cpl)
              nvisbmi_cpl(i) = adjvisbmd(i) * (one-ocalvisbm_cpl)
              nvisdfi_cpl(i) = adjvisdfd(i) * (one-ocalvisdf_cpl)
            else
              nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
              nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
              nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
              nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
            endif
            nswsfci_cpl(i) = nnirbmi_cpl(i) + nnirdfi_cpl(i)   &
                           + nvisbmi_cpl(i) + nvisdfi_cpl(i)
            nswsfc_cpl(i)  = nswsfc_cpl(i)  + nswsfci_cpl(i)*dtf
            nnirbm_cpl(i)  = nnirbm_cpl(i)  + nnirbmi_cpl(i)*dtf
            nnirdf_cpl(i)  = nnirdf_cpl(i)  + nnirdfi_cpl(i)*dtf
            nvisbm_cpl(i)  = nvisbm_cpl(i)  + nvisbmi_cpl(i)*dtf
            nvisdf_cpl(i)  = nvisdf_cpl(i)  + nvisdfi_cpl(i)*dtf
          enddo
        endif

        if (lssav) then
          do i=1,im
            gflux(i)   = gflux(i)  + gflx(i)  * dtf
            evbsa(i)   = evbsa(i)  + evbs(i)  * dtf
            evcwa(i)   = evcwa(i)  + evcw(i)  * dtf
            transa(i)  = transa(i) + trans(i) * dtf
            sbsnoa(i)  = sbsnoa(i) + sbsno(i) * dtf
            snowca(i)  = snowca(i) + snowc(i) * dtf
            snohfa(i)  = snohfa(i) + snohf(i) * dtf
            ep(i)      = ep(i)     + ep1d(i)  * dtf

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
            runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * dtf
            srunoff(i) = srunoff(i) + runof(i) * dtf
          enddo
        endif

!
!  in order to achieve heat storage within canopy layer, in the canopy
!    heat torage parameterization the kinematic sensible heat flux
!    (hflx) as surface boundary forcing to the pbl scheme is
!    reduced in a factor of hffac given as a function of surface roughness &
!    green vegetation fraction (zvfun) 
!
        do i=1,im
          hflxq(i) = hflx(i)
          hffac(i) = 1.0
        enddo
        if (lheatstrg) then
          do i=1,im
            if (dry(i)) then
              if(hflx(i) > 0.) then
                hffac(i) = h0facu * zvfun(i)
              else
                hffac(i) = h0facs * zvfun(i)
              endif
              hffac(i) = 1. + hffac(i)
              hflxq(i) = hflx(i) / hffac(i)
            endif
          enddo
        endif

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post
