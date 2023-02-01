!> \file GFS_surface_generic_pre.F90
!!  Contains code related to running prior to all GFS surface schemes.

      module GFS_surface_generic_pre

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_pre_init, GFS_surface_generic_pre_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

!>\defgroup mod_GFS_surface_generic_pre GFS surface_generic_pre module
!! This module contains code related to running prior to all GFS surface schemes.
!> @{
!> \section arg_table_GFS_surface_generic_pre_init Argument Table
!! \htmlinclude GFS_surface_generic_pre_init.html
!!
      subroutine GFS_surface_generic_pre_init (nthreads, im, slmsk, isot, ivegsrc, stype, vtype, slope, &
                                               vtype_save, stype_save, slope_save, errmsg, errflg)

        implicit none

        ! Interface variables
        integer,                       intent(in)    :: nthreads, im, isot, ivegsrc
        real(kind_phys), dimension(:), intent(in)    :: slmsk
        integer,         dimension(:), intent(inout) :: vtype, stype, slope
        integer,         dimension(:), intent(out)   :: vtype_save, stype_save, slope_save

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer, dimension(1:im) :: islmsk
        integer :: i

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        islmsk = nint(slmsk)

        ! Save current values of vegetation, soil and slope type
        vtype_save(:) = vtype(:)
        stype_save(:) = stype(:)
        slope_save(:) = slope(:)

        call update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype, slope)

      end subroutine GFS_surface_generic_pre_init

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! \htmlinclude GFS_surface_generic_pre_run.html
!!
      subroutine GFS_surface_generic_pre_run (nthreads, im, levs, vfrac, islmsk, isot, ivegsrc, stype, vtype, slope, &
                          prsik_1, prslk_1, tsfc, phil, con_g, sigmaf, work3, zlvl,                        &
                          drain_cpl, dsnow_cpl, rain_cpl, snow_cpl, lndp_type, n_var_lndp, sfc_wts,        &
                          lndp_var_list, lndp_prt_list,                                                    &
                          sfcemis , dlwflx  , snet    , tg3     , cm      , ch      ,   & ! JP add
                          prsl1   , land    , shdmin  , shdmax  , snoalb  , sfalb   ,   & ! JP add
                          weasd , snwdph, tskin , tprcp , srflag, smc   , stc   , slc   , canopy, trans   ,& ! JP add
                          z0rl  , ustar,                                                                   & ! JP add
                          pgr   , tgrs_1, qgrs_1,                                                          & ! JP add                           
                          z01d, zt1d, bexp1d, xlai1d, vegf1d, lndp_vgf,                                    &
                          cplflx, flag_cice, islmsk_cice, slimskin_cpl,                                    &
                          wind, u1, v1, cnvwind, smcwlt2, smcref2, vtype_save, stype_save, slope_save,     &
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
        integer, intent(in) :: nthreads, im, levs, isot, ivegsrc
        integer, dimension(:), intent(in) :: islmsk

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(:), intent(in) :: vfrac, prsik_1, prslk_1
        integer, dimension(:), intent(inout) :: vtype, stype, slope
        integer, dimension(:), intent(out)   :: vtype_save(:), stype_save(:), slope_save(:)

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

        ! Save current values of vegetation, soil and slope type
        vtype_save(:) = vtype(:)
        stype_save(:) = stype(:)
        slope_save(:) = slope(:)

        call update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype, slope)

        do i=1,im
          sigmaf(i) = max(vfrac(i), 0.01_kind_phys)
          islmsk_cice(i) = islmsk(i)

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
          soiltyp_cpl   (i) = stype(i)
          vegtype_cpl   (i) = vtype(i)
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
          slopetyp_cpl  (i) = slope(i)
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

      subroutine update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype, slope)

        implicit none

        integer, intent(in)    :: nthreads, im, isot, ivegsrc, islmsk(:)
        integer, intent(inout) :: vtype(:), stype(:), slope(:)
        integer :: i

!$OMP  parallel do num_threads(nthreads) default(none) private(i) &
!$OMP      shared(im, isot, ivegsrc, islmsk, vtype, stype, slope)
        do i=1,im
          if (islmsk(i) == 2) then
            if (isot == 1) then
              stype(i) = 16
            else
              stype(i) = 9
            endif
            if (ivegsrc == 0 .or. ivegsrc == 4) then
              vtype(i) = 24
            elseif (ivegsrc == 1) then
              vtype(i) = 15
            elseif (ivegsrc == 2) then
              vtype(i) = 13
            elseif (ivegsrc == 3 .or. ivegsrc == 5) then
              vtype(i) = 15
            endif
            slope(i)  = 9
          else
            if (vtype(i)  < 1) vtype(i)  = 17
            if (slope(i) < 1) slope(i) = 1
          endif
        enddo
!$OMP end parallel do

      end subroutine update_vegetation_soil_slope_type
!> @}

      end module GFS_surface_generic_pre
