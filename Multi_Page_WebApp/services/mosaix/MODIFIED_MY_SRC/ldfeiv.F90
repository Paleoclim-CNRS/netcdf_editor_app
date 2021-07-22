MODULE ldfeiv
   !!======================================================================
   !!                     ***  MODULE  ldfeiv  ***
   !! Ocean physics:  variable eddy induced velocity coefficients
   !!======================================================================
   !! History :  OPA  ! 1999-03  (G. Madec, A. Jouzeau)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  Free form, F90
   !!----------------------------------------------------------------------
#if   defined key_traldf_eiv   &&   defined key_traldf_c2d
   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'      and                     eddy induced velocity
   !!   'key_traldf_c2d'                    2D tracer lateral  mixing coef.
   !!----------------------------------------------------------------------
   !!   ldf_eiv      : compute the eddy induced velocity coefficients
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcrnf          ! river runoffs
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE phycst          ! physical constants
   USE ldfslp          ! iso-neutral slopes
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! I/O library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   ldf_eiv    ! routine called by step.F90
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfeiv.F90 6631 2016-05-27 07:12:30Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_eiv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!              growth rate of baroclinic instability.
      !!
      !! ** Method :
      !!
      !! ** Action : - uslp , vslp  : i- and j-slopes of neutral surfaces at u- & v-points
      !!             - wslpi, wslpj : i- and j-slopes of neutral surfaces at w-points. 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step inedx
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zfw, ze3w, zn2, zf20, zaht, zaht_min      ! temporary scalars
      REAL(wp), DIMENSION(:,:), POINTER ::   zn, zah, zhw, zross   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_eiv')
      !
      CALL wrk_alloc( jpi,jpj, zn, zah, zhw, zross )

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF
      
      ! 0. Local initialization
      ! -----------------------
      zn   (:,:) = 0._wp
      zhw  (:,:) = 5._wp
      zah  (:,:) = 0._wp
      zross(:,:) = 0._wp


      ! 1. Compute lateral diffusive coefficient 
      ! ----------------------------------------
      IF( ln_traldf_grif ) THEN
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ! Take the max of N^2 and zero then take the vertical sum 
                  ! of the square root of the resulting N^2 ( required to compute 
                  ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
                  zn2 = MAX( rn2b(ji,jj,jk), 0._wp )
                  zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * fse3w(ji,jj,jk)
                  ! Compute elements required for the inverse time scale of baroclinic
                  ! eddies using the isopycnal slopes calculated in ldfslp.F : 
                  ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
                  ze3w = fse3w(ji,jj,jk) * tmask(ji,jj,jk)
                  zah(ji,jj) = zah(ji,jj) + zn2 * wslp2(ji,jj,jk) * ze3w
                  zhw(ji,jj) = zhw(ji,jj) + ze3w
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ! Take the max of N^2 and zero then take the vertical sum 
                  ! of the square root of the resulting N^2 ( required to compute 
                  ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
                  zn2 = MAX( rn2b(ji,jj,jk), 0._wp )
                  zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * fse3w(ji,jj,jk)
                  ! Compute elements required for the inverse time scale of baroclinic
                  ! eddies using the isopycnal slopes calculated in ldfslp.F : 
                  ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
                  ze3w = fse3w(ji,jj,jk) * tmask(ji,jj,jk)
                  zah(ji,jj) = zah(ji,jj) + zn2 * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)   &
                     &                            + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) ) * ze3w
                  zhw(ji,jj) = zhw(ji,jj) + ze3w
               END DO
            END DO
         END DO
      END IF

      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
            ! Rossby radius at w-point taken < 40km and  > 2km
            zross(ji,jj) = MAX( MIN( .4 * zn(ji,jj) / zfw, 40.e3 ), 2.e3 )
            ! Compute aeiw by multiplying Ro^2 and T^-1
            aeiw(ji,jj) = zross(ji,jj) * zross(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * ssmask(ji,jj)
         END DO
      END DO

      ! Add PALEORCA2 configuration -- JBL 08.02.2017
      IF( ( cp_cfg == "orca" .OR. cp_cfg == "paleorca" ) .AND. jp_cfg == 2 ) THEN   ! ORCA R2
         DO jj = 2, jpjm1
!CDIR NOVERRCHK 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! Take the minimum between aeiw and 1000 m2/s over shelves (depth shallower than 650 m)
               IF( mbkt(ji,jj) <= 20 )   aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
            END DO
         END DO
      ENDIF

      ! Decrease the coefficient in the tropics (20N-20S) 
      zf20 = 2._wp * omega * sin( rad * 20._wp )
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiw(ji,jj) = MIN( 1., ABS( ff(ji,jj) / zf20 ) ) * aeiw(ji,jj)
         END DO
      END DO

      ! ORCA R05: Take the minimum between aeiw  and aeiv0
      IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               aeiw(ji,jj) = MIN( aeiw(ji,jj), aeiv0 )
            END DO
         END DO
      ENDIF

      ! ORCA R1: Take the minimum between aeiw  and aeiv0
      IF( cp_cfg == "orca" .AND. jp_cfg == 1 ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               aeiw(ji,jj) = MIN( aeiw(ji,jj), aeiv0 )
            END DO
         END DO
      ENDIF

      CALL lbc_lnk( aeiw, 'W', 1. )      ! lateral boundary condition on aeiw 


      ! Average the diffusive coefficient at u- v- points 
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiu(ji,jj) = 0.5_wp * ( aeiw(ji,jj) + aeiw(ji+1,jj  ) )
            aeiv(ji,jj) = 0.5_wp * ( aeiw(ji,jj) + aeiw(ji  ,jj+1) )
         END DO 
      END DO 
      CALL lbc_lnk( aeiu, 'U', 1. )   ;   CALL lbc_lnk( aeiv, 'V', 1. )      ! lateral boundary condition


      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=aeiu, clinfo1=' eiv  - u: ', ovlap=1)
         CALL prt_ctl(tab2d_1=aeiv, clinfo1=' eiv  - v: ', ovlap=1)
      ENDIF

      ! ORCA R05: add a space variation on aht (=aeiv except at the equator and river mouth)
      IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
         zf20     = 2._wp * omega * SIN( rad * 20._wp )
         zaht_min = 100._wp                           ! minimum value for aht
         DO jj = 1, jpj
            DO ji = 1, jpi
               zaht      = ( 1._wp -  MIN( 1._wp , ABS( ff(ji,jj) / zf20 ) ) ) * ( aht0 - zaht_min )  &
                  &      + aht0 * rnfmsk(ji,jj)                          ! enhanced near river mouths
               ahtu(ji,jj) = MAX( MAX( zaht_min, aeiu(ji,jj) ) + zaht, aht0 )
               ahtv(ji,jj) = MAX( MAX( zaht_min, aeiv(ji,jj) ) + zaht, aht0 )
               ahtw(ji,jj) = MAX( MAX( zaht_min, aeiw(ji,jj) ) + zaht, aht0 )
            END DO
         END DO
         IF(ln_ctl) THEN
            CALL prt_ctl(tab2d_1=ahtu, clinfo1=' aht  - u: ', ovlap=1)
            CALL prt_ctl(tab2d_1=ahtv, clinfo1=' aht  - v: ', ovlap=1)
            CALL prt_ctl(tab2d_1=ahtw, clinfo1=' aht  - w: ', ovlap=1)
         ENDIF
      ENDIF

      IF( aeiv0 == 0._wp ) THEN
         aeiu(:,:) = 0._wp
         aeiv(:,:) = 0._wp
         aeiw(:,:) = 0._wp
      ENDIF

      CALL iom_put( "aht2d"    , ahtw )   ! lateral eddy diffusivity
      CALL iom_put( "aht2d_eiv", aeiw )   ! EIV lateral eddy diffusivity
      !  
      CALL wrk_dealloc( jpi,jpj, zn, zah, zhw, zross )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_eiv')
      !
   END SUBROUTINE ldf_eiv

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ldf_eiv( kt )       ! Empty routine
      INTEGER :: kt
      WRITE(*,*) 'ldf_eiv: You should not have seen this print! error?', kt
   END SUBROUTINE ldf_eiv
#endif

   !!======================================================================
END MODULE ldfeiv
