MODULE diafwb
   !!======================================================================
   !!                       ***  MODULE  diafwb  ***
   !! Ocean diagnostics: freshwater budget
   !!======================================================================
   !! History :  8.2  !  01-02  (E. Durand)  Original code
   !!            8.5  !  02-06  (G. Madec)  F90: Free form and module
   !!            9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
   !!---------------------------------------------------------------------- 
   !!----------------------------------------------------------------------
   !!   Only for ORCA2 ORCA1 and ORCA025
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dia_fwb     : freshwater budget for global ocean configurations
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! ???
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC dia_fwb    ! routine called by step.F90

   REAL(wp)               ::   a_fwf ,          &
      &                        a_sshb, a_sshn, a_salb, a_saln
   REAL(wp), DIMENSION(4) ::   a_flxi, a_flxo, a_temi, a_temo, a_sali, a_salo

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diafwb.F90 5506 2015-06-29 15:19:38Z clevy $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_fwb( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_fwb  ***
      !!     
      !! ** Purpose :
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  :: inum             ! temporary logical unit
      INTEGER  :: ji, jj, jk, jt   ! dummy loop indices
      INTEGER  :: ii0, ii1, ij0, ij1
      INTEGER  :: isrow         ! index for ORCA1 starting row
      REAL(wp) :: zarea, zvol, zwei
      REAL(wp) :: ztemi(4), ztemo(4), zsali(4), zsalo(4), zflxi(4), zflxo(4)
      REAL(wp) :: zt, zs, zu  
      REAL(wp) :: zsm0, zfwfnew
      IF( ( cp_cfg == "orca" .OR. cp_cfg == "paleorca" ) .AND. jp_cfg == 1 .OR. jp_cfg == 2 .OR. jp_cfg == 4 ) THEN
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_fwb')

      ! Mean global salinity
      zsm0 = 34.72654

      ! To compute fwf mean value mean fwf

      IF( kt == nit000 ) THEN

         a_fwf    = 0.e0
         a_sshb   = 0.e0 ! valeur de ssh au debut de la simulation
         a_salb   = 0.e0 ! valeur de sal au debut de la simulation
         ! sshb used because diafwb called after tranxt (i.e. after the swap)
         a_sshb = SUM( e1t(:,:) * e2t(:,:) * sshb(:,:) * tmask_i(:,:) )
         IF( lk_mpp )   CALL mpp_sum( a_sshb )      ! sum over the global domain

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwei  = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
                  a_salb = a_salb + ( tsb(ji,jj,jk,jp_sal) - zsm0 ) * zwei
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum( a_salb )      ! sum over the global domain
      ENDIF
      
      a_fwf    = SUM( e1t(:,:) * e2t(:,:) * ( emp(:,:)-rnf(:,:) ) * tmask_i(:,:) ) 
      IF( lk_mpp )   CALL mpp_sum( a_fwf    )       ! sum over the global domain

      IF( kt == nitend ) THEN
         a_sshn = 0.e0
         a_saln = 0.e0
         zarea = 0.e0
         zvol  = 0.e0
         zfwfnew = 0.e0
         ! Mean sea level at nitend
         a_sshn = SUM( e1t(:,:) * e2t(:,:) * sshn(:,:) * tmask_i(:,:) )
         IF( lk_mpp )   CALL mpp_sum( a_sshn )      ! sum over the global domain
         zarea  = SUM( e1t(:,:) * e2t(:,:) *             tmask_i(:,:) )
         IF( lk_mpp )   CALL mpp_sum( zarea  )      ! sum over the global domain
         
         DO jk = 1, jpkm1   
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwei  = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
                  a_saln = a_saln + ( tsn(ji,jj,jk,jp_sal) - zsm0 ) * zwei
                  zvol  = zvol  + zwei
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum( a_saln )      ! sum over the global domain
         IF( lk_mpp )   CALL mpp_sum( zvol )      ! sum over the global domain
         
         ! Conversion in m3
         a_fwf    = a_fwf * rdttra(1) * 1.e-3 
         
         ! fwf correction to bring back the mean ssh to zero
         zfwfnew = a_sshn / ( ( nitend - nit000 + 1 ) * rdt ) * 1.e3 / zarea

      ENDIF


      ! Calcul des termes de transport 
      ! ------------------------------
      
      ! 1 --> Gibraltar 
      ! 2 --> Cadiz 
      ! 3 --> Red Sea
      ! 4 --> Baltic Sea

      IF( kt == nit000 ) THEN
         a_flxi(:) = 0.e0
         a_flxo(:) = 0.e0
         a_temi(:) = 0.e0
         a_temo(:) = 0.e0
         a_sali(:) = 0.e0
         a_salo(:) = 0.e0
      ENDIF

      zflxi(:) = 0.e0
      zflxo(:) = 0.e0
      ztemi(:) = 0.e0
      ztemo(:) = 0.e0
      zsali(:) = 0.e0
      zsalo(:) = 0.e0

      ! Mean flow at Gibraltar

      IF( cp_cfg == "orca" ) THEN  
               
         SELECT CASE ( jp_cfg )
         !                                           ! =======================
         CASE ( 4 )                                  !  ORCA_R4 configuration
            !                                        ! =======================
            ii0 = 70   ;   ii1 = 70
            ij0 = 52   ;   ij1 = 52
            !                                        ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration
            !                                        ! =======================
            ii0 = 140   ;   ii1 = 140
            ij0 = 102   ;   ij1 = 102
            !                                        ! =======================
         CASE ( 1 )                                  !  ORCA_R1 configurations
            !                                        ! =======================
            ! This dirty section will be suppressed by simplification process:
            ! all this will come back in input files
            ! Currently these hard-wired indices relate to configuration with
            ! extend grid (jpjglo=332)
            isrow = 332 - jpjglo
            !
            ii0 = 283           ;   ii1 = 283
            ij0 = 241 - isrow   ;   ij1 = 241 - isrow
            !                                        ! =======================
         CASE DEFAULT                                !    ORCA R05 or R025
            !                                        ! =======================
            CALL ctl_stop( ' dia_fwb Not yet implemented in ORCA_R05 or R025' )
            ! 
         END SELECT
         ! 
         DO ji = mi0(ii0), MIN(mi1(ii1),jpim1)
            DO jj = mj0(ij0), mj1(ij1)
               DO jk = 1, jpk 
                  zt = 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
                  zs = 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
                  zu = un(ji,jj,jk) * fse3t(ji,jj,jk) * e2u(ji,jj) * tmask_i(ji,jj)

                  IF( un(ji,jj,jk) > 0.e0 ) THEN 
                     zflxi(1) = zflxi(1) +    zu
                     ztemi(1) = ztemi(1) + zt*zu
                     zsali(1) = zsali(1) + zs*zu
                  ELSE
                     zflxo(1) = zflxo(1) +    zu
                     ztemo(1) = ztemo(1) + zt*zu
                     zsalo(1) = zsalo(1) + zs*zu
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      
      ! Mean flow at Cadiz
      IF( cp_cfg == "orca" ) THEN
               
         SELECT CASE ( jp_cfg )
         !                                           ! =======================
         CASE ( 4 )                                  !  ORCA_R4 configuration
            !                                        ! =======================
            ii0 = 69   ;   ii1 = 69
            ij0 = 52   ;   ij1 = 52
            !                                        ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration
            !                                        ! =======================
            ii0 = 137   ;   ii1 = 137
            ij0 = 101   ;   ij1 = 102
            !                                        ! =======================
         CASE ( 1 )                                  !  ORCA_R1 configurations
            !                                        ! =======================
            ! This dirty section will be suppressed by simplification process:
            ! all this will come back in input files
            ! Currently these hard-wired indices relate to configuration with
            ! extend grid (jpjglo=332)
            isrow = 332 - jpjglo
            ii0 = 282           ;   ii1 = 282
            ij0 = 240 - isrow   ;   ij1 = 240 - isrow
            !                                        ! =======================
         CASE DEFAULT                                !    ORCA R05 or R025
            !                                        ! =======================
            CALL ctl_stop( ' dia_fwb Not yet implemented in ORCA_R05 or R025' )
            ! 
         END SELECT
         ! 
         DO ji = mi0(ii0), MIN(mi1(ii1),jpim1)
            DO jj = mj0(ij0), mj1(ij1)
               DO jk = 1, jpk 
                  zt = 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
                  zs = 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
                  zu = un(ji,jj,jk) * fse3t(ji,jj,jk) * e2u(ji,jj) * tmask_i(ji,jj)
                  
                  IF( un(ji,jj,jk) > 0.e0 ) THEN 
                     zflxi(2) = zflxi(2) +    zu
                     ztemi(2) = ztemi(2) + zt*zu
                     zsali(2) = zsali(2) + zs*zu
                  ELSE
                     zflxo(2) = zflxo(2) +    zu
                     ztemo(2) = ztemo(2) + zt*zu
                     zsalo(2) = zsalo(2) + zs*zu
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      ! Mean flow at Red Sea entrance
      IF( cp_cfg == "orca" ) THEN
               
         SELECT CASE ( jp_cfg )
         !                                           ! =======================
         CASE ( 4 )                                  !  ORCA_R4 configuration
            !                                        ! =======================
            ii0 = 83   ;   ii1 = 83
            ij0 = 45   ;   ij1 = 45
            !                                        ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration
            !                                        ! =======================
            ii0 = 160   ;   ii1 = 160
            ij0 = 88    ;   ij1 = 88 
            !                                        ! =======================
         CASE ( 1 )                                  !  ORCA_R1 configurations
            !                                        ! =======================
            ! This dirty section will be suppressed by simplification process:
            ! all this will come back in input files
            ! Currently these hard-wired indices relate to configuration with
            ! extend grid (jpjglo=332)
            isrow = 332 - jpjglo
            ii0 = 331           ;   ii1 = 331
            ij0 = 215 - isrow   ;   ij1 = 215 - isrow
            !                                        ! =======================
         CASE DEFAULT                                !    ORCA R05 or R025
            !                                        ! =======================
            CALL ctl_stop( ' dia_fwb Not yet implemented in ORCA_R05 or R025' )
            ! 
         END SELECT
         ! 
         DO ji = mi0(ii0), MIN(mi1(ii1),jpim1)
            DO jj = mj0(ij0), mj1(ij1)
               DO jk = 1, jpk 
                  zt = 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
                  zs = 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
                  zu = un(ji,jj,jk) * fse3t(ji,jj,jk) * e2u(ji,jj) * tmask_i(ji,jj)
                  
                  IF( un(ji,jj,jk) > 0.e0 ) THEN 
                     zflxi(3) = zflxi(3) +    zu
                     ztemi(3) = ztemi(3) + zt*zu
                     zsali(3) = zsali(3) + zs*zu
                  ELSE
                     zflxo(3) = zflxo(3) +    zu
                     ztemo(3) = ztemo(3) + zt*zu
                     zsalo(3) = zsalo(3) + zs*zu
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      ! Mean flow at Baltic Sea entrance
      IF( cp_cfg == "orca" ) THEN
               
         SELECT CASE ( jp_cfg )
         !                                           ! =======================
         CASE ( 4 )                                  !  ORCA_R4 configuration
            !                                        ! =======================
            ii0 = 1     ;   ii1 = 1  
            ij0 = 1     ;   ij1 = 1  
            !                                        ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration
            !                                        ! =======================
            ii0 = 146   ;   ii1 = 146 
            ij0 = 116   ;   ij1 = 116
            !                                        ! =======================
         CASE ( 1 )                                  !  ORCA_R1 configurations
            !                                        ! =======================
            ! This dirty section will be suppressed by simplification process:
            ! all this will come back in input files
            ! Currently these hard-wired indices relate to configuration with
            ! extend grid (jpjglo=332)
            isrow = 332 - jpjglo
            ii0 = 297           ;   ii1 = 297
            ij0 = 269 - isrow   ;   ij1 = 269 - isrow
            !                                        ! =======================
         CASE DEFAULT                                !    ORCA R05 or R025
            !                                        ! =======================
            CALL ctl_stop( ' dia_fwb Not yet implemented in ORCA_R05 or R025' )
            ! 
         END SELECT
         ! 
         DO ji = mi0(ii0), MIN(mi1(ii1),jpim1)
            DO jj = mj0(ij0), mj1(ij1)
               DO jk = 1, jpk
                  zt = 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
                  zs = 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
                  zu = un(ji,jj,jk) * fse3t(ji,jj,jk) * e2u(ji,jj) * tmask_i(ji,jj)
                  
                  IF( un(ji,jj,jk) > 0.e0 ) THEN 
                     zflxi(4) = zflxi(4) +    zu
                     ztemi(4) = ztemi(4) + zt*zu
                     zsali(4) = zsali(4) + zs*zu
                  ELSE
                     zflxo(4) = zflxo(4) +    zu
                     ztemo(4) = ztemo(4) + zt*zu
                     zsalo(4) = zsalo(4) + zs*zu
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      ! Sum at each time-step
      DO jt = 1, 4 
         !
         IF( zflxi(jt) /= 0.e0 ) THEN
            a_flxi(jt) = a_flxi(jt) + zflxi(jt)
            a_temi(jt) = a_temi(jt) + ztemi(jt)/zflxi(jt)
            a_sali(jt) = a_sali(jt) + zsali(jt)/zflxi(jt)
         ENDIF
         !
         IF( zflxo(jt) /= 0.e0 ) THEN
            a_flxo(jt) = a_flxo(jt) + zflxo(jt)
            a_temo(jt) = a_temo(jt) + ztemo(jt)/zflxo(jt)
            a_salo(jt) = a_salo(jt) + zsalo(jt)/zflxo(jt)
         ENDIF
         !
      END DO

      IF( kt == nitend ) THEN
         DO jt = 1, 4 
            a_flxi(jt) = a_flxi(jt) / ( FLOAT( nitend - nit000 + 1 ) * 1.e6 )
            a_temi(jt) = a_temi(jt) /   FLOAT( nitend - nit000 + 1 )
            a_sali(jt) = a_sali(jt) /   FLOAT( nitend - nit000 + 1 )
            a_flxo(jt) = a_flxo(jt) / ( FLOAT( nitend - nit000 + 1 ) * 1.e6 )
            a_temo(jt) = a_temo(jt) /   FLOAT( nitend - nit000 + 1 )
            a_salo(jt) = a_salo(jt) /   FLOAT( nitend - nit000 + 1 )
         END DO
         IF( lk_mpp ) THEN
            CALL mpp_sum( a_flxi, 4 )      ! sum over the global domain
            CALL mpp_sum( a_temi, 4 )      ! sum over the global domain
            CALL mpp_sum( a_sali, 4 )      ! sum over the global domain

            CALL mpp_sum( a_flxo, 4 )      ! sum over the global domain
            CALL mpp_sum( a_temo, 4 )      ! sum over the global domain
            CALL mpp_sum( a_salo, 4 )      ! sum over the global domain
         ENDIF
      ENDIF


      ! Ecriture des diagnostiques 
      ! --------------------------

      IF ( kt == nitend .AND. ( cp_cfg == "orca" .OR.  cp_cfg == "paleorca" ) .AND. lwp ) THEN

         CALL ctl_opn( inum, 'STRAIT.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         WRITE(inum,*)
         WRITE(inum,*)    'Net freshwater budget '
         WRITE(inum,9010) '  fwf    = ',a_fwf,   ' m3 =', a_fwf   /(FLOAT(nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
         WRITE(inum,*)
         WRITE(inum,9010) '  zarea =',zarea
         WRITE(inum,9010) '  zvol  =',zvol
         WRITE(inum,*)
         WRITE(inum,*)    'Mean sea level : '
         WRITE(inum,9010) '  at nit000 = ',a_sshb        ,' m3 '
         WRITE(inum,9010) '  at nitend = ',a_sshn        ,' m3 '
         WRITE(inum,9010) '  diff      = ',(a_sshn-a_sshb),' m3 =', (a_sshn-a_sshb)/(FLOAT(nitend-nit000+1)*rdt) * 1.e-6,' Sv'
         WRITE(inum,9020) '  mean sea level elevation    =', a_sshn/zarea,' m'
         WRITE(inum,*)
         WRITE(inum,*)    'Anomaly of salinity content : '
         WRITE(inum,9010) '  at nit000 = ',a_salb        ,' psu.m3 '
         WRITE(inum,9010) '  at nitend = ',a_saln        ,' psu.m3 '
         WRITE(inum,9010) '  diff      = ',(a_saln-a_salb),' psu.m3'
         WRITE(inum,*)
         WRITE(inum,*)    'Mean salinity : '
         WRITE(inum,9020) '  at nit000 =',a_salb/zvol+zsm0   ,' psu '
         WRITE(inum,9020) '  at nitend =',a_saln/zvol+zsm0   ,' psu '
         WRITE(inum,9020) '  diff      =',(a_saln-a_salb)/zvol,' psu'
         WRITE(inum,9020) '  S-SLevitus=',a_saln/zvol,' psu'
         WRITE(inum,*)
         WRITE(inum,*)    'Gibraltar : '
         WRITE(inum,9030) '  Flux entrant (Sv) :', a_flxi(1)
         WRITE(inum,9030) '  Flux sortant (Sv) :', a_flxo(1)
         WRITE(inum,9030) '  T entrant (deg)   :', a_temi(1)
         WRITE(inum,9030) '  T sortant (deg)   :', a_temo(1)
         WRITE(inum,9030) '  S entrant (psu)   :', a_sali(1)
         WRITE(inum,9030) '  S sortant (psu)   :', a_salo(1)
         WRITE(inum,*)
         WRITE(inum,*)    'Cadiz : '
         WRITE(inum,9030) '  Flux entrant (Sv) :', a_flxi(2)
         WRITE(inum,9030) '  Flux sortant (Sv) :', a_flxo(2)
         WRITE(inum,9030) '  T entrant (deg)   :', a_temi(2)
         WRITE(inum,9030) '  T sortant (deg)   :', a_temo(2)
         WRITE(inum,9030) '  S entrant (psu)   :', a_sali(2)
         WRITE(inum,9030) '  S sortant (psu)   :', a_salo(2)
         WRITE(inum,*)
         WRITE(inum,*)    'Bab el Mandeb : '
         WRITE(inum,9030) '  Flux entrant (Sv) :', a_flxi(3)
         WRITE(inum,9030) '  Flux sortant (Sv) :', a_flxo(3)
         WRITE(inum,9030) '  T entrant (deg)   :', a_temi(3)
         WRITE(inum,9030) '  T sortant (deg)   :', a_temo(3)
         WRITE(inum,9030) '  S entrant (psu)   :', a_sali(3)
         WRITE(inum,9030) '  S sortant (psu)   :', a_salo(3)
         WRITE(inum,*)
         WRITE(inum,*)    'Baltic : '
         WRITE(inum,9030) '  Flux entrant (Sv) :', a_flxi(4)
         WRITE(inum,9030) '  Flux sortant (Sv) :', a_flxo(4)
         WRITE(inum,9030) '  T entrant (deg)   :', a_temi(4)
         WRITE(inum,9030) '  T sortant (deg)   :', a_temo(4)
         WRITE(inum,9030) '  S entrant (psu)   :', a_sali(4)
         WRITE(inum,9030) '  S sortant (psu)   :', a_salo(4)
         CLOSE(inum)
      ENDIF

      IF( nn_timing == 1 )   CALL timing_start('dia_fwb')

 9005 FORMAT(1X,A,ES24.16)
 9010 FORMAT(1X,A,ES12.5,A,F10.5,A)
 9020 FORMAT(1X,A,F10.5,A)
 9030 FORMAT(1X,A,F9.4,A)
  
      ENDIF 

   END SUBROUTINE dia_fwb

   !!======================================================================
END MODULE diafwb
