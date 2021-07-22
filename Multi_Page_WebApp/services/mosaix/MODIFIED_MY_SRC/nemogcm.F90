MODULE nemogcm
   !!======================================================================
   !!                       ***  MODULE nemogcm   ***
   !! Ocean system   : NEMO GCM (ocean dynamics, on-line tracers, biochemistry and sea-ice)
   !!======================================================================
   !! History :  OPA  ! 1990-10  (C. Levy, G. Madec)  Original code
   !!            7.0  ! 1991-11  (M. Imbard, C. Levy, G. Madec)
   !!            7.1  ! 1993-03  (M. Imbard, C. Levy, G. Madec, O. Marti, M. Guyon, A. Lazar,
   !!                             P. Delecluse, C. Perigaud, G. Caniaux, B. Colot, C. Maes) release 7.1
   !!             -   ! 1992-06  (L.Terray)  coupling implementation
   !!             -   ! 1993-11  (M.A. Filiberti) IGLOO sea-ice
   !!            8.0  ! 1996-03  (M. Imbard, C. Levy, G. Madec, O. Marti, M. Guyon, A. Lazar,
   !!                             P. Delecluse, L.Terray, M.A. Filiberti, J. Vialar, A.M. Treguier, M. Levy) release 8.0
   !!            8.1  ! 1997-06  (M. Imbard, G. Madec)
   !!            8.2  ! 1999-11  (M. Imbard, H. Goosse)  LIM sea-ice model
   !!                 ! 1999-12  (V. Thierry, A-M. Treguier, M. Imbard, M-A. Foujols)  OPEN-MP
   !!                 ! 2000-07  (J-M Molines, M. Imbard)  Open Boundary Conditions  (CLIPPER)
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and modules
   !!             -   ! 2004-06  (R. Redler, NEC CCRLE, Germany) add OASIS[3/4] coupled interfaces
   !!             -   ! 2004-08  (C. Talandier) New trends organization
   !!             -   ! 2005-06  (C. Ethe) Add the 1D configuration possibility
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!             -   ! 2006-03  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   ! 2006-04  (G. Madec, R. Benshila)  Step reorganization
   !!             -   ! 2007-07  (J. Chanut, A. Sellar) Unstructured open boundaries (BDY)
   !!            3.2  ! 2009-08  (S. Masson)  open/write in the listing file in mpp
   !!            3.3  ! 2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.4  ! 2011-11  (C. Harris) decomposition changes for running with CICE
   !!                 ! 2012-05  (C. Calone, J. Simeon, G. Madec, C. Ethe) Add grid coarsening 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   nemo_gcm       : solve ocean dynamics, tracer, biogeochemistry and/or sea-ice
   !!   nemo_init      : initialization of the NEMO system
   !!   nemo_ctl       : initialisation of the contol print
   !!   nemo_closefile : close remaining open files
   !!   nemo_alloc     : dynamical allocation
   !!   nemo_partition : calculate MPP domain decomposition
   !!   factorise      : calculate the factors of the no. of MPI processes
   !!----------------------------------------------------------------------
   USE step_oce        ! module used in the ocean time stepping module
   USE cla             ! cross land advection               (tra_cla routine)
   USE domcfg          ! domain configuration               (dom_cfg routine)
   USE mppini          ! shared/distributed memory setting (mpp_init routine)
   USE domain          ! domain initialization             (dom_init routine)
#if defined key_nemocice_decomp
   USE ice_domain_size, only: nx_global, ny_global
#endif
   USE tideini         ! tidal components initialization   (tide_ini routine)
   USE bdyini          ! open boundary cond. setting       (bdy_init routine)
   USE bdydta          ! open boundary cond. setting   (bdy_dta_init routine)
   USE bdytides        ! open boundary cond. setting   (bdytide_init routine)
   USE istate          ! initial state setting          (istate_init routine)
   USE ldfdyn          ! lateral viscosity setting      (ldfdyn_init routine)
   USE ldftra          ! lateral diffusivity setting    (ldftra_init routine)
   USE zdfini          ! vertical physics setting          (zdf_init routine)
   USE phycst          ! physical constant                  (par_cst routine)
   USE trdini          ! dyn/tra trends initialization     (trd_init routine)
   USE asminc          ! assimilation increments     
   USE asmbkg          ! writing out state trajectory
   USE diaptr          ! poleward transports           (dia_ptr_init routine)
   USE diadct          ! sections transports           (dia_dct_init routine)
   USE diaobs          ! Observation diagnostics       (dia_obs_init routine)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE step            ! NEMO time-stepping                 (stp     routine)
   USE icbini          ! handle bergs, initialisation
   USE icbstp          ! handle bergs, calving, themodynamics and transport
   USE cpl_oasis3      ! OASIS3 coupling
   USE c1d             ! 1D configuration
   USE step_c1d        ! Time stepping loop for the 1D configuration
   USE dyndmp          ! Momentum damping
#if defined key_top
   USE trcini          ! passive tracer initialisation
#endif
   USE lib_mpp         ! distributed memory computing
#if defined key_iomput
   USE xios
#endif
   USE sbctide, ONLY: lk_tide
   USE crsini          ! initialise grid coarsening utility
   USE lbcnfd, ONLY: isendto, nsndto, nfsloop, nfeloop ! Setup of north fold exchanges 
   USE sbc_oce, ONLY: lk_oasis
   USE stopar
   USE stopts

   IMPLICIT NONE
   PRIVATE

   PUBLIC   nemo_gcm    ! called by model.F90
   PUBLIC   nemo_init   ! needed by AGRIF
   PUBLIC   nemo_alloc  ! needed by TAM

   CHARACTER(lc) ::   cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: nemogcm.F90 8566 2017-09-27 13:15:25Z mchekki $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_gcm
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_gcm  ***
      !!
      !! ** Purpose :   NEMO solves the primitive equations on an orthogonal
      !!              curvilinear mesh on the sphere.
      !!
      !! ** Method  : - model general initialization
      !!              - launch the time-stepping (stp routine)
      !!              - finalize the run by closing files and communications
      !!
      !! References : Madec, Delecluse, Imbard, and Levy, 1997:  internal report, IPSL.
      !!              Madec, 2008, internal report, IPSL.
      !!----------------------------------------------------------------------
      INTEGER ::   istp       ! time step index
      !!----------------------------------------------------------------------
      !
#if defined key_agrif
      CALL Agrif_Init_Grids()      ! AGRIF: set the meshes
#endif

      !                            !-----------------------!
      CALL nemo_init               !==  Initialisations  ==!
      !                            !-----------------------!
#if defined key_agrif
      CALL Agrif_Declare_Var_dom   ! AGRIF: set the meshes for DOM
      CALL Agrif_Declare_Var       !  "      "   "   "      "  DYN/TRA 
# if defined key_top
      CALL Agrif_Declare_Var_top   !  "      "   "   "      "  TOP
# endif
# if defined key_lim2
      CALL Agrif_Declare_Var_lim2  !  "      "   "   "      "  LIM
# endif
#endif
      ! check that all process are still there... If some process have an error,
      ! they will never enter in step and other processes will wait until the end of the cpu time!
      IF( lk_mpp )   CALL mpp_max( nstop )

      IF(lwp) WRITE(numout,cform_aaa)   ! Flag AAAAAAA

      !                            !-----------------------!
      !                            !==   time stepping   ==!
      !                            !-----------------------!
      istp = nit000
#if defined key_c1d
         DO WHILE ( istp <= nitend .AND. nstop == 0 )
            CALL stp_c1d( istp )
            istp = istp + 1
         END DO
#else
          IF( lk_asminc ) THEN
             IF( ln_bkgwri ) CALL asm_bkg_wri( nit000 - 1 )    ! Output background fields
             IF( ln_asmdin ) THEN                        ! Direct initialization
                IF( ln_trainc ) CALL tra_asm_inc( nit000 - 1 )    ! Tracers
                IF( ln_dyninc ) CALL dyn_asm_inc( nit000 - 1 )    ! Dynamics
                IF( ln_sshinc ) CALL ssh_asm_inc( nit000 - 1 )    ! SSH
             ENDIF
          ENDIF

#if defined key_agrif
          CALL Agrif_Regrid()
#endif

         DO WHILE ( istp <= nit000 .AND. nstop == 0 )
#if defined key_agrif
            CALL stp                         ! AGRIF: time stepping
#else
            CALL stp( istp )                 ! standard time stepping
#endif
            istp = istp + 1
            IF( lk_mpp )   CALL mpp_max( nstop )
         END DO
#endif

      IF( lk_diaobs   )   CALL dia_obs_wri
      !
      IF( ln_icebergs )   CALL icb_end( nitend )

      !                            !------------------------!
      !                            !==  finalize the run  ==!
      !                            !------------------------!
      IF(lwp) WRITE(numout,cform_aaa)   ! Flag AAAAAAA
      !
      IF( nstop /= 0 .AND. lwp ) THEN   ! error print
         WRITE(numout,cform_err)
         WRITE(numout,*) nstop, ' error have been found'
      ENDIF
      !
#if defined key_agrif
      CALL Agrif_ParentGrid_To_ChildGrid()
      IF( lk_diaobs ) CALL dia_obs_wri
      IF( nn_timing == 1 )   CALL timing_finalize
      CALL Agrif_ChildGrid_To_ParentGrid()
#endif
      IF( nn_timing == 1 )   CALL timing_finalize
      !
      CALL nemo_closefile
      !
#if defined key_iomput
      CALL xios_finalize                ! end mpp communications with xios
      IF( lk_oasis ) CALL cpl_finalize    ! end coupling and mpp communications with OASIS
#else
      IF( lk_oasis ) THEN 
         CALL cpl_finalize              ! end coupling and mpp communications with OASIS
      ELSE
         IF( lk_mpp )   CALL mppstop    ! end mpp communications
      ENDIF
#endif
      !
   END SUBROUTINE nemo_gcm


   SUBROUTINE nemo_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_init  ***
      !!
      !! ** Purpose :   initialization of the NEMO GCM
      !!----------------------------------------------------------------------
      INTEGER ::   ji            ! dummy loop indices
      INTEGER ::   ilocal_comm   ! local integer
      INTEGER ::   ios
      CHARACTER(len=80), DIMENSION(16) ::   cltxt
      !
      NAMELIST/namctl/ ln_ctl  , nn_print, nn_ictls, nn_ictle,   &
         &             nn_isplt, nn_jsplt, nn_jctls, nn_jctle,   &
         &             nn_bench, nn_timing
      NAMELIST/namcfg/ cp_cfg, cp_cfz, jp_cfg, jpidta, jpjdta, jpkdta, jpiglo, jpjglo, &
         &             jpizoom, jpjzoom, jperio, ln_use_jattr
      !!----------------------------------------------------------------------
      !
      cltxt = ''
      cxios_context = 'nemo'
      !
      !                             ! Open reference namelist and configuration namelist files
      CALL ctl_opn( numnam_ref, 'namelist_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
      CALL ctl_opn( numnam_cfg, 'namelist_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
      !
      REWIND( numnam_ref )              ! Namelist namctl in reference namelist : Control prints & Benchmark
      READ  ( numnam_ref, namctl, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namctl in reference namelist', .TRUE. )

      REWIND( numnam_cfg )              ! Namelist namctl in confguration namelist : Control prints & Benchmark
      READ  ( numnam_cfg, namctl, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namctl in configuration namelist', .TRUE. )

      !
      REWIND( numnam_ref )              ! Namelist namcfg in reference namelist : Control prints & Benchmark
      READ  ( numnam_ref, namcfg, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcfg in reference namelist', .TRUE. )

      REWIND( numnam_cfg )              ! Namelist namcfg in confguration namelist : Control prints & Benchmark
      READ  ( numnam_cfg, namcfg, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcfg in configuration namelist', .TRUE. )   

! Force values for AGRIF zoom (cf. agrif_user.F90)
#if defined key_agrif
   IF( .NOT. Agrif_Root() ) THEN
      jpiglo  = nbcellsx + 2 + 2*nbghostcells
      jpjglo  = nbcellsy + 2 + 2*nbghostcells
      jpi     = ( jpiglo-2*jpreci + (jpni-1+0) ) / jpni + 2*jpreci
      jpj     = ( jpjglo-2*jprecj + (jpnj-1+0) ) / jpnj + 2*jprecj
      jpidta  = jpiglo
      jpjdta  = jpjglo
      jpizoom = 1
      jpjzoom = 1
      nperio  = 0
      jperio  = 0
      ln_use_jattr = .false.
   ENDIF
#endif
      !
      !                             !--------------------------------------------!
      !                             !  set communicator & select the local node  !
      !                             !  NB: mynode also opens output.namelist.dyn !
      !                             !      on unit number numond on first proc   !
      !                             !--------------------------------------------!
#if defined key_iomput
      IF( Agrif_Root() ) THEN
         IF( lk_oasis ) THEN
            CALL cpl_init( "oceanx", ilocal_comm )                     ! nemo local communicator given by oasis
            CALL xios_initialize( "not used",local_comm=ilocal_comm )    ! send nemo communicator to xios
         ELSE
            CALL  xios_initialize( "for_xios_mpi_id",return_comm=ilocal_comm )    ! nemo local communicator given by xios
         ENDIF
      ENDIF
      ! Nodes selection (control print return in cltxt)
      narea = mynode( cltxt, 'output.namelist.dyn', numnam_ref, numnam_cfg, numond , nstop, ilocal_comm )
#else
      IF( lk_oasis ) THEN
         IF( Agrif_Root() ) THEN
            CALL cpl_init( "oceanx", ilocal_comm )                      ! nemo local communicator given by oasis
         ENDIF
         ! Nodes selection (control print return in cltxt)
         narea = mynode( cltxt, 'output.namelist.dyn', numnam_ref, numnam_cfg, numond , nstop, ilocal_comm )
      ELSE
         ilocal_comm = 0
         ! Nodes selection (control print return in cltxt)
         narea = mynode( cltxt, 'output.namelist.dyn', numnam_ref, numnam_cfg, numond , nstop )
      ENDIF
#endif
      narea = narea + 1                                     ! mynode return the rank of proc (0 --> jpnij -1 )

      lwm = (narea == 1)                                    ! control of output namelists
      lwp = (narea == 1) .OR. ln_ctl                        ! control of all listing output print

      IF(lwm) THEN
         ! write merged namelists from earlier to output namelist now that the
         ! file has been opened in call to mynode. nammpp has already been
         ! written in mynode (if lk_mpp_mpi)
         WRITE( numond, namctl )
         WRITE( numond, namcfg )
      ENDIF

      ! If dimensions of processor grid weren't specified in the namelist file
      ! then we calculate them here now that we have our communicator size
      IF( (jpni < 1) .OR. (jpnj < 1) )THEN
#if   defined key_mpp_mpi
         IF( Agrif_Root() ) CALL nemo_partition(mppsize)
#else
         jpni  = 1
         jpnj  = 1
         jpnij = jpni*jpnj
#endif
      END IF

      ! Calculate domain dimensions given calculated jpni and jpnj
      ! This used to be done in par_oce.F90 when they were parameters rather
      ! than variables
      IF( Agrif_Root() ) THEN
#if defined key_nemocice_decomp
         jpi = ( nx_global+2-2*jpreci + (jpni-1) ) / jpni + 2*jpreci ! first  dim.
         jpj = ( ny_global+2-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj ! second dim. 
#else
         jpi = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci   ! first  dim.
         jpj = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj   ! second dim.
#endif
      ENDIF         
         jpk = jpkdta                                             ! third dim
#if defined key_agrif
         ! simple trick to use same vertical grid as parent
         ! but different number of levels: 
         ! Save maximum number of levels in jpkdta, then define all vertical grids
         ! with this number.
         ! Suppress once vertical online interpolation is ok
         IF(.NOT.Agrif_Root()) jpkdta = Agrif_Parent(jpkdta)
#endif
         jpim1 = jpi-1                                            ! inner domain indices
         jpjm1 = jpj-1                                            !   "           "
         jpkm1 = jpk-1                                            !   "           "
         jpij  = jpi*jpj                                          !  jpi x j

      IF(lwp) THEN                            ! open listing units
         !
         CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
         !
         WRITE(numout,*)
         WRITE(numout,*) '   CNRS - NERC - Met OFFICE - MERCATOR-ocean - INGV - CMCC'
         WRITE(numout,*) '                       NEMO team'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '                  version 3.6  (2015) '
         WRITE(numout,*)
         WRITE(numout,*)
         DO ji = 1, SIZE(cltxt)
            IF( TRIM(cltxt(ji)) /= '' )   WRITE(numout,*) cltxt(ji)      ! control print of mynode
         END DO
         WRITE(numout,cform_aaa)                                         ! Flag AAAAAAA
         !
      ENDIF

      ! Now we know the dimensions of the grid and numout has been set we can
      ! allocate arrays
      CALL nemo_alloc()

      !                             !-------------------------------!
      !                             !  NEMO general initialization  !
      !                             !-------------------------------!

      CALL nemo_ctl                          ! Control prints & Benchmark

      !                                      ! Domain decomposition
      IF( jpni*jpnj == jpnij ) THEN   ;   CALL mpp_init      ! standard cutting out
      ELSE                            ;   CALL mpp_init2     ! eliminate land processors
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_init
      !
      !                                      ! General initialization
                            CALL     phy_cst    ! Physical constants
                            CALL     eos_init   ! Equation of state
      IF( lk_c1d        )   CALL     c1d_init   ! 1D column configuration
                            CALL     dom_cfg    ! Domain configuration
                            CALL     dom_init   ! Domain

      IF( ln_nnogather )    CALL nemo_northcomms   ! Initialise the northfold neighbour lists (must be done after the masks are defined)

      IF( ln_ctl        )   CALL prt_ctl_init   ! Print control

                            CALL  istate_init   ! ocean initial state (Dynamics and tracers)

      IF( lk_tide       )   CALL    tide_init( nit000 )    ! Initialisation of the tidal harmonics

                            CALL     sbc_init   ! Forcings : surface module (clem: moved here for bdy purpose)

      IF( lk_bdy        )   CALL     bdy_init   ! Open boundaries initialisation
      IF( lk_bdy        )   CALL bdy_dta_init   ! Open boundaries initialisation of external data arrays
      IF( lk_bdy .AND. lk_tide )   &
         &                  CALL bdytide_init   ! Open boundaries initialisation of tidal harmonic forcing

                            CALL dyn_nept_init  ! simplified form of Neptune effect
      !     
      IF( ln_crs        )   CALL     crs_init   ! Domain initialization of coarsened grid
      !
                                ! Ocean physics
      !                                         ! Vertical physics
                            CALL     zdf_init      ! namelist read
                            CALL zdf_bfr_init      ! bottom friction
      IF( lk_zdfric     )   CALL zdf_ric_init      ! Richardson number dependent Kz
      IF( lk_zdftke     )   CALL zdf_tke_init      ! TKE closure scheme
      IF( lk_zdfgls     )   CALL zdf_gls_init      ! GLS closure scheme
      IF( lk_zdfkpp     )   CALL zdf_kpp_init      ! KPP closure scheme
      IF( lk_zdftmx     )   CALL zdf_tmx_init      ! tidal vertical mixing
      IF( lk_zdfddm .AND. .NOT. lk_zdfkpp )   &
         &                  CALL zdf_ddm_init      ! double diffusive mixing
      !                                         ! Lateral physics
                            CALL ldf_tra_init      ! Lateral ocean tracer physics
                            CALL ldf_dyn_init      ! Lateral ocean momentum physics
      IF( lk_ldfslp     )   CALL ldf_slp_init      ! slope of lateral mixing

      !                                     ! Active tracers
                            CALL tra_qsr_init   ! penetrative solar radiation qsr
                            CALL tra_bbc_init   ! bottom heat flux
      IF( lk_trabbl     )   CALL tra_bbl_init   ! advective (and/or diffusive) bottom boundary layer scheme
                            CALL tra_dmp_init   ! internal damping trends- tracers
                            CALL tra_adv_init   ! horizontal & vertical advection
                            CALL tra_ldf_init   ! lateral mixing
                            CALL tra_zdf_init   ! vertical mixing and after tracer fields

      !                                     ! Dynamics
      IF( lk_c1d        )   CALL dyn_dmp_init   ! internal damping trends- momentum
                            CALL dyn_adv_init   ! advection (vector or flux form)
                            CALL dyn_vor_init   ! vorticity term including Coriolis
                            CALL dyn_ldf_init   ! lateral mixing
                            CALL dyn_hpg_init   ! horizontal gradient of Hydrostatic pressure
                            CALL dyn_zdf_init   ! vertical diffusion
                            CALL dyn_spg_init   ! surface pressure gradient

      !                                     ! Misc. options
      ! ADD PALEORCA2 configuration -- JBL 07.02.2017
      IF( nn_cla == 1 .AND. cp_cfg == 'orca' .AND. jp_cfg == 2 )   CALL cla_init       ! Cross Land Advection
      IF( nn_cla == 1 .AND. cp_cfg == 'paleorca' .AND. jp_cfg == 2 )   CALL cla_init       ! Cross Land Advection if needed
                            CALL icb_init( rdt, nit000)   ! initialise icebergs instance
                            CALL sto_par_init   ! Stochastic parametrization
      IF( ln_sto_eos     )  CALL sto_pts_init   ! RRandom T/S fluctuations
     
#if defined key_top
      !                                     ! Passive tracers
                            CALL     trc_init
#endif
      !                                     ! Diagnostics
      IF( lk_floats     )   CALL     flo_init   ! drifting Floats
                            CALL dia_ptr_init   ! Poleward TRansports initialization
      IF( lk_diadct     )   CALL dia_dct_init   ! Sections tranports
                            CALL dia_hsb_init   ! heat content, salt content and volume budgets
                            CALL     trd_init   ! Mixed-layer/Vorticity/Integral constraints trends
      IF( lk_diaobs     ) THEN                  ! Observation & model comparison
                            CALL dia_obs_init            ! Initialize observational data
                            CALL dia_obs( nit000 - 1 )   ! Observation operator for restart
      ENDIF

      !                                     ! Assimilation increments
      IF( lk_asminc     )   CALL asm_inc_init   ! Initialize assimilation increments
      IF(lwp) WRITE(numout,*) 'Euler time step switch is ', neuler
      !
   END SUBROUTINE nemo_init


   SUBROUTINE nemo_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_ctl  ***
      !!
      !! ** Purpose :   control print setting
      !!
      !! ** Method  : - print namctl information and check some consistencies
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'nemo_ctl: Control prints & Benchmark'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namctl'
         WRITE(numout,*) '      run control (for debugging)     ln_ctl     = ', ln_ctl
         WRITE(numout,*) '      level of print                  nn_print   = ', nn_print
         WRITE(numout,*) '      Start i indice for SUM control  nn_ictls   = ', nn_ictls
         WRITE(numout,*) '      End i indice for SUM control    nn_ictle   = ', nn_ictle
         WRITE(numout,*) '      Start j indice for SUM control  nn_jctls   = ', nn_jctls
         WRITE(numout,*) '      End j indice for SUM control    nn_jctle   = ', nn_jctle
         WRITE(numout,*) '      number of proc. following i     nn_isplt   = ', nn_isplt
         WRITE(numout,*) '      number of proc. following j     nn_jsplt   = ', nn_jsplt
         WRITE(numout,*) '      benchmark parameter (0/1)       nn_bench   = ', nn_bench
         WRITE(numout,*) '      timing activated    (0/1)       nn_timing  = ', nn_timing
      ENDIF
      !
      nprint    = nn_print          ! convert DOCTOR namelist names into OLD names
      nictls    = nn_ictls
      nictle    = nn_ictle
      njctls    = nn_jctls
      njctle    = nn_jctle
      isplt     = nn_isplt
      jsplt     = nn_jsplt
      nbench    = nn_bench

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'namcfg  : configuration initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namcfg'
         WRITE(numout,*) '      configuration name              cp_cfg      = ', TRIM(cp_cfg)
         WRITE(numout,*) '      configuration zoom name         cp_cfz      = ', TRIM(cp_cfz)
         WRITE(numout,*) '      configuration resolution        jp_cfg      = ', jp_cfg
         WRITE(numout,*) '      1st lateral dimension ( >= jpi ) jpidta     = ', jpidta
         WRITE(numout,*) '      2nd    "         "    ( >= jpj ) jpjdta     = ', jpjdta
         WRITE(numout,*) '      3nd    "         "               jpkdta     = ', jpkdta
         WRITE(numout,*) '      1st dimension of global domain in i jpiglo  = ', jpiglo
         WRITE(numout,*) '      2nd    -                  -    in j jpjglo  = ', jpjglo
         WRITE(numout,*) '      left bottom i index of the zoom (in data domain) jpizoom = ', jpizoom
         WRITE(numout,*) '      left bottom j index of the zoom (in data domain) jpizoom = ', jpjzoom
         WRITE(numout,*) '      lateral cond. type (between 0 and 6) jperio = ', jperio   
         WRITE(numout,*) '      use file attribute if exists as i/p j-start ln_use_jattr = ', ln_use_jattr
      ENDIF
      !                             ! Parameter control
      !
      IF( ln_ctl ) THEN                 ! sub-domain area indices for the control prints
         IF( lk_mpp .AND. jpnij > 1 ) THEN
            isplt = jpni   ;   jsplt = jpnj   ;   ijsplt = jpni*jpnj   ! the domain is forced to the real split domain
         ELSE
            IF( isplt == 1 .AND. jsplt == 1  ) THEN
               CALL ctl_warn( ' - isplt & jsplt are equal to 1',   &
                  &           ' - the print control will be done over the whole domain' )
            ENDIF
            ijsplt = isplt * jsplt            ! total number of processors ijsplt
         ENDIF
         IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
         IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt
         !
         !                              ! indices used for the SUM control
         IF( nictls+nictle+njctls+njctle == 0 )   THEN    ! print control done over the default area
            lsp_area = .FALSE.
         ELSE                                             ! print control done over a specific  area
            lsp_area = .TRUE.
            IF( nictls < 1 .OR. nictls > jpiglo )   THEN
               CALL ctl_warn( '          - nictls must be 1<=nictls>=jpiglo, it is forced to 1' )
               nictls = 1
            ENDIF
            IF( nictle < 1 .OR. nictle > jpiglo )   THEN
               CALL ctl_warn( '          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo' )
               nictle = jpiglo
            ENDIF
            IF( njctls < 1 .OR. njctls > jpjglo )   THEN
               CALL ctl_warn( '          - njctls must be 1<=njctls>=jpjglo, it is forced to 1' )
               njctls = 1
            ENDIF
            IF( njctle < 1 .OR. njctle > jpjglo )   THEN
               CALL ctl_warn( '          - njctle must be 1<=njctle>=jpjglo, it is forced to jpjglo' )
               njctle = jpjglo
            ENDIF
         ENDIF
      ENDIF
      !
      IF( nbench == 1 ) THEN              ! Benchmark
         SELECT CASE ( cp_cfg )
         CASE ( 'gyre' )   ;   CALL ctl_warn( ' The Benchmark is activated ' )
         CASE DEFAULT      ;   CALL ctl_stop( ' The Benchmark is based on the GYRE configuration:',   &
            &                                 ' cp_cfg = "gyre" in namelist &namcfg or set nbench = 0' )
         END SELECT
      ENDIF
      !
      IF( 1_wp /= SIGN(1._wp,-0._wp)  )   CALL ctl_stop( 'nemo_ctl: The intrinsec SIGN function follows ',  &
         &                                               'f2003 standard. '                              ,  &
         &                                               'Compile with key_nosignedzero enabled' )
      !
   END SUBROUTINE nemo_ctl


   SUBROUTINE nemo_closefile
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_closefile  ***
      !!
      !! ** Purpose :   Close the files
      !!----------------------------------------------------------------------
      !
      IF( lk_mpp )   CALL mppsync
      !
      CALL iom_close                                 ! close all input/output files managed by iom_*
      !
      IF( numstp          /= -1 )   CLOSE( numstp          )   ! time-step file
      IF( numsol          /= -1 )   CLOSE( numsol          )   ! solver file
      IF( numnam_ref      /= -1 )   CLOSE( numnam_ref      )   ! oce reference namelist
      IF( numnam_cfg      /= -1 )   CLOSE( numnam_cfg      )   ! oce configuration namelist
      IF( lwm.AND.numond  /= -1 )   CLOSE( numond          )   ! oce output namelist
      IF( numnam_ice_ref  /= -1 )   CLOSE( numnam_ice_ref  )   ! ice reference namelist
      IF( numnam_ice_cfg  /= -1 )   CLOSE( numnam_ice_cfg  )   ! ice configuration namelist
      IF( lwm.AND.numoni  /= -1 )   CLOSE( numoni          )   ! ice output namelist
      IF( numevo_ice      /= -1 )   CLOSE( numevo_ice      )   ! ice variables (temp. evolution)
      IF( numout          /=  6 )   CLOSE( numout          )   ! standard model output file
      IF( numdct_vol      /= -1 )   CLOSE( numdct_vol      )   ! volume transports
      IF( numdct_heat     /= -1 )   CLOSE( numdct_heat     )   ! heat transports
      IF( numdct_salt     /= -1 )   CLOSE( numdct_salt     )   ! salt transports

      !
      numout = 6                                     ! redefine numout in case it is used after this point...
      !
   END SUBROUTINE nemo_closefile


   SUBROUTINE nemo_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE diawri    , ONLY: dia_wri_alloc
      USE dom_oce   , ONLY: dom_oce_alloc
      USE ldfdyn_oce, ONLY: ldfdyn_oce_alloc
      USE ldftra_oce, ONLY: ldftra_oce_alloc
      USE trc_oce   , ONLY: trc_oce_alloc
#if defined key_diadct 
      USE diadct    , ONLY: diadct_alloc 
#endif 
#if defined key_bdy
      USE bdy_oce   , ONLY: bdy_oce_alloc
#endif
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        oce_alloc       ()          ! ocean
      ierr = ierr + dia_wri_alloc   ()
      ierr = ierr + dom_oce_alloc   ()          ! ocean domain
      ierr = ierr + ldfdyn_oce_alloc()          ! ocean lateral  physics : dynamics
      ierr = ierr + ldftra_oce_alloc()          ! ocean lateral  physics : tracers
      ierr = ierr + zdf_oce_alloc   ()          ! ocean vertical physics
      !
      ierr = ierr + trc_oce_alloc   ()          ! shared TRC / TRA arrays
      !
#if defined key_diadct 
      ierr = ierr + diadct_alloc    ()          ! 
#endif 
#if defined key_bdy
      ierr = ierr + bdy_oce_alloc   ()          ! bdy masks (incl. initialization)
#endif
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'nemo_alloc : unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE nemo_alloc


   SUBROUTINE nemo_partition( num_pes )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE nemo_partition  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   num_pes   ! The number of MPI processes we have
      !
      INTEGER, PARAMETER :: nfactmax = 20
      INTEGER :: nfact ! The no. of factors returned
      INTEGER :: ierr  ! Error flag
      INTEGER :: ji
      INTEGER :: idiff, mindiff, imin ! For choosing pair of factors that are closest in value
      INTEGER, DIMENSION(nfactmax) :: ifact ! Array of factors
      !!----------------------------------------------------------------------
      !
      ierr = 0
      !
      CALL factorise( ifact, nfactmax, nfact, num_pes, ierr )
      !
      IF( nfact <= 1 ) THEN
         WRITE (numout, *) 'WARNING: factorisation of number of PEs failed'
         WRITE (numout, *) '       : using grid of ',num_pes,' x 1'
         jpnj = 1
         jpni = num_pes
      ELSE
         ! Search through factors for the pair that are closest in value
         mindiff = 1000000
         imin    = 1
         DO ji = 1, nfact-1, 2
            idiff = ABS( ifact(ji) - ifact(ji+1) )
            IF( idiff < mindiff ) THEN
               mindiff = idiff
               imin = ji
            ENDIF
         END DO
         jpnj = ifact(imin)
         jpni = ifact(imin + 1)
      ENDIF
      !
      jpnij = jpni*jpnj
      !
   END SUBROUTINE nemo_partition


   SUBROUTINE factorise( kfax, kmaxfax, knfax, kn, kerr )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE factorise  ***
      !!
      !! ** Purpose :   return the prime factors of n.
      !!                knfax factors are returned in array kfax which is of
      !!                maximum dimension kmaxfax.
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kn, kmaxfax
      INTEGER                    , INTENT(  out) ::   kerr, knfax
      INTEGER, DIMENSION(kmaxfax), INTENT(  out) ::   kfax
      !
      INTEGER :: ifac, jl, inu
      INTEGER, PARAMETER :: ntest = 14
      INTEGER, DIMENSION(ntest) :: ilfax
      !
      ! ilfax contains the set of allowed factors.
      ilfax(:) = (/(2**jl,jl=ntest,1,-1)/)

      ! Clear the error flag and initialise output vars
      kerr = 0
      kfax = 1
      knfax = 0

      ! Find the factors of n.
      IF( kn == 1 )   GOTO 20

      ! nu holds the unfactorised part of the number.
      ! knfax holds the number of factors found.
      ! l points to the allowed factor list.
      ! ifac holds the current factor.

      inu   = kn
      knfax = 0

      DO jl = ntest, 1, -1
         !
         ifac = ilfax(jl)
         IF( ifac > inu )   CYCLE

         ! Test whether the factor will divide.

         IF( MOD(inu,ifac) == 0 ) THEN
            !
            knfax = knfax + 1            ! Add the factor to the list
            IF( knfax > kmaxfax ) THEN
               kerr = 6
               write (*,*) 'FACTOR: insufficient space in factor array ', knfax
               return
            ENDIF
            kfax(knfax) = ifac
            ! Store the other factor that goes with this one
            knfax = knfax + 1
            kfax(knfax) = inu / ifac
            !WRITE (*,*) 'ARPDBG, factors ',knfax-1,' & ',knfax,' are ', kfax(knfax-1),' and ',kfax(knfax)
         ENDIF
         !
      END DO

   20 CONTINUE      ! Label 20 is the exit point from the factor search loop.
      !
   END SUBROUTINE factorise

#if defined key_mpp_mpi

   SUBROUTINE nemo_northcomms
      !!======================================================================
      !!                     ***  ROUTINE  nemo_northcomms  ***
      !! nemo_northcomms    :  Setup for north fold exchanges with explicit 
      !!                       point-to-point messaging
      !!=====================================================================
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   Initialization of the northern neighbours lists.
      !!----------------------------------------------------------------------
      !!    1.0  ! 2011-10  (A. C. Coward, NOCS & J. Donners, PRACE)
      !!    2.0  ! 2013-06 Setup avoiding MPI communication (I. Epicoco, S. Mocavero, CMCC) 
      !!----------------------------------------------------------------------

      INTEGER  ::   sxM, dxM, sxT, dxT, jn
      INTEGER  ::   njmppmax

      njmppmax = MAXVAL( njmppt )
    
      !initializes the north-fold communication variables
      isendto(:) = 0
      nsndto = 0

      !if I am a process in the north
      IF ( njmpp == njmppmax ) THEN
          !sxM is the first point (in the global domain) needed to compute the
          !north-fold for the current process
          sxM = jpiglo - nimppt(narea) - nlcit(narea) + 1
          !dxM is the last point (in the global domain) needed to compute the
          !north-fold for the current process
          dxM = jpiglo - nimppt(narea) + 2

          !loop over the other north-fold processes to find the processes
          !managing the points belonging to the sxT-dxT range
  
          DO jn = 1, jpni
                !sxT is the first point (in the global domain) of the jn
                !process
                sxT = nfiimpp(jn, jpnj)
                !dxT is the last point (in the global domain) of the jn
                !process
                dxT = nfiimpp(jn, jpnj) + nfilcit(jn, jpnj) - 1
                IF ((sxM .gt. sxT) .AND. (sxM .lt. dxT)) THEN
                   nsndto = nsndto + 1
                     isendto(nsndto) = jn
                ELSEIF ((sxM .le. sxT) .AND. (dxM .ge. dxT)) THEN
                   nsndto = nsndto + 1
                     isendto(nsndto) = jn
                ELSEIF ((dxM .lt. dxT) .AND. (sxT .lt. dxM)) THEN
                   nsndto = nsndto + 1
                     isendto(nsndto) = jn
                END IF
          END DO
          nfsloop = 1
          nfeloop = nlci
          DO jn = 2,jpni-1
           IF(nfipproc(jn,jpnj) .eq. (narea - 1)) THEN
              IF (nfipproc(jn - 1 ,jpnj) .eq. -1) THEN
                 nfsloop = nldi
              ENDIF
              IF (nfipproc(jn + 1,jpnj) .eq. -1) THEN
                 nfeloop = nlei
              ENDIF
           ENDIF
        END DO

      ENDIF
      l_north_nogather = .TRUE.
   END SUBROUTINE nemo_northcomms
#else
   SUBROUTINE nemo_northcomms      ! Dummy routine
      WRITE(*,*) 'nemo_northcomms: You should not have seen this print! error?'
   END SUBROUTINE nemo_northcomms
#endif

   !!======================================================================
END MODULE nemogcm


