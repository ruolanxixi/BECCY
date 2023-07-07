! Description: COSMO orbital forcing calculations for paleo-climate
!              (and present-day)
!
! Compilation:
! f2py -c --fcompiler=gnu95 --opt=-O3 cosmo_paleo.f95 -m cosmo_paleo
!
! Use in Python: import cosmo_paleo
!                 cosmo_paleo.orbit()
!
! Create conda environment with Fortran compiler:
! conda create -n fortran \
!   numpy matplotlib skyfield fortran-compiler ipython -c conda-forge
!
! Author: Christian Steger, July 2023

!##############################################################################
!# Earth-sun orbit
!##############################################################################

SUBROUTINE orbit(jj, itaja, zstunde, &
                lpalorb, itype_calendar, recc, robld, rlonp, &
                zdek, zsocof, zsct_save, zeit0)

!! Input:
! jj:               year
! itaja:            day of the year
! zstunde           actual hour of the day (UTC)
! lpalorb:          if .TRUE., paleo orbit routine is used
! itype_calendar:   0: gregorian calendar (default), 1: 360 days, 2: 365 days
! recc:             Eccentricitity for paleo orbit
! robld:            Obliquity for paleo orbit
! rlonp:            Longitude of perihelion

!! Output:
! zdek:             solar declination angle [rad] (?)
! zsocof:           solar constant scaling factor [-]
! zsct_save:        solar constant depending on Earth-Sun-distance [W m-2]
! zeit0:            

IMPLICIT NONE

INTEGER :: jj, itaja, itype_calendar
REAL :: zstunde
LOGICAL :: lpalorb
REAL :: yearl, api, zve00, zvebas, zyearve, zyeardif, zone, zvetim
REAL :: zecc, zobld, zlonp
REAL :: recc, robld, rlonp
REAL :: zoblr, zlonpr, zsqecc, zeps, zeve, ztlonpr, zztime, zeold, zenew
INTEGER :: iter
REAL :: zzeps
CHARACTER (LEN=80)   yzerrmsg      ! for error message
REAL :: zcose
REAL :: zsocof, ztgean, znu, zlambda, zsinde, zdek
REAL :: zdeksin_save, zdekcos_save, zsct_save
REAL :: solc, zsct_h
INTEGER :: nz_zsct
REAL :: ztwo, ztho, zdtzgl_save
REAL :: pi, z_1d7200, dt, zeit0

!! Constants
solc     =  1368.0  ! solar constant [W m-2]
pi = 4.0 * ATAN (1.0)
z_1d7200 = 1.0/7200.0
dt = 0.0  ! model time step [seconds]

! Directives for f2py (ignored by compiler but read by f2py) 
!f2py threadsafe
!f2py intent(in) jj, itaja, zstunde, lpalorb, itype_calendar, recc, robld, rlonp
!f2py intent(out) zdek, zsocof, zsct_save, zeit0

          !! PL 07/2020 start orbital routine for paleo
          IF (lpalorb) THEN
             if (itype_calendar==0) then
                yearl = 365 + IABS( MOD(jj,4) - 4) / 4
             else	
                yearl = 360
             end if

             api    = 2.*ASIN(1.)  !Pi
             zve00 = 80.5          !Day of the vernal equinox in 1900
             zvebas = yearl-zve00  !Remainder of the year
             zyearve = float(itaja) - 1. + zvebas   !itaja: actual day of the year; zyearve: time of the year from last years vernal equinox in degrees
             !     zvetim = MOD(zyearve/yearl,1.)*2.*api  !Time of the year from the vernal equinox in radians
             zyeardif=zyearve/yearl
             zone=1.
             zvetim = MOD(zyeardif,zone)*2.*api

             ! these are the values for the preindustrial (present-day) case:
             !ZECC    = 0.016715
             !ZOBLD   = 23.441
             !ZLONP   = 282.7
             zecc=recc
             zobld=robld
             zlonp=rlonp
             !
             zoblr=zobld*api/180.
             zlonpr=zlonp*api/180.
             zsqecc=SQRT((1+zecc)/(1-zecc))
             !
             zeps=1.E-6
             !
             !     CALCULATION OF ECCENTRIC ANOMALY OF VERNAL EQUINOX  !
             zeve=2.*ATAN(TAN(0.5*zlonpr)/zsqecc)
             !
             !     CALCULATION OF TIME ANGLE IN RADIANS OF LONGITUDE OF PERIHELION
             !       FROM VERNAL EQUINOX
             !
             ztlonpr= zeve - zecc*SIN(zeve)
             !
             !     CALCULATE ECCENTRIC ANOMALY:
             !     USE FIRST DERIVATIVE (NEWTON) TO CALCULATE SOLUTION FOR
             !     EQUATION OF ECCENTRIC ANOMALY *ZENEW*
             !
             zztime=zvetim-ztlonpr
             !
             zeold=zztime/(1.-zecc)
             zenew=zztime
             iter=0
             !
             DO
                zzeps=zeold-zenew
                IF (iter.GE.10) THEN
                   yzerrmsg = 'PALEO ORBIT:  - eccentric anomaly not found!'
                   CALL ABORT  !CALL model_abort (my_world_id, 2091, yzerrmsg, 'radiation')
                END IF
                IF (ABS(zzeps).LT.zeps) EXIT
                iter=iter+1
                zeold=zenew
                zcose=COS(zenew)
                zenew=(zztime+zecc*(SIN(zenew)-zenew*zcose))/(1.-zecc*zcose)
             END DO
           !
             zsocof=(1./(1.-zecc*COS(zenew)))**2
             !
             !     CALCULATION OF THE DECLINATION.
             !
             ztgean=TAN(zenew*0.5)
             !
             !     *znu*: TRUE ANOMALY   !            (ACTUAL ANGLE OF EARTH'S POSITION FROM PERIHELION)
             !     *zlambda*: TRUE LONGITUDE OF THE EARTH     !                (ACTUAL ANGLE FROM VERNAL EQUINOX)
             !
             znu=2.*ATAN(zsqecc*ztgean)
             zlambda=znu+zlonpr
             zsinde=SIN(zoblr)*SIN(zlambda)
             zdek=ASIN(zsinde)
             
             zdeksin_save = SIN (zdek)
             zdekcos_save = COS (zdek)
             
             zsct_save = zsocof*solc
             
             zsct_h = zsct_h + zsct_save !mkk da in orig routine
             nz_zsct = nz_zsct + 1 !mkk da in orig routine
             !print*,"ORB",zdeksin_save,zdekcos_save,zsct_save,zsct_h,nz_zsct
             ztwo    = 0.681 + 0.2422*(jj-1949)-(jj-1949)/4 !mkk, da fuer zeit0 benoetigt und in orig schleife
             ztho    = 2.*pi*( REAL(itaja) -1.0 + ztwo )/365.2422 !mkk, da fuer zeit0 benoetigt und in orig schleife 
             zdtzgl_save  = 0.000075 + 0.001868*COS(   ztho) - 0.032077*SIN(   ztho) & !mkk, da fuer zeit0 benoetigt und in orig schleife
                  - 0.014615*COS(2.*ztho) - 0.040849*SIN(2.*ztho)

          ELSE ! end of paleo orbit
 !! PL 07/2020

             ztwo    = 0.681 + 0.2422*(jj-1949)-(jj-1949)/4
             ztho    = 2.*pi*( REAL(itaja) -1.0 + ztwo )/365.2422
             zdtzgl_save  = 0.000075 + 0.001868*COS(   ztho) - 0.032077*SIN(   ztho) &
                  - 0.014615*COS(2.*ztho) - 0.040849*SIN(2.*ztho)
             zdek    = 0.006918 - 0.399912*COS(   ztho) + 0.070257*SIN(   ztho) &
                  - 0.006758*COS(2.*ztho) + 0.000907*SIN(2.*ztho) &
                  - 0.002697*COS(3.*ztho) + 0.001480*SIN(3.*ztho)
             
             zdeksin_save = SIN (zdek)
             zdekcos_save = COS (zdek)
             
             zsocof  = 1.000110 + 0.034221*COS(   ztho) + 0.001280*SIN(   ztho) &
                  + 0.000719*COS(2.*ztho) + 0.000077*SIN(2.*ztho)
             
             zsct_save = zsocof*solc
             zsct_h = zsct_h + zsct_save
             nz_zsct = nz_zsct + 1

          ENDIF

       zstunde = zstunde + dt*z_1d7200 !add half a timestep
       zeit0   = pi*(zstunde-12.0)/12.0 + zdtzgl_save

END SUBROUTINE orbit

!##############################################################################