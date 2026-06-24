
!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   H E A T  E X C H A N G E                                         **
!***********************************************************************************************************************************

SUBROUTINE HEAT_EXCHANGE
  USE GLOBAL; USE GDAYC; USE SURFHE; USE TVDC; USE SHADEC; USE PREC; USE MetFileRegion; USE MAIN, ONLY:CC_SW, I_SW
  IMPLICIT NONE

! Type declaration
  REAL     :: JDAY, LOCAL
  REAL(R8) :: TSUR
  REAL(R8) :: MPS_TO_MPH=2.23714D0,W_M2_TO_BTU_FT2_DAY=7.60796D0,FLUX_BR_TO_FLUX_SI=0.23659D0,BOWEN_CONSTANT=0.47D0
  REAL(R8) :: DEG_F,DEG_C,BTU_FT2_DAY_TO_W_M2=0.1314D0
  REAL(R8) :: STANDARD,HOUR,TAUD,SINAL,A0,TDEW_F,TAIR_F,SRO_BR,WIND_MPH,WIND2M,ACONV,BCONV
  REAL(R8) :: TSTAR,BETA,FW,RA,ETP,EA,ES,TAIRV,DTV,DTVL
  INTEGER  :: IDAY,J

RETURN

!***********************************************************************************************************************************
!**                                            S H O R T  W A V E  R A D I A T I O N                                              **
!***********************************************************************************************************************************

ENTRY SHORT_WAVE_RADIATION (JDAY)
  IF(I_SW==1)THEN
      CALL MEEUS
  ELSE
      
  LOCAL    =  LONGIT(JW)
  STANDARD =  15.0*INT(LONGIT(JW)/15.0)
  HOUR     = (JDAY-INT(JDAY))*24.0
  IDAY     =  JDAY-((INT(JDAY/365))*365)
  !IDAY     =  IDAY+INT(INT(JDAY/365)/4)  ! SR 12/2018 IDAY FIX
  IDAY     =  IDAY-INT(INT(JDAY/365)/4)
  TAUD     = (2.*PI*(IDAY-1))/365.
  EQTNEW   =  0.170*SIN(4.*PI*(IDAY-80)/373.)-0.129*SIN(2.*PI*(IDAY-8)/355.)
  HH(JW)   =  0.261799*(HOUR-(LOCAL-STANDARD)*0.0666667+EQTNEW-12.0)
  DECL(JW) =  0.006918-0.399912*COS(TAUD)+0.070257*SIN(TAUD)-0.006758*COS(2.*TAUD)+0.000907*SIN(2.*TAUD)-0.002697*COS(3.*TAUD)        &
              +0.001480*SIN(3.*TAUD)
  SINAL    =  SIN(LAT(JW)*.0174533)*SIN(DECL(JW))+COS(LAT(JW)*.0174533)*COS(DECL(JW))*COS(HH(JW))
  A00(JW)  =  57.2957795*ASIN(SINAL)
  A0       =  A00(JW)
  if(.not.Met_Regions)then
  IF (A0 > 0.0) THEN
    SRON(JW) = (1.0-CC_SW*CLOUD(JW)*CLOUD(JW))*24.0*(2.044*A0+0.1296*A0**2-1.941E-3*A0**3+7.591E-6*A0**4)*BTU_FT2_DAY_TO_W_M2
  ELSE
    SRON(JW) = 0.0
  END IF
  else
  IF (A0 > 0.0) THEN
    SRON(NMet) = (1.0-CC_SW*CLOUD(NMet)*CLOUD(NMet))*24.0*(2.044*A0+0.1296*A0**2-1.941E-3*A0**3+7.591E-6*A0**4)*BTU_FT2_DAY_TO_W_M2
  ELSE
    SRON(NMet) = 0.0
  END IF
  endif
  ENDIF
  
  
RETURN

!***********************************************************************************************************************************
!**                                          E Q U I L I B R I U M  T E M P E R A T U R E                                         **
!***********************************************************************************************************************************

ENTRY EQUILIBRIUM_TEMPERATURE

! British units
  if(.not.Met_Regions)then
  TDEW_F   = DEG_F(TDEW(JW))
  TAIR_F   = DEG_F(TAIR(JW))
  SRO_BR   = SRON(JW)*W_M2_TO_BTU_FT2_DAY*SHADE(I)
  else
  TDEW_F   = DEG_F(TDEW(NMet))
  TAIR_F   = DEG_F(TAIR(NMet))
  SRO_BR   = SRON(NMet)*W_M2_TO_BTU_FT2_DAY*SHADE(I)   
  endif
  
  !WIND_MPH = WIND(JW)*WSC(I)*MPS_TO_MPH
  !WIND2M   = WIND_MPH*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))+NONZERO     ! SW 11/28/07  old version z0=0.003
  WIND2M=WIND2(I)*MPS_TO_MPH    ! ALREADY COMPUTED IN w2 MAIN PROGRAM
  ACONV    = W_M2_TO_BTU_FT2_DAY
  IF (CFW(JW) == 1.0) BCONV = 3.401062
  IF (CFW(JW) == 2.0) BCONV = 1.520411
  ! SW Issues: CFW not determined for other values of CFW

! Equilibrium temperature and heat exchange coefficient

  ET(I)   =  TDEW_F
  TSTAR   = (ET(I)+TDEW_F)*0.5
  BETA    =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
  FW      =  ACONV*AFW(JW)+BCONV*BFW(JW)*WIND2M**CFW(JW)
  CSHE(I) =  15.7+(0.26+BETA)*FW
  RA      =  3.1872E-08*(TAIR_F+459.67)**4   ! original is missing the effect of cloud cover
  ETP     = (SRO_BR+RA-1801.0)/CSHE(I)+(CSHE(I)-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE(I)*(0.26+BETA))
  J       =  0
  DO WHILE (ABS(ETP-ET(I)) > 0.05 .AND. J < 10)
    ET(I)   =  ETP
    TSTAR   = (ET(I)+TDEW_F)*0.5
    BETA    =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
    CSHE(I) =  15.7+(0.26+BETA)*FW
    ETP     = (SRO_BR+RA-1801.0)/CSHE(I)+(CSHE(I)-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE(I)*(0.26+BETA))
    J       =  J+1
  END DO

! SI units

  ET(I)   = DEG_C(ET(I))
  CSHE(I) = CSHE(I)*FLUX_BR_TO_FLUX_SI/RHOWCP
RETURN

!***********************************************************************************************************************************
!**                                                   S U R F A C E   T E R M S                                                   **
!***********************************************************************************************************************************

ENTRY SURFACE_TERMS (TSUR)

! Partial water vapor pressure of air (mm hg)

IF(.NOT.Met_regions)then
  EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))
else
  EA = DEXP(2.3026D0*(7.5D0*TDEW(NMet)/(TDEW(NMet)+237.3D0)+0.6609D0))
endif


! Partial water vapor pressure at the water surface

  IF(TSUR<0.0)THEN
  ES = DEXP(2.3026D0*(9.5D0*TSUR/(TSUR+265.5D0)+0.6609D0))
  ELSE
  ES = DEXP(2.3026D0*(7.5D0*TSUR/(TSUR+237.3D0)+0.6609D0))
  ENDIF

! Wind function

  IF (RH_EVAP(JW)) THEN
    if(.not.Met_Regions)then
    TAIRV = (TAIR(JW)+273.0D0)/(1.0D0-0.378D0*EA/760.0D0)
    else
    TAIRV = (TAIR(NMet)+273.0D0)/(1.0D0-0.378D0*EA/760.0D0)
    endif
    DTV   = (TSUR+273.0D0)/(1.0D0-0.378D0*ES/760.0D0)-TAIRV
    DTVL  =  0.0084D0*WIND2(I)**3
    IF (DTV < DTVL) DTV = DTVL
    FW = (3.59D0*DTV**0.3333D0+4.26D0*WIND2(I))
  ELSE
    FW = AFW(JW)+BFW(JW)*WIND2(I)**CFW(JW)
  END IF

! Evaporative flux

  RE(I) = FW*(ES-EA)
  !IF(RE(I) < 0.0)RE(I)=0.0     ! SW 6/22/2016  SHOULD WE USE THIS AS AN ANALOG FOR CONDENSATION WHEN LESS THAN 0? TVA(1972) SUGGESTS YOU CAN - SEE SECTION 4 AND 5 but Ryan, Harleman suggest setting RE=0 if less than zero since the condensation coefficients are unknown...

! Conductive flux
 if(.not.Met_Regions)then
  RC(I) = FW*BOWEN_CONSTANT*(TSUR-TAIR(JW))
 else
  RC(I) = FW*BOWEN_CONSTANT*(TSUR-TAIR(NMet))  
 endif
 

! Back radiation flux

  RB(I) = 5.51D-8*(TSUR+273.15D0)**4
  RETURN  
END SUBROUTINE HEAT_EXCHANGE

! Function declaration

   REAL(R8) FUNCTION DEG_F(X)
   USE PREC
   REAL(R8) :: X
   DEG_F =  X*1.8+32.0
   END FUNCTION DEG_F
   REAL(R8) FUNCTION DEG_C(X)
   USE PREC
   REAL(R8) :: X
   DEG_C = (X-32.0)*5.0/9.0
    END FUNCTION DEG_C
    
Module Bird_Meeus
use prec
    REAL(R8) :: B0, B1, BB, Beta, jme,rc1
    REAL(R8) :: A9,B9, jd, t, LO, M, Mrad, sinM, sin2M, sin3M, c, TLO, v, Selev, Te
    REAL(R8) :: omega, lambda, seconds, e0, e, tananum, tanadenom, alpha, sint, theta, ep,long2(100)
    REAL(R8) :: y, sin2LO, cos2LO, sin4LO, Etime, EQT, solarazimuth
    REAL(R8) :: day_decimal, HOUR_M, JDAYM                   ! real HOUR, day_decimal,Jdaym
    REAL(R8) :: Zone, daySavings, hh, mm, ss, timenow
    REAL(R8) :: solarTimeFix, trueSolarTime, hourangle, harad, csz, zenith, azDenom, azRad
    REAL(R8) :: azimuth, exoatmElevation, step1, step2, step3, refractionCorrection
    Real(R8) :: r, solarzen, RToD=57.2957796, DToR=0.017453293, solaralt, Uo     ! RToD = (180.0/PI)      DToR = (PI/180.0)
    Real(R8) :: directHz, diffuseHz, globalHz
    Integer  :: year9,month9                 !Debug_Meeus=1, 
    
End Module Bird_Meeus

SUBROUTINE MEEUS
  USE GLOBAL; USE GDAYC; USE SURFHE; USE TVDC; USE SHADEC; USE MetFileRegion; USE MAIN, ONLY:CC_SW,TZONE
  USE Bird_Meeus;use screenc, only: JDAY
  IMPLICIT NONE

! Calculates Solar position and radiation from Meeus code, currently set up for clear skies only
      !if(debug_Meeus==1.and.NIT==0)then
      !  open(9595,file='meeus_debug.csv',status='unknown')
      !  write(9595,'(A)')'jday,jdaym,t,jd,globalHz,solarzen,RC1,Selev,azimuth,azRad,csz,hourangle,ETime,y,LO,r,ep,EQT,day_decimal,year,imon,jdayg,zenith'
      !  open(9596,file='meeus_debug2.csv',status='unknown')
      !  write(9596,'(A)')'jdaym,Uo,AM,diffuseHz,global_hz,direct_Hz,Ias,toh,tum,tw,taa,am,tr,ta,birdetr,rs,Uw,tdew,sw_BA,sw_albedo,ELWS(DS(BE(JW))),sw_k1'
      !endif
      
    Jdaym=JDAY-TZONE/24.0

    HOUR_M  =  (Jdaym-INT(Jdaym))*24        
    day_decimal = GDAY + HOUR_M/24.0
    
    if(IMON <= 2)then
    year9 = year - 1
    month9 = IMON + 12
    else
    year9=year
    month9=imon
    end if   

    A9 = FLOOR(year9/100.0)   
    B9 = 2.0 - A9 + FLOOR(A9/4.0)   
    jd = FLOOR(365.25*(year9+4716.0)) + FLOOR(30.6001*(month9+1))+day_decimal + B9 - 1524.5         ! Eq 7.1, pg 61

! time in Julian centuries from epoch J2000.0
      t = (jd-2451545.0)/36525.0                                                               ! Eq 25.1, pg 163
! julian ephemeris millennium
      jme  = t/10.0

! calculate the Geometric Mean Longitude of the Sun, degrees
      LO = 280.46646+t*(36000.76983+0.0003032*t)                                               ! Eq 25.2, pg 163
!      LO = 280.4664567+36000.76982779*t+0.0003032028*(t**2)+(t**3)/49931-(t**4)/15300-(t**5)/2000000 ! Eq 28.2, pg 183, difference with line above is negligible
      DO WHILE(LO > 360 .or. LO < 0)
       If(LO > 360.0) LO = LO - 360
       If(LO < 0.0)   LO = LO + 360
      END DO
! calculate the Geometric Mean Anomaly of the Sun, degrees
      M = 357.52911 + t*(35999.05029-0.0001537*t)                                              ! Eq 25.3, pg 163
! calculate the eccentricity of earth's orbit
      e = 0.016708634 - t*(0.000042037+0.0000001267*t)                                         ! Eq 25.4, pg 163
! calculate the equation of center for the sun, degrees
      Mrad = DToR*M
      sinM = Sin(Mrad)
      sin2M = Sin(2*Mrad)
      sin3M = Sin(3*Mrad)
      c = sinM*(1.914602-t*(0.004817+0.000014*t))+sin2M*(0.019993-0.000101*t)+sin3M*0.000289   ! Between Eq 25.4 and Eq 25.5, pg 164
! calculate the true longitude of the sun, degrees
      TLO = LO + c                                                                             ! Between Eq 25.4 and Eq 25.5, pg 164
! calculate the true anomaly of the sun, degrees
      v = M + c                                                                                ! Between Eq 25.4 and Eq 25.5, pg 164
! calculate the distance to the sun in AU
      r = (1.000001018*(1 - e*e))/(1 + e*Cos(DToR*v))                                          ! Eq 25.5, pg 164
! calculate the apparent longitude of the sun, degrees
      omega = 125.04 - 1934.136*t                                                              ! Between Eq 25.5 and Eq 25.6, pg 164
!      omega  = 125.04452 - 1934.136261*t+0.0020708*(t**2)+(t**3)/45000                        ! Equation on page 144, difference with line above is negligible
      lambda = TLO - 0.00569 - 0.00478*Sin(DToR*omega)                                         ! Between Eq 25.5 and Eq 25.6, pg 164
! calculate the mean obliquity of the ecliptic, degrees
      seconds = 21.448 - t*(46.8150 + t*(0.00059 - 0.001813*t))                                ! Eq 22.2, pg 147
      e0 = 23.0 + (26.0 + (seconds/60.0))/60.0                                                 ! Eq 22.2, pg 147
! calculate the corrected obliquity of the ecliptic, degrees
      ep = e0 + 0.00256 * Cos(DToR*omega)
!      ep = e0 + (9.20/3600.0)*Cos(DToR*omega)+(0.57/3600.0)*cos(DToR*2*LO)                    ! Eq on page 144, still truncated though here in code, difference with line above is negligible

! earth heliocentric latitude, BB, degrees, Equations 9 to 12, NREL/TP-560-34302
      !B0   = 280.0*Cos(3.199+84334.662*jme)+102.0*Cos(5.422+5507.553*jme)+  &
      !        80.0*Cos(3.880+ 5223.690*jme)+ 44.0*Cos(3.700+2352.870*jme)+  &
      !        32.0*Cos(4.000+ 1577.340*jme)
      !B1   =   9.0*Cos(3.900+ 5507.550*jme)+  6.0*Cos(1.730+5223.690*jme)
      !BB   = ((B0+B1*jme)/100000000.0)*RToD
! geocentric latitude, Beta, degrees
!      Beta = -1.0*BB
! calculate the right ascension of the sun, degrees
      tananum = Cos(DToR*ep) * Sin(DToR*lambda)
!      tananum   = Cos(DToR*ep)*Sin(DToR*lambda)-tan(DToR*Beta)*Sin(DToR*ep)                   ! Eq 13.3 on page 93, difference with line above is negligible
      tanadenom = Cos(DToR*lambda)                                                             ! Eq 13.3 on page 93
! Javascript and FORTRAN use same format, Atan2(Y,X)
      alpha = RtoD*atan2(tananum, tanadenom)                                                   ! Eq 13.3 on page 93

! calculate the declination of the sun, degrees
      sint = Sin(DToR*ep) * Sin(DToR*lambda)
!      sint = Sin(DToR*ep)*Sin(DToR*lambda)*Cos(DToR*Beta)+Cos(DToR*ep)*Sin(DToR*Beta)         ! Eq 13.4 on page 93, difference with line above is negligible
      theta = RToD*Asin(sint)                                                                  ! Eq 13.4 on page 93

! calculate the equation of time, the difference between true solar time 
! and mean solar time, in minutes of time
      y = (Tan(DToR*ep/2.0))**2
      sin2LO = Sin(2.0*DToR*LO)
      cos2LO = Cos(2.0*DToR*LO)
      sin4LO = Sin(4.0*DToR*LO)
      Etime = y*sin2LO - 2.0*e*sinM + 4.0*e*y*sinM*cos2LO - 0.5*y*y*sin4LO - 1.25*e*e*sin2M    ! Eq 28.3 on page 185
      EQT = RToD*Etime*4.0

! change sign convention for longitude from negative to positive in western hemisphere  
      if(LongIT(JW) < 0.0)then
       Long2(JW) = -LongIT(JW)
      else
       Long2(JW) = LongIT(JW)
      end if
      if(Lat(JW) > 89.8)  Lat(JW) = 89.8
      if(Lat(JW) < -89.8) Lat(JW) = -89.8

      solarTimeFix = EQT - 4.0*Long2(JW) !+ 60.0*TZ(i)  ! in minutes
      trueSolarTime = HOUR_M*60.0 + solarTimeFix   !in minutes

      if(trueSolarTime > 1440)then
       trueSolarTime = trueSolarTime - 1440
      end if
            
      hourangle = trueSolarTime/4.0 - 180.0

!Thanks to Louis Schwarzmayr for the next line:
      if(hourangle < -180.0)hourangle = hourangle + 360.0

!Solar zenith (Solar altitude)
      harad = DToR*hourangle
      csz = Sin(DToR*Lat(JW))*Sin(DToR*theta) + Cos(DToR*Lat(JW))*Cos(DToR*theta)*Cos(harad)
      if(csz > 1.0) Then
       csz = 1.0
      else if(csz < -1.0) Then
       csz = -1.0
      end if
      zenith = RToD*Acos(csz)

!Solar azimuth
      azDenom = Cos(DToR*Lat(JW))*Sin(DToR*zenith)
      if(Abs(azDenom) > 0.001)then
       azRad = ((Sin(DToR*Lat(JW))*Cos(DToR*zenith)) - Sin(DToR*theta))/azDenom
       if(Abs(azRad) > 1.0)then
        if(azRad < 0.0)then
         azRad = -1.0
        else
         azRad = 1.0
        end if
       end if
       azimuth = 180.0 - RToD*Acos(azRad)
       if(hourangle > 0.0)then
        azimuth = -azimuth
       end if
      else
       if (Lat(JW) > 0.0)then
        azimuth = 180.0
       else
        azimuth = 0.0
       end if
      end if

      if(azimuth < 0.0)then
       azimuth = azimuth + 360.0
      end if
                        
! Correction for the effect of atmospheric refraction, NOAA website (SElev is Solar Altitude)
      SElev = 90.0 - zenith
      if(SElev > 85.0)then
       RC1 = 0.0
      else
       Te = Tan(DToR*SElev)
       if(SElev > 5.0 .and. SElev <= 85.0)   RC1 = (58.1/Te - 0.07/(Te**3) + 0.000086/(Te**5))/3600.0
       if(SElev > -0.575 .and. SElev <= 5.0) RC1 = (1735-518.2*SElev+103.4*SElev**2-12.79*SElev**3+0.711*SElev**4)/3600.0
       if(SElev < -0.575)                    RC1 = (-20.774/Te)/3600.0
      end if

       solarzen     = zenith - RC1                    ! solar zenith, corrected, degrees
       solarazimuth = azimuth                        ! solar azimuth, degrees
       solaralt     = 90.0 - solarzen                ! solar altitude, degrees
       solarzen     = solarzen*DToR                  ! solar zenith, radians

       Call Ozone
       Call BirdModel
       
    if(.not.Met_Regions)then
    SRON(JW) = globalHz*(1.0-CC_SW*CLOUD(JW)*CLOUD(JW))*0.94      ! assumes 6% reflection
    else
    SRON(NMet) = globalHz*(1.0-CC_SW*CLOUD(NMet)*CLOUD(NMet))*0.94
    END IF
  
    
    !if(debug_Meeus==1)then
    !    write(9595,'(19(f14.6,","),i10,",",i10,","i10,",",5(f14.6,","))')jday,jdaym,t,jd,globalHz,solarzen,RC1,Selev,azimuth,azRad,csz,hourangle,ETime,y,LO,r,ep,EQT,day_decimal,year,imon,jdayg,zenith
    !endif
    RETURN

    END SUBROUTINE MEEUS

    SUBROUTINE Ozone
    USE Bird_Meeus
! Calculate ozone concentration, Van Heuklon (1979) Eq 4
     USE GLOBAL; USE GDAYC; USE SURFHE; USE TVDC; USE SHADEC; USE MetFileRegion
  IMPLICIT NONE
    Real A7, B7, C7, F7, H7, P7

    if(Lat(jw) >= 0.0)then
      A7 = 150.0; C7 = 40.0; F7 = -30.000; H7 = 3.0; B7 = 1.28
      if(-Longit(jw) >= 0.0)then   ! Convention in W2 is to use a Long W or a positive #, so actual Longit should be "-" in Western hemisphere, hence the minus sign SW
       P7 = 20.0
      else
       P7 =  0.0
      end if
    else
      A7 = 100.0; C7 = 30.0; F7 = 152.625; H7 = 2.0; B7 = 1.50; P7 = -75.0
    end if
! D converts the number of days to a fraction part of 360 degrees, in degrees
!    D  = 360.0/365.25, 0.98562628
! Uo is calculated in atm-cm
    Uo  = (235+(A7+C7*Sin(DToR*0.9856*(int(Jdaym)+F7))+20.0*Sin(DToR*H7*(Longit(JW)+P7)))*(sin(DToR*B7*Lat(jw))**2))/1000.0

END SUBROUTINE Ozone


SUBROUTINE BirdModel
  USE Bird_Meeus; USE GEOMC
  USE GLOBAL; USE GDAYC; USE SURFHE; USE TVDC; USE SHADEC; USE MetFileRegion; USE MAIN, ONLY:SW_ALBEDO,SW_BA,SW_K1,AOD380NM,AOD500NM
  IMPLICIT NONE

    Real AM, tauA, TR, Toh, TUM, Tw, TA, AMP, Xo,uw,xw
    Real TAA, rs, mp, numI, denI, Ias, BirdETR

! Atmospheric Pressure
!     P         = 1013
! Ratio of forward-scattered irradiance to the total scattered, see SERI pages 8 & 9
!     Ba        = 0.84  --> input variable SW_BA
! Constant used in Bird Model associated with aerosol absorptance, see SERI pages 7 & 9
!     K1        = 0.10   --> input variable SW_K1
! Used same albdeo from other code for Methods 1 to 3
!     albedo    = 0.10   --> input variable SW_ALBEDO
! NREL solar constant of 1367 W/m^2
!     SROCON    = 1367.00
! Bras/TVA (Io=Wo(sin(solaralt))/r2, Equation 2.9 on page 26) adjustment of ETR normal to beam for earth sun distance
     birdETR   = 1367.0/(r*r)
! broadband aerosol optical depth from surface in a vertical path (broadband turbidity)
     tauA      = 0.2758*AOD380nm + 0.35*AOD500nm
!     tauA      = 0.05
    if(RToD*solarzen < 89.0)then
! Relative Optical Air Mass, SERI, pg 8
     AM      = (((288.0-0.0065*ELWS(DS(BE(JW))))/288.)**5.256)/(Sin(DToR*solaralt)+0.1500*(solaralt+3.885)**(-1.253))                    ! Equation 23, optical air mass
! Pressure corrected air mass
     AMP       = AM*1013.0/1013                     !Atm Pressure assumed==1013.
! Total amount of ozone in a slanted path (cm)
     Xo        = Uo*AM
! Amount of Precipitable water in a vertical column  from surface (cm), used 1.50 cm for Mt Vernon, range from MLS and USS models 1.42 to 2.93
     if(.not.Met_Regions)then
     Uw        = exp(-0.0592+0.06912*TDEW(JW))           
     else
     Uw        = exp(-0.0592+0.06912*TDEW(NMET))     
     endif
     ! Equation 27, mean hourly precipitable water content, cm
! Total amount of precipitable water in a slanted path (cm)
     Xw        = Uw*AM
! Transmittance of Rayleigh Scattering, SERI, pg 8
     TR        = Exp(-0.0903*(AMP**0.84)*(1+AMP-AMP**1.01))
! Transmittance of ozone absorptance, SERI, pg 8
     Toh       = 1-0.1611*Xo*(1+139.48*Xo)**-0.3035-0.002715*Xo/(1+0.044*Xo+0.0003*Xo*Xo)
! Transmittance of absorptance of uniformly mixed gases (carbon dioxide and oxygen)
     TUM       = Exp(-0.0127*AMP**0.26)
! Transmittance of water vapor absoptance (1-aw)
     Tw        = 1-2.4959*Xw/((1 + 79.034*Xw)**0.6828+6.385*Xw)
! Transmittance of aerosol aborptance and scattering
     TA        = Exp(-(tauA**0.873)*(1+tauA-tauA**0.7088)*AM**0.9108)
! Transmittance of aerosol aborptance
     TAA       = 1-SW_K1*(1-AM+AM**1.06)*(1-TA)
! Sky, or atmospheric, albedo
     rs        = 0.0685+(1-SW_Ba)*(1-(TA/TAA))
! Direct solar irradiance on a horizontal ground surface (W/m^2)
     if(RToD*solarzen < 90.0)then
      directHz = birdETR*Cos(solarzen)*0.9662*TA*Tw*TUM*Toh*TR
     else
      directHz = 0.0
     end if
! Solar irradiance on a horizontal surface from atmospheric scattering (W/m2)
     numI      = 0.5*(1-TR)+ SW_Ba * (1-(TA/TAA))
     denI      = 1-AM+AM**1.02
     Ias       = birdETR*Cos(solarzen)*0.79*Toh*TUM*Tw*TAA*(numI/DenI)
   else
     AM = 0.0; TR = 0.0; Toh = 0.0; TUM = 0.0; Tw = 0.0; TA = 0.0; TAA = 0.0; rs = 0.0; Ias = 0.0; AMP = 0.0
     Xw = 0.0; Xo = 0.0; mp = 0.0
   end if

!global radiation on a horizontal ground surface (W/m^2)
   if(AM > 0.0)then
    globalHz   = (directHz+Ias)/(1-SW_albedo*rs)
   else
    globalHz   = 0.0
   end if

!diffuse radiation on a horizontal ground surface (W/m^2)
   diffuseHz = globalHz - directHz
   
    !   if(debug_Meeus==1)then
    !    write(9596,'(25(f14.6,","))')jdaym,Uo,AM,diffuseHz,globalhz,directHz,Ias,toh,tum,tw,taa,am,tr,ta,birdetr,rs,Uw,tdew(jw),sw_BA,sw_albedo,ELWS(DS(BE(JW))),sw_k1
    !endif


END SUBROUTINE BirdModel
