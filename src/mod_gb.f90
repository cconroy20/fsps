SUBROUTINE MOD_GB(zz,t,age,delt,dell,pagb,redgb,agb,&
     nn,logl,logt,phase,wght)

  !routine to modify TP-AGB stars, HB+RGB, and post-AGB stars. 

  USE sps_vars
  IMPLICIT NONE

  INTEGER,  INTENT(in) :: t, nn,zz
  REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: logl,logt
  REAL(SP), INTENT(in), DIMENSION(nt,nm)    :: phase
  REAL(SP), INTENT(inout), DIMENSION(nm)    :: wght
  REAL(SP), INTENT(in) :: delt, dell, pagb, redgb, agb
  REAL(SP), INTENT(in), DIMENSION(nt) :: age
  INTEGER  :: i
  REAL(SP) :: age8=8.0_sp,age91=9.1_sp,twght

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  DO i=1,nn

     !modify TP-AGB stars
     IF (phase(t,i).EQ.5.0) THEN

        !renormalize the Padova isochrones
        IF (isoc_type.EQ.'pdva'.AND.tpagb_norm_type.NE.0) THEN

           !Conroy & Gunn (2010) normalization
           IF (tpagb_norm_type.EQ.1) THEN

              IF (age(t).GT.age8.AND.age(t).LT.age91) THEN
                 logl(t,i) = logl(t,i) - 1.0+(age(t)-8.)/1.5
                 IF (LOG10(zlegend(zz)/zsol).LT.-0.25) THEN
                    logt(t,i) = logt(t,i) + 0.10
                 ENDIF
              ELSE
                 logl(t,i) = logl(t,i) - &
                      MAX(MIN(0.4,-LOG10(zlegend(zz)/zsol)),0.2)
                 logt(t,i) = logt(t,i) + 0.1 - MIN((age(t)-age91)/1.5,0.2)
              ENDIF

           !Villaume, Conroy, Johnson (2014) normalization
           ELSE IF (tpagb_norm_type.EQ.2) THEN

              twght   = MAX(0.1,10**(-1.+(age(t)-8.0)/2.5))
              wght(i) = twght * wght(i)

           ENDIF

        ENDIF

        ! modify wght of tp-agb stars
        wght(i) = wght(i)*agb

        !add extra fudge factors to the default normalization
        logl(t,i) = logl(t,i) + dell
        logt(t,i) = logt(t,i) + delt

     ENDIF

     !modify post-AGB stars
     IF (phase(t,i).EQ.6.0.AND.pagb.NE.1.0) THEN
        wght(i) = wght(i)*pagb
     ENDIF

     !modify RGB + red clump HB + AGB
     !NB: this feature does not work w/ Padova isochrones, which
     !dont have these phases explicitly defined
     IF (phase(t,i).EQ.2.OR.phase(t,i).EQ.3 &
          .OR.phase(t,i).EQ.4.OR.phase(t,i).EQ.5) THEN
        !logt(t,i) = LOG10(10**logt(t,i) - 100.)
        wght(i) = wght(i)*redgb
     ENDIF

  ENDDO

END SUBROUTINE MOD_GB
