SUBROUTINE MOD_AGB(zz,t,age,delt,dell,pagb,redgb,&
     nn,logl,logt,phase,wght)

  !routine to modify TP-AGB stars and post-AGB stars. 

  USE sps_vars; USE nrtype
  IMPLICIT NONE

  INTEGER,  INTENT(in) :: t, nn,zz
  REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: logl,logt
  REAL(SP), INTENT(in), DIMENSION(nt,nm)    :: phase
  REAL(SP), INTENT(inout), DIMENSION(nm)    :: wght
  REAL(SP), INTENT(in) :: delt, dell, pagb,redgb
  REAL(SP), INTENT(in), DIMENSION(nt) :: age
  INTEGER  :: i

  DO i=1,nn

     !modify the TP-AGB stars
     IF (phase(t,i).EQ.5.0) THEN
        !renormalize the Padova isochrones
        IF (isoc_type.EQ.'pdva') THEN
           IF (age(t).GT.8.0.AND.age(t).LT.9.1) THEN
              logl(t,i) = logl(t,i) - 1.0+(age(t)-8.)/1.5
              IF (log10(zlegend(zz)/0.0190).LT.-0.25) THEN
                 logt(t,i) = logt(t,i) + 0.10
              ENDIF
           ELSE
              logl(t,i) = logl(t,i) - &
                   MAX(MIN(0.4,-log10(zlegend(zz)/0.0190)),0.2)
              logt(t,i) = logt(t,i) + 0.1-MIN((age(t)-9.1)/1.5,0.2)
           ENDIF
        ENDIF
        logl(t,i) = logl(t,i) + dell
        logt(t,i) = logt(t,i) + delt
     ENDIF

     !modify the post-AGB stars
     IF (phase(t,i).EQ.6.0.AND.pagb.NE.1.0) THEN
        wght(i) = wght(i)*pagb
     ENDIF

     !modify the RGB + red clump HB + AGB
     IF (phase(t,i).EQ.2.OR.phase(t,i).EQ.3&
          .OR.phase(t,i).EQ.4) THEN
        wght(i) = wght(i)*redgb
     ENDIF

  ENDDO

END SUBROUTINE MOD_AGB
