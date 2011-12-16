SUBROUTINE WRITE_ISOCHRONE(file,zz)

  !routine to write all isochrones and CMDs at a given metallicity

  USE sps_vars
  IMPLICIT NONE

  INTEGER :: i,tt
  INTEGER, INTENT(in) :: zz
  CHARACTER(100), INTENT(in)  :: file
  REAL(SP), DIMENSION(nspec)  :: spec
  REAL(SP), DIMENSION(nm)     :: wght
  REAL(SP), DIMENSION(nbands) :: mags

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  OPEN(55,FILE=TRIM(SPS_HOME)//'OUTPUTS/'//TRIM(file),STATUS='REPLACE')
  WRITE(55,*) '# log(age/yr), log(Z/Zsol), mini, logl, logt, '//&
       'logg, phase, IMF weight, mags'
      
  !loop on the age of the isochrone
  DO tt=1,nt
     
     !compute IMF-based weights
     !each point at mass M represents a *bin*: M-dM/2<M<M+dM/2
     CALL IMF_WEIGHT(mini_isoc(zz,tt,:),wght,nmass_isoc(zz,tt))
     
     !loop on each star in the isochrone
     DO i=1,nmass_isoc(zz,tt)
        
        !get the spectrum
        CALL GETSPEC(zz,mini_isoc(zz,tt,i),mact_isoc(zz,tt,i),&
             logt_isoc(zz,tt,i),10**logl_isoc(zz,tt,i),&
             phase_isoc(zz,tt,i),ffco_isoc(zz,tt,i),spec)
        !calculate magnitudes
        CALL GETMAGS(0.0,spec,mags)
           
        WRITE(55,10) timestep_isoc(zz,tt),LOG10(zlegend(zz)/0.0190),&
             mini_isoc(zz,tt,i),logl_isoc(zz,tt,i),logt_isoc(zz,tt,i),&
             logg_isoc(zz,tt,i),phase_isoc(zz,tt,i),wght(i),mags
        
     ENDDO
     
  ENDDO
  
  CLOSE(55)

10 FORMAT (7(F7.3,1x),ES10.3,1x,99(F7.3,1x))

END SUBROUTINE WRITE_ISOCHRONE
