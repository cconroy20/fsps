SUBROUTINE WRITE_ISOCHRONE(outfile,pset)

  !routine to write all isochrones and CMDs at a given metallicity

  USE sps_vars; USE sps_utils, ONLY : getmags,getspec,imf_weight
  IMPLICIT NONE

  INTEGER :: i,tt,zz
  TYPE(PARAMS), INTENT(in) :: pset
  CHARACTER(100), INTENT(in)  :: outfile
  CHARACTER(51)  :: fmt
  REAL(SP) :: dz=0.0,loggi
  REAL(SP), DIMENSION(nspec)  :: spec
  REAL(SP), DIMENSION(nm)     :: wght
  REAL(SP), DIMENSION(nbands) :: mags

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  zz = pset%zmet

  fmt = '(F7.4,1x,F8.4,1x,F14.9,1x,6(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(38:40),'(I3,1x,I4)') nbands

  OPEN(40,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.cmd',&
       STATUS='REPLACE')
  WRITE(40,*) '# age log(Z) mass logl logt logg '//&
       'phase composition log(weight) mags'
       
  DO tt=1,nt

     !compute IMF-based weights
     CALL IMF_WEIGHT(mini_isoc(zz,tt,:),wght,nmass_isoc(zz,tt))
     
     DO i=1,nmass_isoc(zz,tt)
        
        !get the spectrum
        CALL GETSPEC(pset,mact_isoc(zz,tt,i),logt_isoc(zz,tt,i),&
             10**logl_isoc(zz,tt,i),logg_isoc(zz,tt,i),&
             phase_isoc(zz,tt,i),ffco_isoc(zz,tt,i),spec)
        !calculate magnitudes
        CALL GETMAGS(dz,spec,mags)

        IF (isoc_type.EQ.'bsti') THEN
           loggi = LOG10( gsig4pi*mact_isoc(zz,tt,i)/&
                logl_isoc(zz,tt,i) ) + 4*logt_isoc(zz,tt,i)
        ELSE
           loggi = logg_isoc(zz,tt,i)
        ENDIF

        !write results to file
        WRITE(40,fmt) timestep_isoc(zz,tt),LOG10(zlegend(zz)),&
             mini_isoc(zz,tt,i),logl_isoc(zz,tt,i),logt_isoc(zz,tt,i),&
             loggi,phase_isoc(zz,tt,i), ffco_isoc(zz,tt,i),&
             LOG10(wght(i)),mags
        
     ENDDO
        
  ENDDO

  CLOSE(40)

END SUBROUTINE WRITE_ISOCHRONE
