PROGRAM CALC_MAGSUN

  !Compute magnitude of Sun, print to file

  !set up modules
  USE sps_vars; USE sps_utils  
  IMPLICIT NONE

  INTEGER :: i
  REAL(SP), DIMENSION(nbands) :: magv, magab

  compute_vega_mags=1
  CALL SPS_SETUP(1)
  magv = magsun

  compute_vega_mags=0
  CALL SPS_SETUP(1)
  magab = magsun

 
  OPEN(1,FILE=TRIM(SPS_HOME)//'data/magsun.dat',&
    ACTION='WRITE',STATUS='REPLACE')
  WRITE(1,*) '# magnitude of the Sun in the filters listed in FILTER_LIST'
  WRITE(1,*) '# columns are: index, AB mag, Vega mag'
  DO i=1,nbands
     WRITE(1,'(I3,1x,2F7.2)') i,magab(i),magv(i)
  ENDDO
  CLOSE(1)

END PROGRAM CALC_MAGSUN
