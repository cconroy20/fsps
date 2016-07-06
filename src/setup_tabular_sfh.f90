! Tabulated SFH defined in a file called sfh.dat that must reside in the data
! directory (or a file specified via the parameter sfh filename; see
! below). The file must contain three columns. The first column is time since
! the Big Bang in Gyr, the second is the SFR in units of solar masses per year,
! the third is the absolute metallicity. An example is provided in the data
! directory. The time grid in this file can be arbitrary (so long as the units
! are correct), but it is up to the user to ensure that the tabulated sfh is
! well-sampled so that the outputs are stable. Obviously, highly oscillatory
! data require dense sampling.
!
! SFRs are clipped to a minimum value of 1e-30 to avoid divide by zero errors.


subroutine setup_tabular_sfh(pset, nzin)

  use sps_vars, only: sfh_tab, ntabsfh, ntabmax, nz, &
                      tiny_number, tiny30, PARAMS, SPS_HOME
  implicit none
  type(PARAMS), intent(in) :: pset
  integer, intent(in) :: nzin
  integer :: stat, n
  
  IF (pset%sfh.EQ.2) THEN

     if (pset%sf_start.gt.tiny_number) then 
        WRITE(*,*) 'COMPSP ERROR: Tabular sfh, but sf_start > 0'
        STOP
     endif
     
     ! Read the sfh file
     IF (TRIM(pset%sfh_filename).EQ.'') THEN
        OPEN(3,FILE=TRIM(SPS_HOME)//'/data/sfh.dat',ACTION='READ',STATUS='OLD')
     ELSE
        OPEN(3,FILE=TRIM(SPS_HOME)//'/data/'//TRIM(pset%sfh_filename),&
             ACTION='READ',STATUS='OLD')
     ENDIF
     DO n=1,ntabmax
        IF (nzin.EQ.nz) THEN
           READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n),sfh_tab(3,n)
        ELSE
           READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n)
           sfh_tab(3,n)=0.0
        ENDIF
        IF (stat.NE.0) GOTO 29
     ENDDO
     WRITE(*,*) 'COMPSP ERROR: didnt finish reading in the sfh file,'
     WRITE(*,*) '     increase ntabmax variable in sps_vars.f90 file'
     STOP
29   CONTINUE
     CLOSE(3)

     
     ntabsfh = n-1
     sfh_tab(1,1:ntabsfh) = sfh_tab(1,1:ntabsfh)*1E9 !convert to yrs        
     
     !special switch to compute only the last time output
     !in the tabulated file
     ! IF (pset%tage.EQ.-99.) imin=imax

  ELSE IF (pset%sfh.EQ.3) THEN 

     !sfh_tab array is supposed to already be filled in, check that it is
     IF (ntabsfh.EQ.0) THEN 
        WRITE(*,*) 'COMPSP ERROR: sfh=3 but sfh_tab array not initialized!'
        STOP
     ENDIF

  ENDIF

  ! clip SFR to a minimum of 1e-30
  do n=1, ntabsfh
     sfh_tab(2, n) = max(sfh_tab(2, n), tiny30)
  enddo
  
end subroutine setup_tabular_sfh
