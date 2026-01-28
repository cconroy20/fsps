PROGRAM SPEC_BIN

  !routine to convert ascii spectral files to binary
  !must be run twice for each value of isoc_type var

  USE sps_vars
  IMPLICIT NONE
  INTEGER  :: z,dumi1,i,j,status
  REAL(SP) :: dumr1,d2,d3
  CHARACTER(6) :: zstype
  CHARACTER(100) :: arg_spec_type

  !----------------------------------------------------------------!

  CALL GETENV('SPS_HOME',SPS_HOME)
  
  ! Read spec_type from command line or default to 'miles'
  CALL GET_COMMAND_ARGUMENT(1, arg_spec_type, STATUS=status)
  IF (status /= 0 .OR. LEN_TRIM(arg_spec_type) == 0) THEN
     spec_type = 'miles'
     WRITE(*,*) 'No spectral library specified, defaulting to: ', spec_type
  ELSE
     spec_type = TRIM(arg_spec_type)
     WRITE(*,*) 'Using spectral library: ', spec_type
  END IF

  ! Set dimensions based on spec_type
  IF (spec_type.EQ.'miles') THEN
     nzinit=5
     nspec=5994
  ELSE IF (spec_type.EQ.'basel') THEN
     nzinit=6
     nspec=1963
  ELSE IF (spec_type(1:5).EQ.'ckc14') THEN
     nzinit=11
     nspec=11149
  ELSE IF (spec_type(1:3).EQ.'c3k') THEN
     nzinit=11
     nspec=11149
  ELSE
     WRITE(*,*) 'SPEC_BIN ERROR: Unknown spec_type: ', spec_type
     STOP
  END IF

  ! Allocate arrays
  ALLOCATE(zlegendinit(nzinit))
  ALLOCATE(speclib(nspec,nzinit,ndim_logt,ndim_logg))
  
  IF (spec_type.EQ.'basel') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/zlegend.dat',&
          STATUS='OLD',ACTION='READ')
  ELSE IF (spec_type.EQ.'miles') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/zlegend.dat',&
          STATUS='OLD',ACTION='READ')
  ELSE IF (spec_type.EQ.'ckc14'.OR.spec_type(1:5).EQ.'ckc14') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/SPECTRA/CKC14/zlegend.dat',&
          STATUS='OLD',ACTION='READ')
  ELSE IF (spec_type(1:3).EQ.'c3k') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/zlegend.dat',&
          STATUS='OLD',ACTION='READ')
  ENDIF
  DO z=1,nzinit
     READ(90,'(F6.4)') zlegendinit(z)
  ENDDO
  CLOSE(90)

  DO z=1,nzinit
     
     WRITE(zstype,'(F6.4)') zlegendinit(z)

     IF (spec_type.EQ.'basel') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_'&
             //basel_str//'_z'//zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_'&
             //basel_str//'_z'//zstype//'.spectra.bin',&
             FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)

     ELSE IF (spec_type.EQ.'miles') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/imiles_z'&
             //zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/imiles_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)

     ELSE IF (spec_type(1:5).EQ.'ckc14') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/CKC14/'//spec_type//'_z'&
             //zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/CKC14/'//spec_type//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
             
     ELSE IF (spec_type(1:3).EQ.'c3k') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/'//spec_type//'_z'&
             //zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/'//spec_type//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
  
    ENDIF

     DO i=1,ndim_logg
        DO j=1,ndim_logt
           READ(92,*) dumi1,d2,d3,dumr1
           READ(92,*) speclib(:,z,j,i)
        ENDDO
     ENDDO

     WRITE(93,rec=1) speclib(:,z,:,:)
     CLOSE(92)
     CLOSE(93)
     
  ENDDO
  
END PROGRAM SPEC_BIN
