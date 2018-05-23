;+
; NAME:
;       READ_SPEC 
;
; PURPOSE:
;      Function to convert a *.spec file into an IDL structure.  
;      Input may be an array of filenames.
;
; CALLING SEQUENCE:
;   res = read_spec(file)
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY: 
;   ? - created by CFC
;   09/09/11 - Updated to read in files that have arbitrary wavelength
;              arrays and where the first full line is the wavelength array
;
;-
;-----------------------------------------------------------------;

FUNCTION READ_SPEC1, file

  openr,lun,file,/get_lun

  ;burn the header
  char = '#'
  WHILE strmid(char,0,1) EQ '#' DO BEGIN
     readf,lun,char
  ENDWHILE 

  ;check if the spec file is of the "new" type, where both 
  ;the number of age steps and the number of spectral elements are included
  char = strsplit(char,' ',/extr)
  IF n_elements(char) GT 1 THEN BEGIN
     nt = long(char[0])
     nl = long(char[1]) 
  ENDIF ELSE BEGIN
     print,'ERROR: you are not passing a properly formatted *spec file'
     STOP
  ENDELSE

  str  = {logz:0.0,agegyr:0.0,logmass:0.0,loglbol:0.0,$
          logsfr:0.0,spec:dblarr(nl),lambda:fltarr(nl)}
  str   = replicate(str,nt)
  tspec = dblarr(nl)
  t = 0.
  m = 0.
  l = 0.
  s = 0.

  ;if the number of spectral elements is passed, then the first
  ;line here is the wavelength array.
  IF n_elements(char) GT 1 THEN BEGIN
     readf,lun,tspec
     lam = tspec
  ENDIF

  FOR i=0,nt-1 DO BEGIN
     
     readf,lun,t,m,l,s
     str[i].agegyr  = 10^t/1E9
     str[i].logmass = m
     str[i].loglbol = l
     str[i].logsfr  = s

     readf,lun,tspec
     str[i].spec   = tspec
     str[i].lambda = lam
     
  ENDFOR
  
  close,lun
  free_lun,lun

  RETURN,str

END

;------------------------------------------------------------;
;------------------------------------------------------------;

FUNCTION READ_SPEC, file,nohead=nohead

  spsdir = getenv('SPS_HOME')
  IF spsdir EQ '' THEN BEGIN
     print,'READ_SPEC ERROR: spsdir environment '+$
           'variable not set, returning...'
     return,0
  ENDIF

  IF NOT(keyword_set(nohead)) THEN $
     infile = spsdir+'/OUTPUTS/'+file ELSE infile=file

  file = findfile(infile,count=ct)

  IF ct EQ 0 THEN BEGIN
     print,'READ_SPEC ERROR: file not found: ',file
     return,0
  ENDIF

  spec = strpos(file[0],'spec')
  IF spec EQ -1 THEN BEGIN
     print,'READ_SPEC ERROR: you did not pass a .spec file, returning'
     return,0
  ENDIF

  astr = read_spec1(file[0])

  IF ct GT 1 THEN BEGIN
     str = replicate(astr[0],n_elements(astr),ct)
     str[*,0] = astr
     FOR i=1,ct-1 DO str[*,i] = read_spec1(file[i])
  ENDIF ELSE str = astr
  
  RETURN, str

END

