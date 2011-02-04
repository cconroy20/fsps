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
; EXAMPLE:
;
; MODIFICATION HISTORY: 
;   ? - created by CFC
;
;-
;-----------------------------------------------------------------;

FUNCTION READ_SPEC1, file,lam,nl

  openr,lun,file,/get_lun

  ;burn the header
  char = '#'
  WHILE strmid(char,0,1) EQ '#' DO BEGIN
     readf,lun,char
  ENDWHILE 

  nt = FIX(char)

  str = {agegyr:0.0,logmass:0.0,loglbol:0.0,logsfr:0.0,spec:fltarr(nl),$
         lambda:fltarr(nl)}
  str  = replicate(str,nt)
  tspec = fltarr(nl)
  t = 0.
  m = 0.
  l = 0.
  s = 0.

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

FUNCTION READ_SPEC, file, miles=miles

  ct = n_elements(file)

  spsdir = getenv('SPS_HOME')
  IF spsdir EQ '' THEN BEGIN
     print,'READ_SPEC ERROR: spsdir environment '+$
           'variable not set, returning...'
     return,0
  ENDIF

  IF keyword_set(miles) THEN BEGIN
     readcol,spsdir+'/MILES/miles.lambda',lam,/silent
  ENDIF ELSE BEGIN
     readcol,spsdir+'/BaSeL3.1/basel.lambda',lam,/silent
  ENDELSE
     
  nl = n_elements(lam)

  ff = findfile(file[0],count=ctf)
  IF ctf EQ 0 THEN BEGIN
     print,'READ_SPEC ERROR: file not found: ',file
     return,0
  ENDIF

  spec = strpos(file[0],'spec')
  IF spec EQ -1 THEN BEGIN
     print,'READ_SPEC ERROR: you did not pass a .spec file, returning'
     return,0
  ENDIF

  astr = read_spec1(file[0],lam,nl)

  IF ct GT 1 THEN BEGIN
     str = replicate(astr[0],n_elements(astr),ct)
     str[*,0] = astr
     FOR i=1,ct-1 DO str[*,i] = read_spec1(file[i],lam,nl)
  ENDIF ELSE str = astr
  
  RETURN, str

END

