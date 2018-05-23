;+
; NAME:
;       READ_INDX 
;
; PURPOSE:
;      Function to convert a *.indx file into an IDL structure.  
;      Input may be an array of filenames.
;
; CALLING SEQUENCE:
;   res = read_indx(file)
;
; KEYWORD PARAMETERS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY: 
;   ? - created by CFC
;   04/11 - added near-IR indices
;
;-
;-----------------------------------------------------------------;

FUNCTION READ_INDX1, file

  res = (read_ascii(file[0],data_start=8)).(0)

  ;extract the metallicity from the header
  openr,lun,file[0],/get_lun
  ss = ''
  readf,lun,ss
  close,lun
  free_lun,lun
  zmet = float((strsplit(ss,':',/extr))[1])
  
  nn  = n_elements(res[0,*])
  str = {logz:0.0,agegyr:0.0,CN1:0.0,CN2:0.0,CA4227:0.0,$
         G4300:0.0,Fe4383:0.0,$
         Ca4455:0.0,Fe4531:0.0,C4668:0.0,Hb:0.0,Fe5015:0.0,Mg1:0.0,$
         Mg2:0.0,Mgb:0.0,Fe5270:0.0,Fe5335:0.0,Fe5406:0.0,Fe5709:0.0,$ 
         Fe5782:0.0,NaD:0.0,TiO1:0.0,TiO2:0.0,HdA:0.0,HgA:0.0,HdF:0.0,$
         HgF:0.0,Dn4000:0.0,CO:0.0,H2O:0.0,cn_ir:0.0,C2:0.0,MgFe:0.0}
  str = replicate(str,nn)

  str.agegyr = reform(res[0,*])
  str.logz   = zmet
  str.cn1    = reform(res[1,*])
  str.cn2    = reform(res[2,*])
  str.ca4227 = reform(res[3,*])
  str.g4300  = reform(res[4,*])
  str.fe4383 = reform(res[5,*])
  str.ca4455 = reform(res[6,*])
  str.fe4531 = reform(res[7,*])
  str.c4668  = reform(res[8,*])
  str.hb     = reform(res[9,*])
  str.fe5015 = reform(res[10,*])
  str.mg1    = reform(res[11,*])
  str.mg2    = reform(res[12,*])
  str.mgb    = reform(res[13,*])
  str.fe5270 = reform(res[14,*])
  str.fe5335 = reform(res[15,*])
  str.fe5406 = reform(res[16,*])
  str.fe5709 = reform(res[17,*])
  str.fe5782 = reform(res[18,*])
  str.nad    = reform(res[19,*])
  str.tio1   = reform(res[20,*])
  str.tio2   = reform(res[21,*])
  str.hda    = reform(res[22,*])
  str.hga    = reform(res[23,*])
  str.hdf    = reform(res[24,*])
  str.hgf    = reform(res[25,*])
  str.dn4000 = reform(res[26,*])
  str.co     = reform(res[27,*])
  str.h2o    = reform(res[28,*])
  str.cn_ir  = reform(res[29,*])
  str.c2     = reform(res[30,*])

  ;a/fe insensitive index (see Thomas et al. 2003)
  str.mgfe = sqrt(str.mgb*(0.72*str.fe5270+0.28*str.fe5335)) 

  RETURN,str

END

;-----------------------------------------------------------------;


FUNCTION READ_INDX,file

  spsdir = getenv('SPS_HOME')
  IF spsdir EQ '' THEN BEGIN
     print,'READ_SPEC ERROR: spsdir environment '+$
           'variable not set, returning...'
     return,0
  ENDIF

  file = spsdir+'/OUTPUTS/'+file

  ff = findfile(file[0],count=ct)
  IF ct EQ 0 THEN BEGIN
     print,'READ_INDX ERROR: file not found: ',file
     return,0
  ENDIF
 
  mags = strpos(file[0],'indx')
  IF mags EQ -1 THEN BEGIN
     print,'READ_INDX ERROR: you did not pass a .indx file, returning'
     return,0
  ENDIF

  str = read_indx1(file[0])

  IF n_elements(file) GT 1 THEN BEGIN
     all = replicate(str[0],n_elements(str),n_elements(file))
     all[*,0] = str
     FOR i=1,n_elements(file)-1 DO all[*,i] = read_indx1(file[i])     
  ENDIF ELSE all = str
    
  RETURN, all
  
END
