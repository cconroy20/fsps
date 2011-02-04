;+
; NAME:
;       READ_MAGS 
;
; PURPOSE:
;      Function to convert a *.mags file into an IDL structure.  Input may be an 
;      array of filenames.
;
; CALLING SEQUENCE:
;   res = read_mags(file)
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

FUNCTION READ_MAGS1, file

  res = (read_ascii(file[0],data_start=8)).(0)

  nn  = n_elements(res[0,*])
  str = {agegyr:-99.,logmass:-99.,loglbol:-99.,logsfr:-99.,u:-99.,$
         b:-99.,b3:-99.,v:-99.,r:-99.,i:-99.,deep_b:-99.,deep_r:-99.,$
         deep_i:-99.,twomass_j:-99.,twomass_h:-99.,twomass_k:-99.,sdss_u:-99.,$
         sdss_g:-99.,sdss_r:-99.,sdss_i:-99.,sdss_z:-99.,wfc_f435w:-99.,$
         wfc_f606w:-99.,wfc_f775w:-99.,wfc_f814w:-99.,$
         wfc_f850lp:-99.,wfpc2_f555w:-99.,wfpc2_f814w:-99.,irac1:-99.,$
         irac2:-99.,irac3:-99.,isaac_j:-99.,$
         isaac_k:-99.,fors_v:-99.,fors_r:-99.,nicmos_f110:-99.,$
         nicmos_f160:-99.,galex_nuv:-99.,galex_fuv:-99.,wfpc2_f606w:-99.,$
         wfcam_z:-99.,wfcam_y:-99.,wfcam_j:-99.,wfcam_h:-99.,wfcam_k:-99.,$
         denis_i:-99.}
  str = replicate(str,nn)

  str.agegyr      = 10^reform(res[0,*])/1E9 ; age in Gyr
  str.logmass     = reform(res[1,*])
  str.loglbol     = reform(res[2,*])
  str.logsfr      = reform(res[3,*])
  str.v           = reform(res[4,*])
  str.u           = reform(res[5,*])
  str.deep_b      = reform(res[6,*])
  str.deep_r      = reform(res[7,*])
  str.deep_i      = reform(res[8,*])
  str.twomass_j   = reform(res[9,*])
  str.twomass_h   = reform(res[10,*])
  str.twomass_k   = reform(res[11,*])
  str.sdss_u      = reform(res[12,*])
  str.sdss_g      = reform(res[13,*])
  str.sdss_r      = reform(res[14,*])
  str.sdss_i      = reform(res[15,*])
  str.sdss_z      = reform(res[16,*])
  str.wfc_f435w   = reform(res[17,*])
  str.wfc_f606w   = reform(res[18,*])
  str.wfc_f775w   = reform(res[19,*])
  str.wfc_f814w   = reform(res[20,*])
  str.wfc_f850lp  = reform(res[21,*])
  str.irac1       = reform(res[22,*])
  str.irac2       = reform(res[23,*])
  str.nicmos_f110 = reform(res[29,*])
  str.nicmos_f160 = reform(res[30,*])
  str.fors_v      = reform(res[27,*])
  str.fors_r      = reform(res[28,*])
  str.galex_nuv   = reform(res[31,*])
  str.galex_fuv   = reform(res[32,*])
  str.wfcam_z     = reform(res[38,*])
  str.wfcam_y     = reform(res[39,*])
  str.wfcam_j     = reform(res[40,*])
  str.wfcam_h     = reform(res[41,*])
  str.wfcam_k     = reform(res[42,*])
  IF n_elements(res[*,0]) GE 46 THEN BEGIN
     str.b = reform(res[43,*])
     str.r = reform(res[44,*])
     str.i = reform(res[45,*])
     IF n_elements(res[*,0]) GE 47 THEN BEGIN
        str.b3 = reform(res[46,*])
        IF n_elements(res[*,0]) GE 48 THEN BEGIN
           str.wfpc2_f555w = reform(res[47,*])
           str.wfpc2_f814w = reform(res[48,*])
        ENDIF
     ENDIF
  ENDIF
 
  RETURN,str

END

;-----------------------------------------------------------------;


FUNCTION READ_MAGS,file

  ff = findfile(file[0],count=ct)
  IF ct EQ 0 THEN BEGIN
     print,'READ_MAGS ERROR: file not found: ',file
     return,0
  ENDIF
 
  mags = strpos(file[0],'mags')
  IF mags EQ -1 THEN BEGIN
     print,'READ_MAGS ERROR: you did not pass a .mags file, returning'
     return,0
  ENDIF

  str = read_mags1(file[0])

  IF n_elements(file) GT 1 THEN BEGIN
     all = replicate(str[0],n_elements(str),n_elements(file))
     all[*,0] = str
     FOR i=1,n_elements(file)-1 DO all[*,i] = read_mags1(file[i])     
  ENDIF ELSE all = str
    
  RETURN, all
  
END
