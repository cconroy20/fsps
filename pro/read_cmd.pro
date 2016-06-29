;+
; NAME:
;       READ_CMD 
;
; PURPOSE:
;      Function to convert a *.cmd file into an IDL structure.  
;      Input may be an array of filenames.
;
; CALLING SEQUENCE:
;   res = read_cmd(file)
;
; KEYWORD PARAMETERS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY: 
;   11/27/13 - created by CFC
;
;-
;-----------------------------------------------------------------;

FUNCTION READ_CMD1, file

  res = (read_ascii(file[0],data_start=1)).(0)
  nn  = n_elements(res[0,*])

  str = {logage:0.0,$
         logZ:0.0,$
         mass:0.d0,$
         loglbol:0.0,$
         logteff:0.0,$
         logg:0.0,$
         phase:0.0,$
         composition:0.0,$
         logimfweight:0.0,$
         v:0.0,$
         u:0.0,$
         b:0.0,$
         b2:0.0,$
         r:0.0,$
         i:0.0,$
         cfht_b:0.0,$
         cfht_r:0.0,$
         cfht_i:0.0,$
         twomass_j:0.0,$
         twomass_h:0.0,$
         twomass_k:0.0,$
         sdss_u:0.0,$
         sdss_g:0.0,$
         sdss_r:0.0,$
         sdss_i:0.0,$
         sdss_z:0.0,$
         wfpc2_f255w:0.0,$
         wfpc2_f300w:0.0,$
         wfpc2_f336w:0.0,$
         wfpc2_f439w:0.0,$
         wfpc2_f450w:0.0,$
         wfpc2_f555w:0.0,$
         wfpc2_f606w:0.0,$
         wfpc2_f814w:0.0,$
         wfpc2_f850lp:0.0,$
         acs_f435w:0.0,$
         acs_f475w:0.0,$
         acs_f555w:0.0,$
         acs_f606w:0.0,$
         acs_f625w:0.0,$
         acs_f775w:0.0,$
         acs_f814w:0.0,$
         acs_f850lp:0.0,$
         uvis_f218w:0.0,$
         uvis_f225w:0.0,$
         uvis_f275w:0.0,$
         uvis_f336w:0.0,$
         uvis_f390w:0.0,$
         uvis_f438w:0.0,$
         uvis_f475w:0.0,$
         uvis_f555w:0.0,$
         uvis_f606w:0.0,$
         uvis_f775w:0.0,$
         uvis_f814w:0.0,$
         uvis_f850lp:0.0,$
         wfc3_f098m:0.0,$
         wfc3_f105w:0.0,$
         wfc3_f110w:0.0,$
         wfc3_f125w:0.0,$
         wfc3_f140w:0.0,$
         wfc3_f160w:0.0,$
         irac1:0.0,$
         irac2:0.0,$
         irac3:0.0,$
         irac4:0.0,$
         isaac_k:0.0,$
         fors_v:0.0,$
         fors_r:0.0,$
         nicmos_f110w:0.0,$
         nicmos_f160w:0.0,$
         galex_fuv:0.0,$
         galex_nuv:0.0,$
         des_g:0.0,$
         des_r:0.0,$
         des_i:0.0,$
         des_z:0.0,$
         des_y:0.0,$
         wfcam_z:0.0,$
         wfcam_y:0.0,$
         wfcam_j:0.0,$
         wfcam_h:0.0,$
         wfcam_k:0.0,$
         steidel_un:0.0,$
         steidel_g:0.0,$
         steidel_rs:0.0,$
         steidel_i:0.0,$
         megacam_u:0.0,$
         megacam_g:0.0,$
         megacam_r:0.0,$
         megacam_i:0.0,$
         megacam_z:0.0,$
         wise_w1:0.0,$
         wise_w2:0.0,$
         wise_w3:0.0,$
         wise_w4:0.0,$
         uvot_w2:0.0,$
         uvot_m2:0.0,$
         uvot_w1:0.0,$
         mips_24:0.0,$
         mips_70:0.0,$
         mips_160:0.0,$
         scuba_450:0.0,$
         scuba_850:0.0,$
         pacs_70:0.0,$
         pacs_100:0.0,$
         pacs_160:0.0,$
         spire_250:0.0,$
         spire_350:0.0,$
         spire_500:0.0,$
         iras_12:0.0,$
         iras_25:0.0,$
         iras_60:0.0,$
         iras_100:0.0,$
         bessell_l:0.0,$
         bessell_l_prime:0.0,$
         bessell_m:0.0,$
         stromgren_u:0.0,$
         stromgren_v:0.0,$
         stromgren_b:0.0,$
         stromgren_y:0.0,$
         m1500:0.0,$
         m2300:0.0,$
         m2800:0.0,$
         jwst_f070w:0.0,$
         jwst_f090w:0.0,$
         jwst_f115w:0.0,$
         jwst_f150w:0.0,$
         jwst_f200w:0.0,$
         jwst_f277w:0.0,$
         jwst_f356w:0.0,$
         jwst_f444w:0.0,$
         newfirm_j1:0.0,$
         newfirm_j2:0.0,$
         newfirm_j3:0.0,$
         newfirm_h1:0.0,$
         newfirm_h2:0.0,$
         newfirm_k:0.0,$
         vista_y:0.0,$
         vista_j:0.0,$
         vista_h:0.0,$
         vista_k:0.0,$
         suprimecam_b:0.0,$
         suprimecam_g:0.0,$
         suprimecam_v:0.0,$
         suprimecam_r:0.0,$
         suprimecam_i:0.0,$
         suprimecam_z:0.0  }

  str = replicate(str,nn)

  IF n_tags(str) NE n_elements(res[*,0]) THEN BEGIN
     print,'READ_CMD ERROR: structure and *cmd file are incompatable!'
     RETURN,0
  ENDIF

  FOR i=0,n_tags(str)-1 DO str.(i) = reform(res[i,*])
 
  ;readcol,file[0],tage,tlogz,tmass,format='(F,F,D)',/sil

  RETURN,str

END

;-----------------------------------------------------------------;
;-----------------------------------------------------------------;

FUNCTION READ_CMD,file,nohead=nohead

  spsdir = getenv('SPS_HOME')
  IF spsdir EQ '' THEN BEGIN
     print,'READ_SPEC ERROR: spsdir environment '+$
           'variable not set, returning...'
     return,0
  ENDIF

  IF NOT(keyword_set(nohead)) THEN $
     file = spsdir+'/OUTPUTS/'+file

  ff = findfile(file[0],count=ct)
  IF ct EQ 0 THEN BEGIN
     print,'READ_CMD ERROR: file not found: ',file
     return,0
  ENDIF
 
  cmd = strpos(file[0],'cmd')
  IF cmd EQ -1 THEN BEGIN
     print,'READ_CMD ERROR: you did not pass a .cmd file, returning'
     return,0
  ENDIF

  str = read_cmd1(file[0])
  IF n_tags(str) EQ 0 THEN RETURN,0

  IF n_elements(file) GT 1 THEN BEGIN
     all = replicate(str[0],n_elements(str),n_elements(file))
     all[*,0] = str
     FOR i=1,n_elements(file)-1 DO all[*,i] = read_cmd1(file[i])     
  ENDIF ELSE all = str
    
  RETURN, all
  
END
