FUNCTION READ_FSPS, file

  sfile = strpos(file[0],'spec')
  IF sfile NE -1 THEN BEGIN
     str = read_spec(file)
     RETURN,str
  ENDIF

  sfile = strpos(file[0],'mags')
  IF sfile NE -1 THEN BEGIN
     str = read_mags(file)
     RETURN,str
  ENDIF

  sfile = strpos(file[0],'indx')
  IF sfile NE -1 THEN BEGIN
     str = read_indx(file)
     RETURN,str
  ENDIF

  sfile = strpos(file[0],'cmd')
  IF sfile NE -1 THEN BEGIN
     str = read_cmd(file)
     RETURN,str
  ENDIF

  print,'READ_FSPS ERROR: unrecognizable input file'

END
