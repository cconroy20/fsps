FUNCTION IMF(mass)

  !define IMFs (dn/dM)

  !the following is just a fast-and-loose way to pass things around
  !if the imf_type var is +10 then we calculate dn/dm*m
  !if the imf_type var is <10 then we calculate dn/dm
  
  USE sps_vars; USE nrtype
  IMPLICIT NONE

  REAL(SP), DIMENSION(:), INTENT(in) :: mass
  REAL(SP), DIMENSION(size(mass)) :: imf
  INTEGER :: i

  !Salpeter (1955) IMF
  IF (MOD(imf_type,10).EQ.0) THEN
     imf = mass**(-salp_ind) 
     IF (imf_type.EQ.10) imf = mass*imf
  ENDIF
  
  !Chabrier (2003) IMF
  IF (MOD(imf_type,10).EQ.1) THEN
     DO i=1,size(mass)
        IF (mass(i).LT.1) THEN
           imf(i) = exp(-(log10(mass(i))-log10(chab_mc))**2&
                /2/chab_sigma2)
        ELSE
           imf(i) = exp(-log10(chab_mc)**2/2./chab_sigma2)*&
                mass(i)**(-chab_ind)
        ENDIF
     ENDDO
     !convert from dn/dlnM to dn/dM
     imf = imf/mass
     IF (imf_type.EQ.11) imf = mass*imf
  ENDIF
  
  !Kroupa (2001) IMF
  IF (MOD(imf_type,10).EQ.2) THEN
     DO i=1,size(mass)
        IF (mass(i).GE.0.08.AND.mass(i).LT.0.5) &
             imf(i) = mass(i)**(-imf_alpha(1))
        IF (mass(i).GE.0.5.AND.mass(i).LT.1.0) &
             imf(i) = 0.5**(-imf_alpha(1)+imf_alpha(2))*&
             mass(i)**(-imf_alpha(2))
        IF (mass(i).GE.1.0) &
             imf(i) = 0.5**(-imf_alpha(1)+imf_alpha(2))*&
             mass(i)**(-imf_alpha(3))
     ENDDO
     IF (imf_type.EQ.12) imf = mass*imf
  ENDIF
  
  !van Dokkum (2008) IMF
  IF (MOD(imf_type,10).EQ.3) THEN
     DO i=1,size(mass)
        IF (mass(i).LE.vd_nc*vdmc) THEN
           imf(i) = vd_al*(0.5*vd_nc*vdmc)**(-vd_ind)*&
                exp(-(log10(mass(i))-log10(vdmc))*&
                (log10(mass(i))-log10(vdmc))/2./vd_sigma2)
        ELSE
           imf(i) = vd_ah*mass(i)**(-vd_ind)
        ENDIF
     ENDDO
     !convert from dn/dlnM to dn/dM
     imf = imf/mass
     IF (imf_type.EQ.13) imf = mass*imf
  ENDIF
  
END FUNCTION IMF
