SUBROUTINE SFHSTAT(pos,model,ssfr,age)

  !compute basic statistics given a parameterized star formation history.

  USE sps_vars; USE nrtype
  USE nrutil, ONLY : assert_eq; USE nr, ONLY : locate
  IMPLICIT NONE
  TYPE(PARAMS), INTENT(in)   :: pos
  TYPE(COMPSPOUT), INTENT(in):: model
  REAL(SP), INTENT(out)      :: ssfr,age

  age = pos%fburst * pos%tburst +&
       pos%const * pos%tage/2 +&
       (1.-pos%const-pos%fburst)*&
       (pos%tau + pos%sf_start)

  ssfr= pos%const/pos%tage +&
       (1.-pos%const-pos%fburst)*&
       exp(-(pos%tage-pos%sf_start)/pos%tau)
  ssfr= ssfr/model%mass_csp

END SUBROUTINE SFHSTAT
