	SUBROUTINE gasdev_s(harvest)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	REAL(SP) :: rsq,v1,v2
	REAL(SP), SAVE :: g
	LOGICAL, SAVE :: gaus_stored=.false.
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		do
			call ran1(v1)
			call ran1(v2)
			v1=2.0_sp*v1-1.0_sp
			v2=2.0_sp*v2-1.0_sp
			rsq=v1**2+v2**2
			if (rsq > 0.0 .and. rsq < 1.0) exit
		end do
		rsq=sqrt(-2.0_sp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
	END SUBROUTINE gasdev_s

	SUBROUTINE gasdev_v(harvest)
	USE nrtype; USE nrutil, ONLY : array_copy
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	REAL(SP), DIMENSION(size(harvest)) :: rsq,v1,v2
	REAL(SP), ALLOCATABLE, DIMENSION(:), SAVE :: g
	INTEGER(I4B) :: n,ng,nn,m
	INTEGER(I4B), SAVE :: last_allocated=0
	LOGICAL, SAVE :: gaus_stored=.false.
	LOGICAL, DIMENSION(size(harvest)) :: mask
	n=size(harvest)
	if (n /= last_allocated) then
		if (last_allocated /= 0) deallocate(g)
		allocate(g(n))
		last_allocated=n
		gaus_stored=.false.
	end if
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		ng=1
		do
			if (ng > n) exit
			call ran1(v1(ng:n))
			call ran1(v2(ng:n))
			v1(ng:n)=2.0_sp*v1(ng:n)-1.0_sp
			v2(ng:n)=2.0_sp*v2(ng:n)-1.0_sp
			rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
			mask(ng:n)=(rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
			call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
			v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
			rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
			ng=ng+nn
		end do
		rsq=sqrt(-2.0_sp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
	END SUBROUTINE gasdev_v
