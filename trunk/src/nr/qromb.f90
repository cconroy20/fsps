	FUNCTION qromb(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromb
	INTERFACE
		FUNCTION func(x)
		USE sps_vars; USE nrtype
                IMPLICIT NONE
		REAL, DIMENSION(:), INTENT(IN) :: x
		REAL, DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
        !REAL(SP), PARAMETER :: EPS=1.0e-6_sp
        REAL(SP), PARAMETER :: EPS=1.0e-5_sp
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_sp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
