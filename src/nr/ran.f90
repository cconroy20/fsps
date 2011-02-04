	FUNCTION ran(idum)
	IMPLICIT NONE
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
	INTEGER(K4B), INTENT(INOUT) :: idum
	REAL :: ran
	INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
	REAL, SAVE :: am
	INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
	if (idum <= 0 .or. iy < 0) then
		am=nearest(1.0,-1.0)/IM
		iy=ior(ieor(888889999,abs(idum)),1)
		ix=ieor(777755555,abs(idum))
		idum=abs(idum)+1
	end if
	ix=ieor(ix,ishft(ix,13))
	ix=ieor(ix,ishft(ix,-17))
	ix=ieor(ix,ishft(ix,5))
	k=iy/IQ
	iy=IA*(iy-k*IQ)-IR*k
	if (iy < 0) iy=iy+IM
	ran=am*ior(iand(IM,ieor(ix,iy)),1)
	END FUNCTION ran
