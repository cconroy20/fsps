	SUBROUTINE ran2_s(harvest)
	USE nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran0,jran0,kran0,nran0,mran0,rans
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	if (lenran < 1) call ran_init(1)
	rans=iran0-kran0
	if (rans < 0) rans=rans+2147483579_k4b
	iran0=jran0
	jran0=kran0
	kran0=rans
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	rans=iand(mran0,65535)
	mran0=ishft(3533*ishft(mran0,-16)+rans,16)+ &
		3533*rans+820265819_k4b
	rans=ieor(nran0,kran0)+mran0
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran2_s

	SUBROUTINE ran2_v(harvest)
	USE nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran,jran,kran,nran,mran,ranv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	ranv(1:n)=iran(1:n)-kran(1:n)
	where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
	iran(1:n)=jran(1:n)
	jran(1:n)=kran(1:n)
	kran(1:n)=ranv(1:n)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	ranv(1:n)=iand(mran(1:n),65535)
	mran(1:n)=ishft(3533*ishft(mran(1:n),-16)+ranv(1:n),16)+ &
		3533*ranv(1:n)+820265819_k4b
	ranv(1:n)=ieor(nran(1:n),kran(1:n))+mran(1:n)
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran2_v
