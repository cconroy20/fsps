	SUBROUTINE ran1_s(harvest)
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
	if (nran0 == 1) nran0=270369_k4b
	mran0=ieor(mran0,ishft(mran0,5))
	mran0=ieor(mran0,ishft(mran0,-13))
	mran0=ieor(mran0,ishft(mran0,6))
	rans=ieor(nran0,rans)+mran0
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran1_s

	SUBROUTINE ran1_v(harvest)
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
	where (nran(1:n) == 1) nran(1:n)=270369_k4b
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
	ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran1_v
