Module hf
  use io
  use integral 
  Implicit None
  Save
	DOUBLE PRECISION, ALLOCATABLE :: xmax(:,:),xmaxt(:,:),P(:,:)
	DOUBLE PRECISION	:: NN, deltaP = 1.D-8
	Integer(kind=4)		:: Ncycle = 100		

	
 Contains 

!----------------------------- 
!----------------------------- 	
!Compute Hcore, X, X+
!----------------------------- 
!----------------------------- 
	SUBROUTINE colect()
	Implicit None	
	Integer(kind=4)	 :: i,j,k,LWORK, INFO
	DOUBLE PRECISION :: A(natom,natom), ssmax(natom,natom)	
	DOUBLE PRECISION :: W(natom)
	DOUBLE PRECISION, ALLOCATABLE :: WORK(:) 
	DOUBLE PRECISION :: d2,d
	
	!-------------------------	
	!CALCULATE THE NUCLEAR-NUCLEAR TERM
	!-------------------------	
		NN =0.d0
		
		Do i=1, natom-1
			Do j=i+1, natom
				d2=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
				d= dsqrt(d2)
				NN=NN+(zc*zc)/d	
			end do
		end do	
    write(*,*) NN
		
	
	!-------------------------	
	!Tranformation matrix X
	!-------------------------
		LWORK = 3*natom-1
		allocate(WORK(LWORK))
		allocate(xmax(natom,natom))
		allocate(xmaxt(natom,natom))

		A = smax  
		!Compute the eigenvalues W (columm) of Overlap matrix (now A),
		!DSYEV replace the A with its eigenvectors 
		call DSYEV('V','U',natom,A,natom,W,WORK,LWORK,INFO)
		if(INFO.ne.0) WRITE(*,*) 'DIAGONALIZATION ERROR 1'
		
		!Construct the s-1/2 matrix from W
		ssmax=0.d0 !(s-1/2 matrix) 
		Do i=1, natom
			  ssmax(i,i)=1.d0/dsqrt(W(i))
		end do	
		
		!Compute tranformation matrix X by U*s-1/2 (U = A)
		!DGEMM compute the anpha*A*B + beta*C. In this case: A = U, B=s-1/2, anpha = 1.0 and C = 0
		!C will be replaced by the "anpha*A*B + beta*C". In this case, C = xmax = A*B = U*s-1/2
		xmax=0.d0  !transformation MATRIX 
        CALL DGEMM('N','N',natom,natom,natom,1.d0,A,natom,ssmax,natom,1.d0,xmax,natom)
        
	!-------------------------	
	!TRANSPOSITION of X
	!-------------------------	
		Do i=1,natom
			do j=1,natom
			   xmaxt(i,j) = xmax(j,i)
			end do
		end do		
 
 
	End SUBROUTINE colect

!----------------------------- 
!----------------------------- 	
!Construct G = P*twoelectron matrix from P
!----------------------------- 
!----------------------------- 
	SUBROUTINE FORMG(Pc,Gc) !give P get G
	Implicit None
	Integer(kind=4)	:: i,j,k,l
	DOUBLE PRECISION :: Pc(natom,natom),Gc(natom,natom)
	
	Gc = 0.d0
	Do i=1,natom
		Do j=1,natom
			Do k=1,natom
				Do l=1,natom
					Gc(i,j)=Gc(i,j)+Pc(k,l)*(twomax(i,j,k,l)-0.5D0*twomax(i,l,k,j))
				end do
			end do
		end do		
	end do
	
	
	End SUBROUTINE FORMG

!----------------------------- 
!----------------------------- 	
!export
!----------------------------- 
!----------------------------- 
	SUBROUTINE export(M)
	Implicit None
	Integer(kind=4)	:: i,j
	DOUBLE PRECISION::	M(natom,natom)	
	
	Do i=1,natom
		Do j=1,natom
			write(*,10,ADVANCE='NO') M(i,j)
			10 FORMAT(TL10, F15.8)
		end do		
		write(*,*)
	end do
	End SUBROUTINE export
	
	
	
	
!----------------------------- 
!----------------------------- 	
!SCF procedure
!----------------------------- 
!----------------------------- 
	SUBROUTINE scf()
	Implicit None		
	Integer(kind=4)	 :: i,j,k,iter,LWORK, INFO,nocc,haha	
	DOUBLE PRECISION :: delta, KE,PE,EE,EN,ET	
	DOUBLE PRECISION :: Pold(natom,natom),G(natom,natom),F(natom,natom),FT(natom,natom),&
									 &C(natom,natom),A(natom,natom)
	DOUBLE PRECISION :: W(natom)
	DOUBLE PRECISION, ALLOCATABLE :: WORK(:) 
	
	
	!call COLECT()
	allocate(P(natom,natom))
	G=0.0d0
	P=0.0d0
	LWORK = 3*natom-1
	allocate(WORK(LWORK))	
	iter = 0

	
10 CONTINUE		
		iter = iter+1
		!Construct Fock matrix
			call FORMG(P,G)
			Do i=1,natom
				Do j=1,natom
					F(i,j)=tmax(i,j) + vmax(i,j)+ G(i,j)
				end do
			end do 
   		
		!Calculate the transformed Fock matrix F'
			FT= 0.d0
		    CALL DGEMM('N','N',natom,natom,natom,1.d0,F,natom,xmax,natom,0.d0,C,natom)
			CALL DGEMM('N','N',natom,natom,natom,1.d0,xmaxt,natom,C,natom,0.d0,FT,natom)

			
		!Diagonalize F' to obtain C' and epxilon. After DSYEV, FT --> C' matrix
			CALL DSYEV('V','U',natom,FT,natom,W,WORK,LWORK,INFO)
			if(INFO.ne.0) WRITE(*,*) 'DIAGONALIZATION ERROR 2'
			CALL DGEMM('N','N',natom,natom,natom,1.d0,xmax,natom,FT,natom,0.d0,C,natom)
		!Construct new P from C
			
			nocc = ele*natom/2
			Do i=1,natom
				Do j=1,natom
					Pold(i,j) = P(i,j)
					P(i,j)=0.d0
					Do k=1,nocc
						P(i,j)=P(i,j)+2.d0*C(i,k)*C(j,k)
					end do
				end do
			end do 		

		!Calculate energy at the nth cycle
			call FORMG(P,G)
			KE=0.d0
			PE=0.d0
			EE=0.d0
			Do i=1,natom
				Do j=1,natom
					KE=KE+0.5D0*P(i,j)*2*tmax(i,j)
					PE=PE+0.5D0*P(i,j)*2*vmax(i,j)
					EE=EE+0.5D0*P(i,j)*G(i,j)
				end do				
			end do	
			EN=KE+PE+EE
			ET=EN+NN
			write (*,*) iter, '   ', KE,PE,EE,NN,ET	
			
			
		!Calculate delta
			delta=0.d0
			Do i=1,natom
				Do j=1,natom
					delta=delta+(P(i,j)-Pold(i,j))**2
				end do				
			end do
			delta=dsqrt(delta/4.0d0)
					
		!Check convergence	
			if(iter.LE.Ncycle .AND. delta.LE.deltaP) GO TO 10

	End SUBROUTINE scf
	
!----------------------------- 
!----------------------------- 	
!RUN HARTREE-FOCK CALCULATION
!----------------------------- 
!----------------------------- 
	SUBROUTINE runhf()
	Implicit None	
		call overlap()
		call kinetic()
		call potential()
		call twoelectron()
		call colect()
		call scf()
	End SUBROUTINE runhf
	
End Module hf