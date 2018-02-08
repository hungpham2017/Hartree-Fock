Module integral
  use io  
  Implicit None
  Save

!	Character(10) 	:: basis, state, element
!	Real(kind=8), DIMENSION(2,6,2)		:: AO
!	Real(kind=8), DIMENSION(10)		:: x,y,z	
!	Integer(kind=4)	:: natom	  
	DOUBLE PRECISION, ALLOCATABLE :: smax(:,:),tmax(:,:),vmax(:,:)
	DOUBLE PRECISION, ALLOCATABLE :: twomax(:,:,:,:)
	DOUBLE PRECISION	::PI = 3.1415926535898D0
  
 Contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!F0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	FUNCTION F0(ARG)
	Implicit None	
	DOUBLE PRECISION	:: F0, ARG
	
	IF (ARG.LT.1.0D-6) GO TO 10
	F0=DSQRT(PI/ARG)*DERF(DSQRT(ARG))/2.0D0
	GO TO 20
10	F0 = 1.0D0-ARG/3.0D0
20  CONTINUE
	RETURN

	End FUNCTION F0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FERROR FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	FUNCTION DERF(ARG)
	Implicit None	
	Integer(kind=8)	:: I	
	DOUBLE PRECISION	:: P,T,TN,POLY,DERF,ARG
	DOUBLE PRECISION, DIMENSION(5)		::	A
	
    DATA P/0.3275911D0/
    DATA A/0.254829592D0,-0.284496736D0,1.421413741D0,-1.453152027D0,1.061405429D0/
    T=1.0D0/(1.0D0+P*ARG)
    TN=T
    POLY=A(1)*TN
    DO 10 I=2,5
    TN=TN*T
    POLY=POLY+A(I)*TN
10  CONTINUE
    DERF=1.0D0-POLY*DEXP(-ARG*ARG)
    RETURN
	
	End FUNCTION DERF	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!CALCULATES OVERLAPS FOR UN-NORMALIZED PRIMITIVES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	FUNCTION S(A,B,RAB2)
	Implicit None	
	DOUBLE PRECISION	:: S,A,B,RAB2
	
    S=(PI/(A+B))**1.5D0*DEXP(-A*B*RAB2/(A+B))
    RETURN	
	
	End FUNCTION S
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!CALCULATES KINETIC ENERGY INTEGRALS FOR UN-NORMALIZED PRIMITIVES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	FUNCTION T(A,B,RAB2)
	Implicit None	
	DOUBLE PRECISION	:: T,A,B,RAB2
	
    T=A*B/(A+B)*(3.0D0-2.0D0*A*B*RAB2/(A+B))*(PI/(A+B))**1.5D0*DEXP(-A*B*RAB2/(A+B))
	RETURN	
	
	End FUNCTION T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!CALCULATES UN-NORMALIZED NUCLEAR ATTRACTION INTEGRALS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	FUNCTION V(A,B,RAB2,RCP2,ZC)
	Implicit None	
	DOUBLE PRECISION	:: V,A,B,RAB2,RCP2,ZC
	
    V=2.0D0*PI/(A+B)*F0((A+B)*RCP2)*DEXP(-A*B*RAB2/(A+B))
    V=-V*ZC
    RETURN	
	
	End FUNCTION V
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!CALCULATES TWO-ELECTRON INTEGRALS FOR UN-NORMALIZED PRIMITIVES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	FUNCTION TWOE(A,B,C,D,RAB2,RCD2,RPQ2)
	Implicit None	
	DOUBLE PRECISION	:: TWOE,A,B,C,D,RAB2,RCD2,RPQ2
	
    TWOE=2.0D0*(PI**2.5D0)/((A+B)*(C+D)*DSQRT(A+B+C+D))&
    &*F0((A+B)*(C+D)*RPQ2/(A+B+C+D))*DEXP(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D))
    RETURN	
	
	End FUNCTION TWOE
		
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!OVERLAP MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	SUBROUTINE overlap()
	Implicit None
	Integer(kind=4)	:: i,j,k,l	
	DOUBLE PRECISION 	::	A1,A2,D1,D2,R2	

	allocate(smax(natom,natom))
	smax=0.0D0
 
    Do k=1,6
		Do l=1,6
			Do i=1,natom
			  Do j=i,natom
				 A1=AO(ele,k,1)
				 A2=AO(ele,l,1)
				 D1=AO(ele,k,2)*((2.0D0*A1/PI)**0.75D0)
				 D2=AO(ele,l,2)*((2.0D0*A2/PI)**0.75D0)		
				 R2= (x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
				 smax(i,j)= smax(i,j)+S(A1,A2,R2)*D1*D2
			  end do
			end do	
		end do
	end do
	
	
	Do i=1,natom
		  Do j=i,natom
			 smax(j,i)=smax(i,j)
		  end do
	end do	
	
	End SUBROUTINE overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!KINETIC MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	SUBROUTINE kinetic()
	Implicit None
	Integer(kind=4)	:: i,j,k,l	
	DOUBLE PRECISION 	::	A1,A2,D1,D2,R2	
	
	allocate(tmax(natom,natom))
	tmax=0.0D0
 
    Do k=1,6
		Do l=1,6
			Do i=1,natom
			  Do j=i,natom
				 A1=AO(ele,k,1)
				 A2=AO(ele,l,1)
				 D1=AO(ele,k,2)*((2.0D0*A1/PI)**0.75D0)
				 D2=AO(ele,l,2)*((2.0D0*A2/PI)**0.75D0)		
				 R2= (x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
				 tmax(i,j)= tmax(i,j)+T(A1,A2,R2)*D1*D2
			  end do
			end do	
		end do
	end do
	
	
	Do i=1,natom
		  Do j=i,natom
			 tmax(j,i)=tmax(i,j)
		  end do
	end do	
	
	End SUBROUTINE kinetic
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!POTENTIAL MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	SUBROUTINE potential()
	Implicit None
	Integer(kind=4)	:: i,j,k,l,m	
	DOUBLE PRECISION 	::	A1,A2,D1,D2,xP,yP,zP,RAB2,RCP2
	
	allocate(vmax(natom,natom))
	vmax=0.0D0
 
 
		Do k=1,6
			Do l=1,6
			    Do m=1,natom
					Do i=1,natom
						Do j=i,natom
                       					  
							 A1=AO(ele,k,1)
							 A2=AO(ele,l,1)
							 D1=AO(ele,k,2)*((2.0D0*A1/PI)**0.75D0)
							 D2=AO(ele,l,2)*((2.0D0*A2/PI)**0.75D0)		
							 RAB2= (x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2 
							 xP = (A1*x(i)+A2*x(j))/(A1+A2)
							 yP = (A1*y(i)+A2*y(j))/(A1+A2)
							 zP = (A1*z(i)+A2*z(j))/(A1+A2)
							 RCP2 = (xP-x(m))**2+(yP-y(m))**2+(zP-z(m))**2
							 vmax(i,j)= vmax(i,j)+V(A1,A2,RAB2,RCP2,zc)*D1*D2
						end do
				    end do	
			    end do
		    end do
	    end do
	
	Do i=1,natom
		  Do j=i,natom
			 vmax(j,i)=vmax(i,j)
		  end do
	end do	
	
	End SUBROUTINE potential	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!TWO-ELECTRON MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
	SUBROUTINE twoelectron()
	Implicit None
	Integer(kind=4)	:: i,j,m,n,k,l,v,u	
	DOUBLE PRECISION 	::	A1,A2,A3,A4,D1,D2,D3,D4,RAB2,RCD2,RPQ2,xP,yP,zP,xQ,yQ,zQ	

	allocate(twomax(natom,natom,natom,natom))
	twomax=0.0D0
 
    Do k=1,6
		Do l=1,6
          	Do v=1,6	
				Do u=1,6
					Do i=1,natom
					  Do j=1,natom
						Do m=1,natom
							Do n=1,natom
								 A1=AO(ele,k,1)
								 A2=AO(ele,l,1)
								 A3=AO(ele,v,1)
								 A4=AO(ele,u,1)
								 D1=AO(ele,k,2)*((2.0D0*A1/PI)**0.75D0)
								 D2=AO(ele,l,2)*((2.0D0*A2/PI)**0.75D0)
								 D3=AO(ele,v,2)*((2.0D0*A3/PI)**0.75D0)
								 D4=AO(ele,u,2)*((2.0D0*A4/PI)**0.75D0)
								 RAB2 = (x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
								 RCD2 = (x(m)-x(n))**2+(y(m)-y(n))**2+(z(m)-z(n))**2
								 xP = (A1*x(i)+A2*x(j))/(A1+A2)
								 yP = (A1*y(i)+A2*y(j))/(A1+A2)
								 zP = (A1*z(i)+A2*z(j))/(A1+A2)	
								 xQ = (A3*x(m)+A4*x(n))/(A3+A4)
								 yQ = (A3*y(m)+A4*y(n))/(A3+A4)
								 zQ = (A3*z(m)+A4*z(n))/(A3+A4)									 
								 RPQ2 = (xQ-xP)**2+(yQ-yP)**2+(zQ-zP)**2
								 twomax(i,j,m,n)= twomax(i,j,m,n)+TWOE(A1,A2,A3,A4,RAB2,RCD2,RPQ2)*D1*D2*D3*D4
							end do	 
						end do		 
					  end do
					end do	
				end do
			end do		
		end do
	end do
	
	
	End SUBROUTINE twoelectron
	
End Module integral