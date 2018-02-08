Module io
  Implicit None
  Save
	Character(10) 	:: basis, state, element
	DOUBLE PRECISION, DIMENSION(2,6,2)		:: AO
	DOUBLE PRECISION, DIMENSION(10)		:: x,y,z
	DOUBLE PRECISION		::	zc, Borh=0.52917721067d0	
	Integer(kind=4)	:: natom,ele
	
	
 Contains 
 
	SUBROUTINE readin()
	Implicit None	
	Integer(kind=4)	:: i	
	
	!-------------------------------------------------------------
	!READ method, atoms and their coordinates
	!-------------------------------------------------------------
		OPEN(UNIT=10,FILE='lucifer.in',STATUS='UNKNOWN')
		
			read (10,*) state
			state = Trim(state)
			read (10,*) element	
			element	 = Trim(element)
			if (element.EQ."H") ele=1 
			if (element.EQ."He") ele=2
			if (element.EQ."H") zc=1.0 
			if (element.EQ."He") zc=2.0

			
			read (10,*) natom		
			do i=1,natom
				read (10, *) x(i), y(i), z(i)
			end do				
		close (10)
	
			do i=1,natom
				x(i)=x(i)/0.52917721067d0
				y(i)=y(i)/0.52917721067d0
				z(i)=z(i)/0.52917721067d0
			end do	
		
	!-------------------------------------------------------------
	!READ BASIS SET
	!-------------------------------------------------------------		
		OPEN(UNIT=10,FILE='basis.in',STATUS='UNKNOWN')
		

			read (10, *)		
			do i=1,6
				read (10, *) AO(1,i,1), AO(1,i,2)
			end do	
			
			read (10, *)		
			do i=1,6
				read (10, *) AO(2,i,1), AO(2,i,2)
			end do	
		close (10)		
		
 
	End SUBROUTINE readin
	

	

End Module io