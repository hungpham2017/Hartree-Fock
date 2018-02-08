!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! H2 eletronic structures under the hydrogenic orbital
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Program main
	use io
	use integral
	use hf
	Implicit NONE
	
	call readin()
		call overlap()
		call kinetic()
		call potential()
		call twoelectron()
		call colect()
		call scf()
	
	END
