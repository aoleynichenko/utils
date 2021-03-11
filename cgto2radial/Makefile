FC=gfortran

all: cgto2radial

cgto2radial:
	gfortran -c gto.f90
	gfortran -c cgto2radial.f90
	gfortran *.o -o cgto2radial.x

clean:
	rm *.mod *.o *.x




