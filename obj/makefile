FILES=../gmres.o
OBJFILES = gmres.o
CXXFLAGS = -g

FORTFILES = ../gmres.f90
FORTPRGFILE = ../testgmres.f90
FILECONVDIFF = ../convdiffmain.f90
CONVDIFFOBJ = convdiff.o

FORT = gfortran

gmrestest: ${FILES}
	mv ../*.o .
	${CXX} ${OBJFILES} ../testgmres.cpp -o testgmres.exe

convdifftest: ${FILES} ../convdiff.o
	mv ../*.o .
	${CXX} ${OBJFILES} convdiff.o ../convdiffmain.cpp -o convdiff.exe

gmresfortran:
	${FORT} -c -fbounds-check ${FORTFILES}
	${FORT} -fbounds-check ${OBJFILES} ${FORTPRGFILE} -o testgmresf90.exe

fortranconvdiff:
	${FORT} -c -fbounds-check ${FORTFILES} ../convdiff.f90
	${FORT} -fbounds-check ${OBJFILES} ${CONVDIFFOBJ} ${FILECONVDIFF} -o f90convdiff.exe

clean:
	rm -rf *.o *.mod *.exe *.dat
