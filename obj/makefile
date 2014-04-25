FILES=../gmres.o
OBJFILES = gmres.o
CXXFLAGS = -g

gmrestest: ${FILES}
	mv ../*.o .
	${CXX} ${OBJFILES} ../testgmres.cpp -o testgmres.exe

convdifftest: ${FILES} ../convdiff.o
	mv ../*.o .
	${CXX} ${OBJFILES} convdiff.o ../convdiffmain.cpp -o convdiff.exe

clean:
	rm -rf *.o
