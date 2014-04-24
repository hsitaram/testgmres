FILES=../gmres.o
OBJFILES = gmres.o
CXXFLAGS = -g

gmrestest: ${FILES}
	mv ../*.o .
	${CXX} ${OBJFILES} ../testgmres.cpp -o testgmres

clean:
	rm -rf *.o
