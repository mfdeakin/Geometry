
CXXFLAGS=-std=c++14 -g -Wall
LDLIBS=-lmpfr -lgtest

test: test.cpp quadrics.hpp point.hpp line.hpp bsgtree.hpp accurate_math.hpp vector.hpp
	${CXX} ${CXXFLAGS} ${LDLIBS} test.cpp -o test

clean:
	rm -f test

rebuild:
	make clean
	make test
