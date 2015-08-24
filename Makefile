
CXXFLAGS=-std=c++14 -g -Wall
LDLIBS=-lmpfr

test: test.cpp quadrics.hpp point.hpp line.hpp bsgtree.hpp accurate_math.hpp vector.hpp

clean:
	rm -f test

rebuild:
	make clean
	make test
