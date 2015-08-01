
CXXFLAGS=-std=c++11 -g -Wall

test: test.cpp quadrics.hpp point.hpp line.hpp bsgtree.hpp accurate_math.hpp

clean:
	rm -f test
