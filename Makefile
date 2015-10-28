
CXXFLAGS=-std=c++14 -g -Wall
LDLIBS=-lmpfr -lgtest

test: test.cpp test_quadric.cpp test_intersections.cpp test_classification.cpp test_vector.cpp test_line.cpp geometry.hpp quadrics.hpp quadric_classify.hpp origin.hpp point.hpp line.hpp bsgtree.hpp accurate_intersections.hpp accurate_math.hpp vector.hpp polynomial.hpp genericfp.hpp mathfuncs.hpp mpreal.hpp timer.o Makefile
	${CXX} ${CXXFLAGS} ${LDLIBS} test.cpp test_classification.cpp test_quadric.cpp test_vector.cpp test_line.cpp test_intersections.cpp timer.o -o test

timer.o: timer.cpp timer.hpp Makefile

clean:
	rm -f test

rebuild:
	make clean
	make test
