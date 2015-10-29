
CXXFLAGS=-std=c++14 -O3 -floop-interchange -floop-strip-mine -floop-block -Wall
LDLIBS=-lmpfr
TESTLIBS=-lgtest

HEADERS=geometry.hpp quadrics.hpp quadric_classify.hpp origin.hpp point.hpp line.hpp bsgtree.hpp accurate_intersections.hpp accurate_math.hpp vector.hpp polynomial.hpp genericfp.hpp mathfuncs.hpp mpreal.hpp

all: test timing_intersections

test: test.cpp test_quadric.cpp test_intersections.cpp test_classification.cpp test_vector.cpp test_line.cpp ${HEADERS} Makefile
	${CXX} ${CXXFLAGS} ${LDLIBS} ${TESTLIBS} test.cpp test_classification.cpp test_quadric.cpp test_vector.cpp test_line.cpp test_intersections.cpp -o test

timing_intersections: timing_intersections.cpp ${HEADERS} timer.o Makefile
	${CXX} ${CXXFLAGS} ${LDLIBS} timing_intersections.cpp timer.o -o timing_intersections

timer.o: timer.cpp timer.hpp Makefile

clean:
	rm -f test timer.o timing_intersections

rebuild:
	make clean
	make test
