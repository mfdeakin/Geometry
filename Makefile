
CXX=nvcc
CXXFLAGS=-std=c++11 -g
LDLIBS=-lmpfr
TESTLIBS=-lgtest -lpthread

HEADERS=geometry.hpp quadrics.hpp quadric_classify.hpp origin.hpp point.hpp line.hpp bsgtree.hpp accurate_intersections.hpp accurate_math.hpp vector.hpp polynomial.hpp genericfp.hpp mathfuncs.hpp mpreal.hpp

all: test timing_intersections

test: test.cpp test_quadric.cpp test_intersections.cpp test_classification.cpp test_vector.cpp test_line.cpp ${HEADERS} Makefile
	${CXX} ${CXXFLAGS} test.cpp test_classification.cpp test_quadric.cpp test_vector.cpp test_line.cpp test_intersections.cpp -o test ${LDLIBS} ${TESTLIBS}

timing_intersections: timing_intersections.cpp ${HEADERS} timer.o Makefile
	${CXX} ${CXXFLAGS} timing_intersections.cpp timer.o -o timing_intersections ${LDLIBS}

timer.o: timer.cpp timer.hpp Makefile

clean:
	rm -f test timer.o timing_intersections

rebuild:
	make clean
	make test
