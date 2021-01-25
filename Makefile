FLAGS = -Wunused -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -O3


all: a.out

a.out: functions.o put_matrix.o solution.o test.o
	g++ functions.o put_matrix.o solution.o test.o -o a.out 

functions.o: functions.cpp functions.h
	g++ -c $(FLAGS) functions.cpp

put_matrix.o: put_matrix.cpp functions.h
	g++ -c $(FLAGS) put_matrix.cpp

solution.o: solution.cpp functions.h
	g++ -c $(FLAGS) solution.cpp

test.o: test.cpp functions.h
	g++ -c $(FLAGS) test.cpp

clean:
	rm -rf *.0 a.out
