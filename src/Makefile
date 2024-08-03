CC =gcc
CFLAGS =-Wall -Werror -Wextra -std=c11 -g
GCOV_FLAGS =-fprofile-arcs -ftest-coverage
FLAGS_CHECK = $(shell pkg-config --cflags --libs check)

ifeq ($(shell uname -s),Linux)
	TEST_FLAGS+=-lsubunit
endif

all: s21_matrix.a 

s21_matrix.a:
	$(CC) $(CFLAGS) -c s21_matrix.c
	ar rc s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a
	@-rm -rf *.o

test:
	$(CC) $(CFLAGS) -c test.c
	$(CC)  *.o s21_matrix.a -o test $(FLAGS_CHECK) -lcheck -lm -lsubunit -g
	./test

gcov_report:
	gcc $(GCOV_FLAGS) *.c -o report.out $(FLAGS_CHECK) -lcheck -lm -lsubunit
	./report.out
	lcov -t "report" -o report.info -c -d .
	genhtml -o ./report report.info
	open  ./report/index.html

clean:
	@-rm *.gcda *.gcno *.info *.out *.o *.a test *.log
	@rm -rf report

rebuild: clean s21_matrix.a

clang:
	cp ../materials/linters/.clang-format .
	clang-format -n *.c
	rm .clang-format

format:
	cp ../materials/linters/.clang-format .
	clang-format *.c -style=Google -i
	clang-format *.h -style=Google -i
	rm .clang-format

testing:
	@gcc -Wall -Werror -Wextra -std=c11 main.c s21_matrix.c -o test -lm
	@./test
	@rm test
 
valgrind:
	@valgrind --leak-check=full --show-leak-kinds=all -s --track-origins=yes --trace-children=yes --track-fds=yes --log-file=valgrind.log ./test
	@grep -F "ERROR SUMMARY" valgrind.log
