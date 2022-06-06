# Define compiler 
CC=clang++

# Define flags used by the compiler
CFLAGS=-g -Wall -Wextra -Wpedantic

# Define all files needed for this
SRC=main.cpp viewshed.cpp
HDR=viewshed.hpp
OBJ=main.o viewshed.o
BIN=viewshed

# $@ - the name of the target being generated
# $< - the first prerequisite
# $^ - list of all prerequisites


all: $(BIN) tidy

# To make the target from object files
# 	compile w/ flags -o target and first prerequisite
$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

# To make all object files fom source files
# 	compile w/ flags -c all prerequisites
%.o: %.c
	$(CC) $(CFLAGS) -c $<

tidy:
	rm -f ${OBJ}

clean: tidy
	rm -f ${BIN}
	

