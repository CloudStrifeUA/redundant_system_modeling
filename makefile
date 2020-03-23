CC=g++
INCLUDE= node.h node.cpp
CFLAGS = -g -Wall
default:
	$(CC) $(CFLAGS) -I $(INCLUDE) main.cpp -o main.o
clean:
	rm -f main.o