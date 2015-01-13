CC=g++
CXXFLAGS=-I/usr/local/include -I./src
LDFLAGS=-L/usr/local/lib
LIBS=-lhts
OBJS = bam_utils.o test.o

all: test

test: $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o test $(LDFLAGS) $(LIBS)


test.o: src/test.cpp src/bam_utils.h
	$(CC) $(CXXFLAGS) -c src/test.cpp

bam_utils.o: src/bam_utils.cpp
	$(CC) $(CXXFLAGS) -c src/bam_utils.cpp
