CC=g++
CXXFLAGS=-I/usr/local/include -I./src

NAME := test
CXX_SRCS := $(wildcard src/*.cpp)
CXX_OBJS := ${CXX_SRCS:.cpp=.o}
OBJS :=  $(CXX_OBJS)
INCLUDE_DIRS := /usr/local/include
LIBRARY_DIRS= /usr/local/lib .
LIBS=hts sqlite3 boost_regex
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library, $(LIBS), -l$(library))

.PHONY: all compile

all: $(NAME)

$(NAME): $(OBJS)
	$(CC)  -o $@ $^ $(LDFLAGS)

%.o: src/%.cpp
	$(CC) $(CXXFLAGS) -c -o $@ $<

clean:
	@- $(RM) $(NAME)
	@- $(RM) $(OBJS)
#bam.o: src/bam_utils.cpp src/bam_sql.cpp src/test.cpp
#	$(CC) $(CXXFLAGS) -c src/bam_utils.cpp src/bam_sql.cpp src/test.cpp
