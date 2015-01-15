CC=g++

NAME := load2db
CXX_SRCS := $(filter-out src/test.cpp, $(wildcard src/*.cpp))

INCLUDE_DIRS := /usr/local/include ./src
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += -std=c++11


LIBRARY_DIRS := /usr/local/lib
LIBS += hts sqlite3 boost_regex boost_program_options

ifdef TEST
CXX_SRCS := $(filter-out src/load2db.cpp, $(wildcard src/*.cpp))
NAME := test
LIBS += profiler
CPPFLAGS += -pg -fno-omit-frame-pointer
endif

LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library, $(LIBS), -l$(library))

OBJS := ${CXX_SRCS:.cpp=.o}

.PHONY: all compile

all: $(NAME)

$(NAME): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: src/%.cpp
	$(CC) $(CXXFLAGS) -c -o $@ $<

clean:
	@- $(RM) $(NAME)
	@- $(RM) $(OBJS)
