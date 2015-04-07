CC=g++
TARGET_NAME = load2db
TEST_NAME = test_se
TEST_PE_NAME = test_pe

TARGET_SRCS = src/load2db.cpp
TEST_SRCS = src/test_se.cpp
TEST_PE_SRCS = src/test_pe.cpp

COMMON_SRCS = $(filter-out $(TARGET_SRCS) $(TEST_SRCS) $(TEST_PE_SRCS), $(wildcard src/*.cpp))

COMMON_OBJS := ${COMMON_SRCS:.cpp=.o}
TARGET_OBJ = ${TARGET_SRCS:.cpp=.o}
TEST_OBJ = ${TEST_SRCS:.cpp=.o}
TEST_PE_OBJ = ${TEST_PE_SRCS:.cpp=.o}

INCLUDE_DIRS := /usr/local/include ./src
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += -std=c++11 -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -pg -fno-omit-frame-pointer

LIBRARY_DIRS := /usr/local/lib
LIBS += hts sqlite3 boost_regex boost_program_options rt
LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library, $(LIBS), -l$(library))
TEST_LDFLAGS = -lprofiler

.PHONY: target test test_pe

target: $(TARGET_NAME)

test: $(TEST_NAME)

test_pe: CPPFLAGS += -DDEBUG -g -O0
test_pe: $(TEST_PE_NAME)

$(TARGET_NAME): $(TARGET_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(TEST_NAME): $(TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_PE_NAME): $(TEST_PE_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

%.o: src/%.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	@- $(RM) $(TARGET_NAME) $(TEST_NAME) $(TEST_PE_NAME)
	@- $(RM) $(COMMON_OBJS) $(TARGET_OBJ) $(TEST_OBJ) $(TEST_PE_OBJ)
