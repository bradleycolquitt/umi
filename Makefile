CC=g++
TARGET_NAME = load2db
TEST_NAME = test_se
TEST_PE_NAME = test_pe

## for 2pass version
SERIAL_NAME = load2db_serial
TEST_SERIAL_NAME = test_serial

## for merge version
MERGE_NAME = load2db_merge
TEST_MERGE_NAME = test_merge

TARGET_SRCS = src/load2db.cpp src/bam_sql.cpp
TEST_SRCS = src/test_se.cpp src/bam_sql.cpp
TEST_PE_SRCS = src/test_pe.cpp src/bam_sql.cpp

SERIAL_SRCS = src/load2db_serial.cpp src/bam_sql_serial.cpp
TEST_SERIAL_SRCS = src/test_serial.cpp src/bam_sql_serial.cpp

MERGE_SRCS = src/load2db_merge.cpp src/bam_sql_serial2.cpp
TEST_MERGE_SRCS = src/test_merge.cpp src/bam_sql_serial2.cpp

COMMON_SRCS = src/bam_utils.cpp src/dbrecord.cpp src/sql_utils.cpp src/string_utils.cpp src/sqlite_wrapper.cpp src/timing.cpp

COMMON_OBJS := ${COMMON_SRCS:.cpp=.o}
TARGET_OBJ = ${TARGET_SRCS:.cpp=.o}
SERIAL_OBJ = ${SERIAL_SRCS:.cpp=.o}
TEST_OBJ = ${TEST_SRCS:.cpp=.o}
TEST_PE_OBJ = ${TEST_PE_SRCS:.cpp=.o}
TEST_SERIAL_OBJ = ${TEST_SERIAL_SRCS:.cpp=.o}

MERGE_OBJ = ${MERGE_SRCS:.cpp=.o}
TEST_MERGE_OBJ = ${TEST_MERGE_SRCS:.cpp=.o}

INCLUDE_DIRS := /usr/local/include ./src
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += -std=c++11 -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -pg -fno-omit-frame-pointer -DSQLITE_ENABLE_RTREE=1

LIBRARY_DIRS := /usr/local/lib
LIBS += hts sqlite3 boost_regex boost_program_options boost_system boost_filesystem rt
LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library, $(LIBS), -l$(library))
TEST_LDFLAGS = -lprofiler -ltcmalloc

.PHONY: target test test_pe serial test_serial merge test_merge

target: $(TARGET_NAME)

serial: COMMON_SRCS += src/load2db_serial.cpp src/bam_sql_serial.cpp
serial: $(SERIAL_NAME)

test: $(TEST_NAME)

test_pe: CPPFLAGS += -DDEBUG -g -O0
test_pe: $(TEST_PE_NAME)

test_serial: CPPFLAGS += -DDEBUG -g
test_serial: COMMON_SRCS += src/test_serial.cpp src/bam_sql_serial.cpp
test_serial: $(TEST_SERIAL_NAME)

merge: COMMON_SRCS += src/load2db_merge.cpp src/bam_sql_serial2.cpp
merge: $(MERGE_NAME)

test_merge: CPPFLAGS += -DDEBUG -g
test_merge: COMMON_SRCS += src/test_merge.cpp src/bam_sql_serial2.cpp
test_merge: $(TEST_MERGE_NAME)

$(TARGET_NAME): $(TARGET_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(SERIAL_NAME): $(SERIAL_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(MERGE_NAME): $(MERGE_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(TEST_NAME): $(TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_PE_NAME): $(TEST_PE_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_SERIAL_NAME): $(TEST_SERIAL_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_MERGE_NAME): $(TEST_MERGE_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

%.o: src/%.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	@- $(RM) $(TARGET_NAME) $(SERIAL_NAME) $(MERGE_NAME) $(TEST_NAME) $(TEST_PE_NAME) $(TEST_SERIAL_NAME) $(TEST_MERGE_NAME)
	@- $(RM) $(COMMON_OBJS) $(SERIAL_OBJ) $(MERGE_OBJ) $(TARGET_OBJ) $(TEST_OBJ) $(TEST_PE_OBJ) $(TEST_SERIAL_OBJ) $(TEST_MERGE_OBJ)
