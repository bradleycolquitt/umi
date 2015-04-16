CC=g++
TARGET_NAME = ../load2db
TEST_NAME = ../test
TEST_BAM_UTILS_NAME = ../test_bam_utils
TEST_THREAD_NAME = ../test_thread

TARGET_SRCS = load2db_merge.cpp
TEST_SRCS = test_merge.cpp
TEST_BAM_UTILS_SRCS = test_bam_utils.cpp
TEST_THREAD_SRCS = threaded_sql.cpp

COMMON_SRCS = bam_utils.cpp dbrecord.cpp sql_utils.cpp string_utils.cpp sqlite_wrapper.cpp timing.cpp bam_sql_merge.cpp

COMMON_OBJS := ${COMMON_SRCS:.cpp=.o}
TARGET_OBJ = ${TARGET_SRCS:.cpp=.o}
TEST_OBJ = ${TEST_SRCS:.cpp=.o}
TEST_BAM_UTILS_OBJ = ${TEST_BAM_UTILS_SRCS:.cpp=.o}
TEST_THREAD_OBJ = ${TEST_THREAD_SRCS:.cpp=.o}


INCLUDE_DIRS := /usr/local/include /home/brad/src/umi/include
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += -std=c++11 -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -pg -fno-omit-frame-pointer -DSQLITE_ENABLE_RTREE=1

LIBRARY_DIRS := /usr/local/lib
LIBS += hts sqlite3 boost_regex boost_program_options boost_system boost_filesystem boost_thread rt
LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library, $(LIBS), -l$(library))
TEST_LDFLAGS = -lprofiler -ltcmalloc

.PHONY: target test thread

target: COMMON_SRCS += $(TARGET_SRCS)
target: $(TARGET_NAME)

test: CPPFLAGS += -DDEBUG -g
test: COMMON_SRCS += $(TEST_SRCS)
test: $(TEST_NAME)

test_bam_utils: COMMON_SRCS += $(TEST_BAM_UTILS_SRCS)
test_bam_utils: $(TEST_BAM_UTILS_NAME)

thread: CPPFLAGS += -DDEBUG -g
thread: COMMON_SRCS += $(TEST_THREAD_SRCS)
thread: $(TEST_THREAD_NAME)

$(TARGET_NAME): $(TARGET_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(TEST_NAME): $(TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_BAM_UTILS_NAME): $(TEST_BAM_UTILS_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_THREAD_NAME): $(TEST_THREAD_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	@- $(RM) $(TARGET_NAME) $(TEST_NAME)
	@- $(RM) $(COMMON_OBJS) $(TARGET_OBJ) $(TEST_OBJ)
