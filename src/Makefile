CC=g++
TARGET_NAME = ../load2db
HASH_NAME = ../bam_hash
HASH_TEST_NAME = ../test_hash
HASH2_NAME = ../bam_hash2
HASH2_TEST_NAME = ../test_hash2
FILTER_NAME = ../bam_filter
FILTER_TEST_NAME = ../test_filter
TEST_NAME = ../test
TEST_BAM_UTILS_NAME = ../test_bam_utils
TEST_THREAD_NAME = ../test_thread

TARGET_SRCS = load2db_merge.cpp
HASH_SRCS = bam_hash_main.cpp bam_hash.cpp bamrecord.cpp
HASH_TEST_SRCS = test_hash.cpp bam_hash.cpp bamrecord.cpp
HASH2_SRCS = bam_hash2_main.cpp bam_hash2.cpp bamrecord.cpp
HASH2_TEST_SRCS = test_hash2.cpp bam_hash2.cpp bamrecord.cpp
FILTER_SRCS = bam_filter_main.cpp bam_filter.cpp bamrecord.cpp
FILTER_TEST_SRCS = test_filter.cpp bam_filter.cpp bamrecord.cpp
TEST_SRCS = test_merge.cpp
TEST_BAM_UTILS_SRCS = test_bam_utils.cpp
TEST_THREAD_SRCS = threaded_sql.cpp

COMMON_SRCS = bam_utils.cpp dbrecord.cpp sql_utils.cpp string_utils.cpp sqlite_wrapper.cpp timing.cpp bam_sql_merge.cpp

COMMON_OBJS := ${COMMON_SRCS:.cpp=.o}
TARGET_OBJ = ${TARGET_SRCS:.cpp=.o}
HASH_OBJ = ${HASH_SRCS:.cpp=.o}
HASH_TEST_OBJ = ${HASH_TEST_SRCS:.cpp=.o}
HASH2_OBJ = ${HASH2_SRCS:.cpp=.o}
HASH2_TEST_OBJ = ${HASH2_TEST_SRCS:.cpp=.o}
FILTER_OBJ = ${FILTER_SRCS:.cpp=.o}
FILTER_TEST_OBJ = ${FILTER_TEST_SRCS:.cpp=.o}
TEST_OBJ = ${TEST_SRCS:.cpp=.o}
TEST_BAM_UTILS_OBJ = ${TEST_BAM_UTILS_SRCS:.cpp=.o}
TEST_THREAD_OBJ = ${TEST_THREAD_SRCS:.cpp=.o}

INCLUDE_DIRS := /usr/local/include /home/brad/src/umi/include
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += -std=c++11 -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -fno-omit-frame-pointer

LIBRARY_DIRS := /usr/local/lib
LIBS += hts sqlite3 boost_regex boost_program_options boost_system boost_filesystem boost_thread rt z
LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library, $(LIBS), -l$(library))
TEST_LDFLAGS = -lprofiler -ltcmalloc

.PHONY: target test thread hash test_hash

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

hash: CPPFLAGS += -O3
hash: COMMON_SRCS += $(HASH_SRCS)
hash: $(HASH_NAME)

test_hash: CPPFLAGS += -DDEBUG -g -O3
test_hash: COMMON_SRCS += $(HASH_TEST_SRCS)
test_hash: $(HASH_TEST_NAME)

hash2: CPPFLAGS += -O3 -g
hash2: COMMON_SRCS += $(HASH2_SRCS)
hash2: $(HASH2_NAME)

test_hash2: CPPFLAGS += -DDEBUG -g -O3
test_hash2: COMMON_SRCS += $(HASH2_TEST_SRCS)
test_hash2: $(HASH2_TEST_NAME)

filter: CPPFLAGS += -O3
filter: COMMON_SRCS += $(FILTER_SRCS)
filter: $(FILTER_NAME)

test_filter: CPPFLAGS += -DDEBUG -g -O3
test_filter: COMMON_SRCS += $(FILTER_TEST_SRCS)
test_filter: $(FILTER_TEST_NAME)

$(TARGET_NAME): $(TARGET_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(TEST_NAME): $(TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_BAM_UTILS_NAME): $(TEST_BAM_UTILS_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(TEST_THREAD_NAME): $(TEST_THREAD_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(HASH_NAME): $(HASH_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(HASH_TEST_NAME): $(HASH_TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(HASH2_NAME): $(HASH2_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(HASH2_TEST_NAME): $(HASH2_TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

$(FILTER_NAME): $(FILTER_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(FILTER_TEST_NAME): $(FILTER_TEST_OBJ) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	@- $(RM) $(TARGET_NAME) $(TEST_NAME) $(TEST_FILTER_NAME) $(HASH2_NAME) $(TEST_HASH2_NAME)
	@- $(RM) $(COMMON_OBJS) $(TARGET_OBJ) $(TEST_OBJ) $(TEST_FILTER_OBJ) $(HASH2_OBJ) $(TEST_HASH2_OBJ)
