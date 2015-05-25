CC=g++
SRC_DIR = src

.PHONY: target test test_bam_utils thread hash test_hash filter test_filter

target:
	${MAKE} -C ${SRC_DIR} target

test:
	${MAKE} -C ${SRC_DIR} test

test_bam_utils:
	${MAKE} -C ${SRC_DIR} test_bam_utils

thread:
	${MAKE} -C ${SRC_DIR} thread

hash:
	${MAKE} -C ${SRC_DIR} hash

test_hash:
	${MAKE} -C ${SRC_DIR} test_hash

filter:
	${MAKE} -C ${SRC_DIR} filter

test_filter:
	${MAKE} -C ${SRC_DIR} test_filter

clean:
	${MAKE} -C ${SRC_DIR} clean
