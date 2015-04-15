CC=g++
SRC_DIR = src

.PHONY: target test test_bam_utils

target:
	${MAKE} -C ${SRC_DIR} target

test:
	${MAKE} -C ${SRC_DIR} test

test_bam_utils:
	${MAKE} -C ${SRC_DIR} test_bam_utils

clean:
	@- $(RM) $(TARGET_NAME) $(SERIAL_NAME) $(MERGE_NAME) $(TEST_NAME) $(TEST_PE_NAME) $(TEST_SERIAL_NAME) $(TEST_MERGE_NAME)
	@- $(RM) $(COMMON_OBJS) $(SERIAL_OBJ) $(MERGE_OBJ) $(TARGET_OBJ) $(TEST_OBJ) $(TEST_PE_OBJ) $(TEST_SERIAL_OBJ) $(TEST_MERGE_OBJ)
