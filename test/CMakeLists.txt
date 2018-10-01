set(TEST_NAMES
    test_collision
    )

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
             WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test")
endforeach(TEST_NAME)
