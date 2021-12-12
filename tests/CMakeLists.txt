add_executable(sakana-test "sakana-test.cpp")
target_link_libraries(sakana-test PRIVATE volpt)

add_executable(test_hetero "test_hetero.cpp")
target_link_libraries(test_hetero PRIVATE volpt)