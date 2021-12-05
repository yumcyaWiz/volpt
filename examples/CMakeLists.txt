add_executable(main "main.cpp")
target_link_libraries(main PRIVATE volpt)

add_executable(sakana-test "sakana-test.cpp")
target_link_libraries(sakana-test PRIVATE volpt)