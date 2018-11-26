add_executable(yodaStruct_test
               ${PROJECT_SOURCE_DIR}/tests/main.cpp
)

include(Conan)
conan_cmake_run(CONANFILE ../conanfile.txt  # or relative build/conanfile.txt
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)

target_link_libraries(yodaStruct_test
                      CONAN_PKG::catch2
                      )
# Run unit tests
add_test(NAME yodaStruct-unit-tests COMMAND $<TARGET_FILE:yodaStruct_test>)
