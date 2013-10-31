 cmake_minimum_required(VERSION 2.8)

project(matrix_test)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/FindPackageHints.cmake")
  include("${CMAKE_CURRENT_LIST_DIR}/FindPackageHints.cmake")
endif()

find_package(Qt5Gui REQUIRED HINTS ${Qt5Tests_PREFIX_PATH})
find_package(Qt5Test REQUIRED HINTS ${Qt5Tests_PREFIX_PATH})
find_package(Qt5Widgets REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS}")
include_directories(
    ${Qt5Gui_INCLUDE_DIRS}
    ${Qt5Test_INCLUDE_DIRS}
)

add_definitions(
    ${Qt5Gui_DEFINITIONS}
    ${Qt5Test_DEFINITIONS}
)

set(main_file "test_matrix.cpp")
set(moc_file "${CMAKE_CURRENT_BINARY_DIR}/test_matrix.moc")

qt5_generate_moc("${main_file}" "${moc_file}")

ADD_EXECUTABLE( test_matrix test_matrix.cpp ../bem_matrix.cpp "${moc_file}")
target_link_libraries(test_matrix ${Qt5Widgets_LIBRARIES} ${Qt5Test_LIBRARIES} lapacke)

ADD_TEST( NAME test_matrix COMMAND test_matrix )
