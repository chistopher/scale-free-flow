
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)
project(flows VERSION 1.0.0)

find_package(girgs REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

set (CMAKE_CXX_STANDARD 20)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-march=core-avx2)
endif()

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-march=native)
    add_compile_options(-Wall -Wextra -Wshadow -Wfloat-equal -Wconversion -Wno-sign-compare)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
                -fsanitize=address,undefined
                -fno-omit-frame-pointer
                -fno-sanitize-recover=all
        )
        add_link_options(
                -fsanitize=address,undefined
                -fno-omit-frame-pointer)
        # add_compile_definitions(_GLIBCXX_DEBUG) # gcc only
    endif()
endif()

set(HEADERS
        include/utils.h
        include/ghtree.h
        include/Dinics.h
        include/DinicsStats.h
        include/minimalDinics.h
        include/PushRelabel.h
        include/BoykovKolmogorov.h

        include/ScopedTimer.hpp
        include/ProgressIndicator.hpp
)

set(SOURCES
        source/utils.cpp
        source/ghtree.cpp
        source/Dinics.cpp
        source/DinicsStats.cpp
        source/minimalDinics.cpp
        source/PushRelabel.cpp
        source/BoykovKolmogorov.cpp
)
configure_file(include/common.h.in include/common.h)

add_library(framework ${HEADERS} ${SOURCES})
target_include_directories(framework PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include)
target_include_directories(framework PRIVATE ${Boost_INCLUDE_DIR})

add_executable(details details.cpp)
target_link_libraries(details framework)

add_executable(pairs pairs.cpp)
target_link_libraries(pairs framework)

add_executable(runtime runtime.cpp)
target_link_libraries(runtime framework)

add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox framework girgs::girgs)

add_executable(sandbox2 sandbox2.cpp)
target_link_libraries(sandbox2 framework girgs::girgs)

add_executable(scalingGHonGIRGs scalingGHonGIRGs.cpp)
target_link_libraries(scalingGHonGIRGs framework girgs::girgs)

add_executable(show_search show_search.cpp)
target_link_libraries(show_search PRIVATE girgs::girgs)

add_executable(gh-profile gh-profile.cpp)
target_link_libraries(gh-profile framework)

add_executable(clustering clustering.cpp)
target_link_libraries(clustering framework girgs::girgs)

add_executable(other_applications other_applications.cpp)
target_link_libraries(other_applications framework)

enable_testing()
add_subdirectory(tests)
