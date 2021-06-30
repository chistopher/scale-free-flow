
add_subdirectory(googletest)

set(target tests)
add_executable(${target}
        global_network_data.hpp
        test_utils.cpp
        test_dinics.cpp
        test_push_relabel.cpp
        test_boykov_kolmogorov.cpp
        test_ghtree.cpp
        test_rational_weights.cpp
)

target_link_libraries(${target} gtest_main framework )

include(GoogleTest)
gtest_discover_tests(${target})
