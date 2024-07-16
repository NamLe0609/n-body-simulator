cmake_minimum_required(VERSION 3.14)

project(n-body-simulator)

set(SOURCES 
    src/compute.cpp
)

set(WARNING_FLAGS
    -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 
    -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual
    -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5
    -Wswitch-default -Wundef -Werror -Wno-unused
)

add_executable(run ${SOURCES})
target_compile_features(run PUBLIC cxx_std_20)

target_compile_options(run PRIVATE
    $<$<CONFIG:Debug>:-g ${WARNING_FLAGS}>
    $<$<CONFIG:Release>:-O3 ${WARNING_FLAGS}>
)