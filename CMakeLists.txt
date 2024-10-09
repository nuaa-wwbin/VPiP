cmake_minimum_required(VERSION 3.13)

project(PAILLIER_ATOMIC VERSION 1.0 LANGUAGES CXX)

# 设置指定的C++编译器版本是必须的，如果不设置，或者为OFF，则指定版本不可用时，会使用上一版本。
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 指定为C++17 版本
set(CMAKE_CXX_STANDARD 17)

# If not called from root CMakeLists.txt
if(NOT DEFINED PAILLIER_ATOMIC_BUILD)
    set(PAILLIER_ATOMIC_BUILD ON)

    # Import OpherLib
    list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")
    find_package(OPHELib REQUIRED)
    find_package(Threads REQUIRED)
    find_package(OpenMP)
    # find_package(NTL REQUIRED)
    # find_package(GMP REQUIRED)
    # find_package(GMPXX REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()


if(PAILLIER_ATOMIC_BUILD)
    add_executable(paillier_atomic)
    # target_compile_features(ophelib PRIVATE cxx_std_11)
    target_sources(paillier_atomic
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/main.cpp
            ${CMAKE_CURRENT_LIST_DIR}/dotProduct/dotProduct.cpp
            ${CMAKE_CURRENT_LIST_DIR}/paillier/newPaillier.cpp
            ${CMAKE_CURRENT_LIST_DIR}/matrixMul/matMul.cpp
            ${CMAKE_CURRENT_LIST_DIR}/convolution/convolution.cpp
            ${CMAKE_CURRENT_LIST_DIR}/utils/myRandom.cpp
            ${CMAKE_CURRENT_LIST_DIR}/benchmarks/benchmark.cpp
    )

    target_link_libraries(paillier_atomic 
                            PRIVATE 
                            OPHELib::ophelib 
                            -lntl
                            -lm
                            -lgmp
                            Threads::Threads)
                            # ${LIBNTL_LIBRARIES}
                            # ${LIBGMP_LIBRARIES}
                            # ${LIBGMPXX_LIBRARIES}
endif()