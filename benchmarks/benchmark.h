#pragma once
#include "../matrixMul/matMul.h"

using namespace opaillierlib;


/*benchmark需要的头文件和自定义*/
#include <chrono>
#define TIMER_START(name) auto start_##name = std::chrono::steady_clock::now();
#define TIMER_STOP(name) \
    auto stop_##name = std::chrono::steady_clock::now(); \
    std::chrono::duration<double, std::milli> time_##name = stop_##name - start_##name; \
    printf(#name " time: %f ms\n", time_##name.count());
    
#define TIMER_INIT(name,total) auto start_##name = std::chrono::steady_clock::now();\
                                auto total_##name = std::chrono::steady_clock::now() - std::chrono::steady_clock::now();\
                                auto pause_##name = std::chrono::steady_clock::now(); \
                                auto times_##name = total;

#define TIMER_CONTINUE(name) start_##name = std::chrono::steady_clock::now();

#define TIMER_PAUSE(name) \
    pause_##name = std::chrono::steady_clock::now(); \
    total_##name += pause_##name - start_##name;

#define TIMER_PRINT(name) \
    std::chrono::duration<double, std::milli> time_##name = total_##name; \
    printf(#name " time: %f ms\n", time_##name.count() / times_##name); 

#define TIMER_TOTAL(name1,name2,name3) \
    printf("total time: %f ms\n", (time_##name1.count() + time_##name2.count() + time_##name3.count())/ times_##name1); 


#define TIMER_TOTAL_CC(name1,name2,name3,name4) \
    printf("total time: %f ms\n", (time_##name1.count() + time_##name2.count() + time_##name3.count() + time_##name4.count())/ times_##name1); 

void benchmark_PaillierVM(size_t times, size_t key_size, size_t vecSize, size_t cipherBits, size_t plaintBits);

void benchamrk_VM(size_t times, size_t key_size, size_t vecSize, size_t cipherBits, size_t plaintBits);


void benchmark_SMP(size_t times, size_t key_size, size_t mat_n1, size_t mat_n2, size_t mat_n3, size_t cipherBits, size_t plaintBits);


void benchmark_conv(size_t times, size_t key_size,
                     size_t image_channels, size_t image_height, size_t image_width,
                     size_t filter_num, size_t filter_channels, size_t filter_height, size_t filter_width,
                     size_t stride_height, size_t stride_width, 
                     size_t cipherBits, size_t plaintBits);