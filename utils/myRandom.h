#pragma once

#include "../convolution/convolution.h"
using namespace std;
using namespace opaillierlib;

void initTime();

/****************gen Vector APIs********************/
void genVector_one(Vec<Integer> &vec, size_t size);

void genVector_zero(Vec<Integer> &vec, size_t size);

void genVector_binary(Vec<bool> &vec, size_t size);

void genVector_absOne(Vec<short> &vec, size_t size);

void genVector_signed(Vec<Integer> &vec, size_t size, size_t bits);

void genVector_unsigned(Vec<Integer> &vec, size_t size, size_t bits);

/****************gen Vector APIs********************/

void genMatrix_binary(Vec<Vec<Integer>> &mat, size_t height, size_t width);

void genMatrix_absOne(Vec<Vec<Integer>> &mat, size_t height, size_t width);

void genMatrix_signed(Vec<Vec<Integer>> &mat, size_t height, size_t width, size_t bits);

void genMatrix_unsigned(Vec<Vec<Integer>> &mat, size_t height, size_t width, size_t bits);

/****************gen Image APIs********************/

void genImg_binary(PaillierNN::Image &image, size_t channels, size_t height, size_t width);

void genImg_absOne(PaillierNN::Image &image, size_t channels, size_t height, size_t width);

void genImg_unsigned(PaillierNN::Image &image, size_t channels, size_t height, size_t width, size_t bits);

void genImg_signed(PaillierNN::Image &image, size_t channels, size_t height, size_t width, size_t bits);

/****************gen Filters APIs********************/
void genFilters_binary(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width);

void genFilters_absOne(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width);

void genFilters_unsigned(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width, size_t bits);

void genFilters_signed(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width, size_t bits);