#pragma once
#include "ophelib/packing.h"
#include "../matrixMul/matMul.h"
using namespace opaillierlib;

namespace opaillierlib
{
    class PaillierNN
    {
        public:
            struct MatrixSize
            {
                size_t height;
                size_t width;
                size_t channels;
            };

            struct Image{
                size_t height;
                size_t width;
                size_t channels;
                Vec<Vec<Vec<Integer>>> data;//data[channels][height][width]
            };

            typedef Image Filter;

            // template <typename T>
            // struct Filter{
            //     size_t height;
            //     size_t width;
            //     size_t channels;
            //     Vec<Vec<Vec< T >>> data;//data[channels][height][width]
            // };

            struct Pad{
                size_t top;
                size_t left;
                size_t bottom;
                size_t right;
            };

            struct Stride
            {
                size_t height;
                size_t width;
            };

            struct FilterSize
            {
                size_t channels;
                size_t height;
                size_t width;
                size_t filterNum;
            };

            // template <typename T>
            typedef Vec<Filter> Filters;

            enum PaddingType {VALID = 0, SAME = 1};//"VALID" = without padding; "SAME" = with zero padding,

            typedef Vec<Vec<PackedCiphertext>> IMGPackedCiphertext;

    private:
        PaillierMat paillierMat;

    public:
        PaillierNN(const size_t key_size_bits);
        PaillierMat& getPaillierMat();
        /**
         * @brief 本函数加密输入数据Image,
         * 注:在调用本函数之前,需要先用`encode`进行编码
         * 
         * @param Image 
         * @return IMGPackedCiphertext 
         */
        PaillierNN::IMGPackedCiphertext encrypt(Vec<Vec<Integer>> Image, size_t slotBits);

        Vec<Vec<Vec<Integer>>> decrypt(PaillierNN::IMGPackedCiphertext ciphertexts, PaillierNN::MatrixSize outMatSize);

        /**
         * @brief 本函数会根据filter, 步长(stride), padding等信息自动对image进行编码. 方便用于后面的加密
         * 
         * @param image 
         * @param FilterSize 
         * @param stride 步长:stride(width,height)
         * @param padType "VALID" = without padding; "SAME" = with zero padding
         * @param slotBits 明文槽的bits
         * @return Vec<Vec<Integer>> rows:filter.width * filter.height; cols:filter.width * filter.height
         */
        Vec<Vec<Integer>> encode(Image &image, const FilterSize filterSize, const Stride stride, PaddingType padType, const size_t slotBits);
        
        PaillierNN::IMGPackedCiphertext convolution(PaillierNN::IMGPackedCiphertext imgCiphertexts, Filters filters);
        PaillierNN::IMGPackedCiphertext convolutionFast(PaillierNN::IMGPackedCiphertext imgCiphertexts, Filters filters, size_t plaintextBits);
        PaillierNN::IMGPackedCiphertext convolution_Binary(PaillierNN::IMGPackedCiphertext imgCiphertexts, Filters filters);
        PaillierNN::IMGPackedCiphertext convolution_AbsOne(PaillierNN::IMGPackedCiphertext imgCiphertexts, Filters filters);
        

        PaillierNN::MatrixSize getOutSize_conv(PaillierNN::Image &image, const PaillierNN::FilterSize filterSize, const PaillierNN::Stride stride, PaillierNN::PaddingType padType);
    };
    
}