#pragma once
#include "ophelib/packing.h"
#include "../dotProduct/dotProduct.h"

using namespace opaillierlib;

namespace opaillierlib
{
    class PaillierMat
    {
    private:
        opaillierlib::PaillierDP paillierDP;
    public:
        /**
         * Initialize a new NewPaillier instance.
         * @param key_size_bits possible values are 1024, 2048, 3072, 4096, 7680.
         *        See PaillierFast detailed description for more
         *        information. **1024 bits should not be used in
         *        production!**
         */
        PaillierMat(const size_t key_size_bits);
        /**
         * @brief checkSlotSize 检测slot结算结果是否会益处,若益处将会导致解密失败(解密出没有意义的数)
         * 
         * @param slot_bits 每个slot分配的bits
         * @param ciphertext_bits 加密密文的每个数的bits
         * @param plaintext_bits 明文的每个数的bits
         * @param plaintext_num 明文数量->关系到:将会有多少个乘法和累加
         * @return true 说明slot不会存在益处,检测通过
         * @return false 说明slot可能会存在溢出,检测不通过,但是实际计算过程中,可能也能得到正确结果.因为这里只是检测slot的bits是否大于结果的最大bits
         */
        bool checkSlotSize(const size_t slot_bits, const size_t ciphertext_bits, const size_t plaintext_bits,const size_t plaintext_num);
        NewPaillier& getPaillier();

        /*------------------------- 一般矩阵乘法：明文量不超过slots ------------------------------*/
        /**
         * @brief encrypt: 加密一个整型矩阵，该函数使用打包技术
         * 假设矩阵plaintexts为：{{1,2,3},{2,5,3},{5,6,8}},那么打包后变成{{1,2,5},{2,5,6},{3,3,8}}
         * 加密后，{1,2,5}将被加密为1条密文，{2,5,6}为1条密文，以此类推。注意：本打包方法不需要填充0.
         * 注：本函数对应的解密函数为`decrypt`
         * 注2：由于ophelib底层加密函数预留1bits作为buffer的原因，最多打包 (keySizeBit / plaintext_bits - 1)个数据
         * 
         * @param plaintexts 待编码和加密的明文数据
         * @param slot_bits 每个slot的空间大小，因此，本参数是控制slot空间大小唯一途径
         * @return PackedCiphertext 
         */
        Vec<PackedCiphertext> encrypt(const Vec< Vec<Integer> > &plaintexts, size_t slot_bits);

        /**
         * @brief decrypt: 解密打包下的点积结果，结果中的每个slot都是对应一个点积结果
         * 
         * @param ciphertexts 
         * @return Vec<Vec<Integer>>
         */
        Vec<Vec<Integer>> decrypt(const Vec<PackedCiphertext> &ciphertexts);

        /**
         * @brief dotProduct_pack paillier打包下的第二套方案点积操作，配合`encrypt_pack`和`decrypt_pack`使用。
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        Vec<PackedCiphertext> matrixMul(const Vec<PackedCiphertext> &ciphertexts, Vec<Vec<Integer> > plaintexts);

        /**
         * @brief dotProduct_pack paillier打包下的矩阵乘法操作，配合`encrypt_pack`和`decrypt_pack`使用。
         * 本函数主要用于plaintexts值为0/1的点积优化实现
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        Vec<PackedCiphertext> matrixMulBinary(const Vec<PackedCiphertext> &ciphertexts, Vec<Vec<bool> > plaintexts);

        /**
         * @brief matrixMulAbsOne paillier打包下的矩阵乘法操作，配合`encrypt_pack`和`decrypt_pack`使用。
         * 本函数主要用于plaintexts绝对值为1(即值为-1/1)的点积优化实现
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        Vec<PackedCiphertext> matrixMulAbsOne(const Vec<PackedCiphertext> &ciphertexts, Vec<Vec<short> > plaintexts);

        /*------------------------- 大矩阵乘法：明文量超过slots ------------------------------*/
        /**
         * @brief encrypt: 加密一个整型矩阵，该函数使用打包技术
         * 假设矩阵plaintexts为：{{1,2,3,4},
         *                      {2,5,3,8},
         *                      {5,6,8,9},
         *                      {6,7,1,4},
         *                      {8,6,5,9}},
         * 假设每条密文只能打包2个消息,
         * 那么打包后变成{ {{1,2},{5,6},{8}},
         *              {{2,5},{6,7},{6}},
         *              {{3,3},{8,1},{5}},
         *              {{4,8},{9,4},{9}} }
         * 注1:本打包方法不需要填充0.
         * 
         * @param plaintexts 待编码和加密的明文数据
         * @param slot_bits 每个slot的空间大小，因此，本参数是控制slot空间大小唯一途径
         * @return PackedCiphertext 
         */
        Vec<Vec<PackedCiphertext>> encrypt_largeMat(const Vec< Vec<Integer> > &plaintexts, size_t slot_bits);

        /**
         * @brief encryptTrans_largeMat: 加密一个整型矩阵，该矩阵加密之前需要先进行转置.该函数使用打包技术
         * 假设矩阵plaintexts为：{{1,2,3,4},
         *                      {2,5,3,8},
         *                      {5,6,8,9},
         *                      {6,7,1,4},
         *                      {8,6,5,9}},
         * 假设每条密文只能打包2个消息,
         * 那么打包后变成{ {{1,2},{5,6},{8}},
         *              {{2,5},{6,7},{6}},
         *              {{3,3},{8,1},{5}},
         *              {{4,8},{9,4},{9}} }
         * 注1:本打包方法不需要填充0.
         * 
         * @param plaintexts 待编码和加密的明文数据
         * @param slot_bits 每个slot的空间大小，因此，本参数是控制slot空间大小唯一途径
         * @return PackedCiphertext 
         */
        Vec<Vec<PackedCiphertext>> encryptTrans_largeMat(const Vec< Vec<Integer> > &plaintexts, size_t slot_bits);

        /**
         * @brief decrypt_largeMat: 解密一个大矩阵密文，该函数使用打包技术
         * 
         * @param plaintexts 待编码和加密的明文数据
         * @param plaintext_bits 每个明文的bits大小，也是每个slot的空间大小，因此，本参数是控制slot空间大小唯一途径
         * @return PackedCiphertext 
         */
        Vec< Vec<Integer> > decrypt_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts);

        /**
         * @brief matrixMul_largeMat 配合`encrypt_largeMat`或`encryptTrans_largeMat`,以及和`decrypt_pack`配合使用。
         * 本方案的性能会优于第一套方案。适用于`矩阵*向量`的点积操作
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        Vec<Vec<PackedCiphertext>> matrixMul_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts);
        Vec<Vec<PackedCiphertext>> matrixMulFast_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts, const size_t plaintextBits);
        Vec<Vec<PackedCiphertext>> matrixMulFast_signed_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts, const size_t plaintextBits);

        /**
         * @brief matrixMulBinary_largeMat 配合`encrypt_largeMat`和`decrypt_pack`使用。
         * 本方案的性能会优于第一套方案。适用于`矩阵*向量`的点积操作
         * 本函数主要用于plaintexts值为0/1的点积优化实现
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        Vec<Vec<PackedCiphertext>> matrixMulBinary_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<bool> > plaintexts);
        Vec<Vec<PackedCiphertext>> matrixMulBinary_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, const Vec<Vec<Integer> > plaintexts);
        
        /**
         * @brief matrixMulAbsOne_largeMat 配合`encrypt_largeMat`和`decrypt_pack`使用。
         * 本方案的性能会优于第一套方案。适用于`矩阵*向量`的点积操作
         * 本函数主要用于plaintexts绝对值为1(即值为-1/1)的点积优化实现
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        Vec<Vec<PackedCiphertext>> matrixMulAbsOne_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<short> > plaintexts);
        Vec<Vec<PackedCiphertext>> matrixMulAbsOne_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts);
    // private:
        Vec<Vec<Integer>> transpose(const Vec<Vec<Integer>> &plaintextMatrix);
    };
}

