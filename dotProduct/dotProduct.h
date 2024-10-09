#pragma once
#include "ophelib/packing.h"
#include "ophelib/paillier_fast.h"
#include "../paillier/newPaillier.h"

using namespace ophelib;

namespace opaillierlib
{
    class PaillierDP
    {
    private:
        // ophelib::PaillierFast paillier;
        NewPaillier paillier;
    public:
        /**
         * @brief Get the pack Count of each ciphertext 
         * 
         * @param plaintext_bits 
         * @return size_t 
         */
        size_t get_packCount(const size_t plaintext_bits);

        /**
         * Initialize a new PaillierFast instance.
         * @param key_size_bits possible values are 1024, 2048, 3072, 4096, 7680.
         *        See PaillierFast detailed description for more
         *        information. **1024 bits should not be used in
         *        production!**
         */
        PaillierDP(const size_t key_size_bits);
        
        NewPaillier& getPaillier();

        // /*------------------ 通用 Paillier APIs -----------------------*/
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

        /*----------------------- 第一套解决方案： 非打包 --------------------------*/
        /**
         * @brief encrypt: 加密一个整型向量，该函数仅进行编码，不使用打包技术
         * 编码明文，前一半单条密文长度slots为0，后一半个slots为真正的明文数据
         * 例如单条密文长度slots为6，plaintexts为{2,3,4},那么编码后{0,0,0,2,3,4}
         * 注1：本函数对应的解密函数为`decrypt`
         * 注2：由于ophelib底层加密函数预留1bits作为buffer的原因，最多打包 (keySizeBit / plaintext_bits) / 2个数据
         * 
         * @param plaintexts 待编码和加密的明文数据
         * @param slot_bits 每个slot的空间大小，因此，本参数是控制slot空间大小唯一途径
         * @return PackedCiphertext 
         */
        PackedCiphertext encrypt(const Vec<Integer> &plaintexts, size_t slot_bits);

        /**
         * @brief decrypt: 解密做完点积后的结果，计算结果保留在密文的第`slots/2`个槽，因此本函数返回下标为`slots/2`位置的值。
         * 注：本函数只针对与解密点积后的计算结果。若需要获取全部解密向量，请使用`decrypt_pack`
         * 
         * @param ciphertext 
         * @return Integer 
         */
        Integer decrypt(const PackedCiphertext &ciphertext, size_t plaintext_length);
        
        /**
         * @brief dotProduct 第一套点积方案，该套方案配套`encrypt`和`decrypt`使用。
         * 本方案主要针对只有一组向量下的点积(即：适用于仅有1组`向量*向量`的点积操作)，需要配合密文`移位`操作实现点积。
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        PackedCiphertext dotProduct(const PackedCiphertext &ciphertext, Vec<Integer> plaintexts);

        /**
         * @brief dotProduct 第一套点积方案，该套方案配套`encrypt`和`decrypt`使用。
         * 本函数主要用于plaintexts值为0/1的点积优化实现
         * 本方案主要针对只有一组向量下的点积(即：适用于仅有1组`向量*向量`的点积操作)，需要配合密文`移位`操作实现点积。
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        PackedCiphertext dotProductBinary(const PackedCiphertext &ciphertext, Vec<bool> plaintexts);

        /**
         * @brief dotProduct 第一套点积方案，该套方案配套`encrypt`和`decrypt`使用。
         * 本函数主要用于plaintexts绝对值为1(即值为-1/1)的点积优化实现
         * 本方案主要针对只有一组向量下的点积(即：适用于仅有1组`向量*向量`的点积操作)，需要配合密文`移位`操作实现点积。
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        PackedCiphertext dotProductAbsOne(const PackedCiphertext &ciphertext, Vec<short> plaintexts);

        /*----------------------- 第二套解决方案： 打包（此方案性能会更优） --------------------------*/
        
        /**
         * @brief encrypt: 加密一个整型矩阵，该函数使用打包技术
         * 假设矩阵plaintexts为：{{1,2,3},{2,5,3},{5,6,8}},那么打包后变成{{1,2,5},{2,5,6},{3,3,8}}
         * 加密后，{1,2,5}将被加密为1条密文，{2,5,6}为1条密文，以此类推。注意：本打包方法不需要填充0.
         * 注：本函数对应的解密函数为`decrypt_pack`
         * 注2：由于ophelib底层加密函数预留1bits作为buffer的原因，最多打包 (keySizeBit / plaintext_bits - 1)个数据
         * 
         * @param plaintexts 待编码和加密的明文数据
         * @param slot_bits 每个slot的空间大小，因此，本参数是控制slot空间大小唯一途径
         * @return PackedCiphertext 
         */
        Vec<PackedCiphertext> encrypt_pack(const Vec< Vec<Integer> > &plaintexts, size_t slot_bits);
        PackedCiphertext encrypt_pack(const Integer *plaintexts_begin, const Integer *plaintexts_end, const size_t slot_bits);

        /**
         * @brief decrypt_pack: 解密打包下的点积结果，结果中的每个slot都是对应一个点积结果
         * 
         * @param ciphertexts 
         * @return Vec<Integer> 
         */
        Vec<Integer> decrypt_pack(const PackedCiphertext &ciphertext);

        /**
         * @brief dotProduct_pack paillier打包下的第二套方案点积操作，配合`encrypt_pack`和`decrypt_pack`使用。
         * 本方案的性能会优于第一套方案。适用于`矩阵*向量`的点积操作
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        PackedCiphertext dotProduct_pack(const Vec<PackedCiphertext> &ciphertexts, Vec<Integer> plaintexts);
        /**
         * 
         * @brief dotProductBinary_pack paillier打包下的第二套方案点积操作，配合`encrypt_pack`和`decrypt_pack`使用。
         * 本函数主要用于plaintexts值为0/1的点积优化实现
         * 本方案的性能会优于第一套方案。适用于`矩阵*向量`的点积操作
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        PackedCiphertext dotProductBinary_pack(const Vec<PackedCiphertext> &ciphertexts, Vec<bool> plaintexts);
        /**
         * @brief dotProductAbsOne_pack paillier打包下的第二套方案点积操作，配合`encrypt_pack`和`decrypt_pack`使用。
         * 本函数主要用于plaintexts绝对值为1(即值为-1/1)的点积优化实现
         * 本方案的性能会优于第一套方案。适用于`矩阵*向量`的点积操作
         * 注：本函数不考虑结果益处情况，因此，请编码时，请考虑并分配好每个slot所需的空间大小。
         * 
         * @param ciphertext 
         * @param plaintexts 
         * @return PackedCiphertext 
         */
        PackedCiphertext dotProductAbsOne_pack(const Vec<PackedCiphertext> &ciphertexts, Vec<short> plaintexts);
    };
    
    
} // namespace opaillierlib
