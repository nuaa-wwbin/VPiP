
#include "dotProduct.h"
using namespace opaillierlib;
using namespace ophelib;

PaillierDP::PaillierDP(const size_t key_size_bits_) : 
              paillier(key_size_bits_){
    paillier.generate_keys();
}

size_t PaillierDP::get_packCount(const size_t plaintext_bits){
    return paillier.get_packCount(plaintext_bits);
}

NewPaillier& PaillierDP::getPaillier(){
    return this->paillier;
}

bool PaillierDP::checkSlotSize(const size_t slot_bits, const size_t ciphertext_bits, const size_t plaintext_bits,const size_t plaintext_num){
    size_t maxBits = ciphertext_bits + plaintext_bits + log(plaintext_num);
    return slot_bits > maxBits? true : false;
}


ophelib::PackedCiphertext opaillierlib::PaillierDP::encrypt(const Vec<Integer> &plaintexts, size_t slot_bits){
    return paillier.encrypt_halfPack(plaintexts, slot_bits);
    // return paillier.encrypt_pack(plaintexts, plaintext_bits);
}


Integer opaillierlib::PaillierDP::decrypt(const PackedCiphertext &ciphertext, size_t plaintext_length) {
    Vec<Integer> ret =paillier.decrypt_pack(ciphertext);
    // std::cout << ret << std::endl;
    return ret[plaintext_length];
}

ophelib::PackedCiphertext opaillierlib::PaillierDP::dotProduct(const ophelib::PackedCiphertext &ciphertext, Vec<Integer> plaintexts){
    Integer base =1;
    size_t shift = 1;
    Ciphertext tmp_mul;
    ophelib::PackedCiphertext result = ophelib::PackedCiphertext(ciphertext.data * plaintexts[0], ciphertext.n_plaintexts, ciphertext.plaintext_bits);
    for (size_t i = 1; i < plaintexts.length(); i++)
    {
        //注:若想进一步优化,这里的`明密文乘法`可以和`移位`结合在一起,换句话说,可以减少一次大整数乘法操作.
        tmp_mul = ciphertext.data * plaintexts[i];
        tmp_mul *= (base<<(ciphertext.plaintext_bits*shift));//移位
        result.data += tmp_mul;
        shift++;
    }
    return result;
}

PackedCiphertext opaillierlib::PaillierDP::dotProductBinary(const PackedCiphertext &ciphertext, Vec<bool> plaintexts){
    Integer base =1;
    size_t shift = 1;
    Ciphertext tmp_mul;
    ophelib::PackedCiphertext result = ophelib::PackedCiphertext(ciphertext.data * plaintexts[0], ciphertext.n_plaintexts, ciphertext.plaintext_bits);
    for (size_t i = 1; i < plaintexts.length(); i++)
    {
        if (plaintexts[i])
        {
            result.data += ciphertext.data * (base<<(ciphertext.plaintext_bits*shift));//移位
        }
        shift++;
    }
    return result;
}

PackedCiphertext opaillierlib::PaillierDP::dotProductAbsOne(const PackedCiphertext &ciphertext, Vec<short> plaintexts){
    Integer base =1;
    size_t shift = 1;
    Ciphertext tmp_mul;
    ophelib::PackedCiphertext result = ophelib::PackedCiphertext(ciphertext.data * plaintexts[0], ciphertext.n_plaintexts, ciphertext.plaintext_bits);
    for (size_t i = 1; i < plaintexts.length(); i++)
    {
        if (plaintexts[i] == 1)
            result.data += ciphertext.data * (base<<(ciphertext.plaintext_bits*shift));//移位
        else if (plaintexts[i] == -1)
            result.data -= ciphertext.data * (base<<(ciphertext.plaintext_bits*shift));//移位
        else
            error_exit("plaintexts value is not 1 or -1 !\n");
        shift++;
    }
    return result;
}

/*---------------------- 第二套方案：打包 ---------------------------*/

Vec<PackedCiphertext> opaillierlib::PaillierDP::encrypt_pack(const Vec< Vec<Integer> > &plaintexts, size_t plaintext_bits){
    Vec<PackedCiphertext> result;
    result.SetLength(plaintexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        Vec<Integer> tmp;
        tmp.SetLength(plaintexts.length());
        for (size_t j = 0; j < plaintexts.length(); j++)
        {
            tmp[j] = plaintexts[j][i];
        }
        
        result[i] = paillier.encrypt_pack(tmp.begin(), tmp.end(), plaintext_bits);
    }
    return result;
}

PackedCiphertext opaillierlib::PaillierDP::encrypt_pack(const Integer *plaintexts_begin, const Integer *plaintexts_end, const size_t plaintext_bits){
    const size_t shift = plaintext_bits;
    const size_t n_plaintexts = (size_t) (plaintexts_end - plaintexts_begin);
    if(n_plaintexts > paillier.get_packCount(plaintext_bits))
        error_exit("trying to pack too many elements!\n");
    Integer sum = n_plaintexts > 0 ? *plaintexts_begin : 0;
    for (auto iter = plaintexts_begin + 1; iter < plaintexts_end; iter++) {
        if(iter->size_bits() > plaintext_bits)
            error_exit("plaintext size too large!\n");
        sum <<= shift;
        sum += *iter;
    }
    return PackedCiphertext(
        paillier.encrypt(sum),
        (const size_t) n_plaintexts,
        plaintext_bits
    );
}


Vec<Integer> opaillierlib::PaillierDP::decrypt_pack(const PackedCiphertext &ciphertext){
    return paillier.decrypt_pack(ciphertext);
}


PackedCiphertext opaillierlib::PaillierDP::dotProduct_pack(const Vec<PackedCiphertext> &ciphertexts, Vec<Integer> plaintexts){
    if (plaintexts.length() != ciphertexts.length())
        error_exit("plaintexts length is not equal to ciphertexts !");
    
    Ciphertext tmp_mul;
    PackedCiphertext result = PackedCiphertext(ciphertexts[0].data * plaintexts[0], ciphertexts[0].n_plaintexts, ciphertexts[0].plaintext_bits);
    for (size_t i = 1; i < plaintexts.length(); i++)
    {
        result.data += ciphertexts[i].data * plaintexts[i];
    }
    return result;
}

PackedCiphertext opaillierlib::PaillierDP::dotProductBinary_pack(const Vec<PackedCiphertext> &ciphertexts, Vec<bool> plaintexts){
    if (plaintexts.length() != ciphertexts.length())
        error_exit("plaintexts length is not equal to ciphertexts !");
    
    Ciphertext tmp_mul;
    PackedCiphertext result = PackedCiphertext(ciphertexts[0].data * plaintexts[0], ciphertexts[0].n_plaintexts, ciphertexts[0].plaintext_bits);
    for (size_t i = 1; i < plaintexts.length(); i++)
    {
        if (plaintexts[i])
        {
            result.data += ciphertexts[i].data;
        }
    }
    return result;
}


PackedCiphertext opaillierlib::PaillierDP::dotProductAbsOne_pack(const Vec<PackedCiphertext> &ciphertexts, Vec<short> plaintexts){
    if (plaintexts.length() != ciphertexts.length())
        error_exit("plaintexts length is not equal to ciphertexts !");
    
    Ciphertext tmp_mul;
    PackedCiphertext result = PackedCiphertext(ciphertexts[0].data * plaintexts[0], ciphertexts[0].n_plaintexts, ciphertexts[0].plaintext_bits);
    for (size_t i = 1; i < plaintexts.length(); i++)
    {
        if (plaintexts[i] == 1)
            result.data += ciphertexts[i].data;
        else if (plaintexts[i] == -1)
            result.data -= ciphertexts[i].data;
        else
            error_exit("plaintexts value is not 1 or -1 !\n");
    }
    return result;
}