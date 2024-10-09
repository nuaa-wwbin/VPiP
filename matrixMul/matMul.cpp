#include "matMul.h"
using namespace opaillierlib;

PaillierMat::PaillierMat(const size_t key_size_bits_) : 
              paillierDP(key_size_bits_){
}

bool PaillierMat::checkSlotSize(const size_t slot_bits, const size_t ciphertext_bits, const size_t plaintext_bits,const size_t plaintext_num){
    return paillierDP.checkSlotSize(slot_bits,ciphertext_bits,plaintext_bits,plaintext_num);
}

NewPaillier& PaillierMat::getPaillier(){
    return paillierDP.getPaillier();
}

Vec<PackedCiphertext> opaillierlib::PaillierMat::encrypt(const Vec< Vec<Integer> > &plaintexts, size_t plaintext_bits){
    return paillierDP.encrypt_pack(plaintexts, plaintext_bits);
}

Vec<PackedCiphertext> opaillierlib::PaillierMat::matrixMul(const Vec<PackedCiphertext> &ciphertexts, Vec<Vec<Integer> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    
    Vec<PackedCiphertext> result;
    result.SetLength(plaintexts[0].length());
    for (size_t i = 0; i < plaintexts[0].length(); i++)
    {
        for (size_t j = 0; j < plaintexts.length(); j++)
        {
            if (0 == j){
                result[i] = PackedCiphertext(ciphertexts[j].data * plaintexts[j][i],ciphertexts[j].n_plaintexts, ciphertexts[j].plaintext_bits);
            }
            else
                result[i].data += ciphertexts[j].data * plaintexts[j][i];
        }
    }
    return result;
}

Vec<PackedCiphertext> opaillierlib::PaillierMat::matrixMulBinary(const Vec<PackedCiphertext> &ciphertexts, Vec<Vec<bool> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    
    Vec<PackedCiphertext> result;
    result.SetLength(plaintexts[0].length());
    for (size_t i = 0; i < plaintexts[0].length(); i++)
    {
        for (size_t j = 0; j < plaintexts.length(); j++)
        {
            if (0 == j)
                result[i] = PackedCiphertext(ciphertexts[j].data * plaintexts[j][i],ciphertexts[j].n_plaintexts, ciphertexts[j].plaintext_bits);
            else if (plaintexts[j][i])
                result[i].data += ciphertexts[j].data;
        }
    }
    return result;
}
Vec<PackedCiphertext> opaillierlib::PaillierMat::matrixMulAbsOne(const Vec<PackedCiphertext> &ciphertexts, Vec<Vec<short> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    
    Vec<PackedCiphertext> result;
    result.SetLength(plaintexts[0].length());
    for (size_t i = 0; i < plaintexts[0].length(); i++)
    {
        for (size_t j = 0; j < plaintexts.length(); j++)
        {
            if (0 == j)
                result[i] = PackedCiphertext(ciphertexts[j].data * plaintexts[j][i],ciphertexts[j].n_plaintexts, ciphertexts[j].plaintext_bits);
            else if (plaintexts[j][i] == 1)
                result[i].data += ciphertexts[j].data;
            else if (plaintexts[j][i] == -1)
                result[i].data -= ciphertexts[j].data;
            else
                error_exit("plaintexts value is not 1 or -1 !\n");
        }
    }
    return result;
}

Vec<Vec<Integer>> opaillierlib::PaillierMat::decrypt(const Vec<PackedCiphertext> &ciphertexts){
    Vec< Vec<Integer> > result;
    result.SetLength(ciphertexts.length());
    for (size_t i = 0; i < ciphertexts.length(); i++)
    {
        result[i] = paillierDP.decrypt_pack(ciphertexts[i]);
    }
    return transpose(result);
}

Vec<Vec<Integer>> opaillierlib::PaillierMat::transpose(const Vec<Vec<Integer>> &plaintextMatrix) {
    size_t rows = plaintextMatrix.length();
    //校验矩阵是否正常
    if (rows<1)//无元素
        error_exit("plaintextMatrix has no elements!\n");
    size_t cols = plaintextMatrix[0].length();
    Vec<Vec<Integer>> result;
    result.SetLength(cols);
    for (size_t i = 0; i < cols; i++)
    {
        result[i].SetLength(rows);
        for (size_t j = 0; j < rows; j++)
        {
            if (i == 0 && plaintextMatrix[j].length() != cols)
                error_exit("The number of matrix elements is not aligned!\n");
            result[i][j] = plaintextMatrix[j][i];
        }
    }
    return result;
}

/*---------------------- 大矩阵乘法 ------------------------*/
Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::encrypt_largeMat(const Vec< Vec<Integer> > &plaintexts, size_t slot_bits){
    const size_t n_plaintexts = paillierDP.get_packCount(slot_bits);
    size_t flag = plaintexts[0].length() % n_plaintexts; // 剩余向量的数量
    size_t loop = plaintexts[0].length() / n_plaintexts; // 矩阵分块数量

    Vec<Vec<PackedCiphertext>> largeMat;
    largeMat.SetLength(plaintexts.length());

    // 转置矩阵
    // 开始加密并push到密文矩阵
    for (size_t i = 0; i < largeMat.length(); i++)
    {
        if (flag)
            largeMat[i].SetLength(loop+1);
        else
            largeMat[i].SetLength(loop);

        // 分片明文矩阵的cols
        for (size_t j = 0; j < loop; j++)
        {
            const Integer* start = plaintexts[i].begin() + j*n_plaintexts;
            largeMat[i][j] = paillierDP.encrypt_pack(start, start+n_plaintexts, slot_bits);
        }
        if (flag)
        {
            const Integer* start = plaintexts[i].begin() + loop*n_plaintexts;
            largeMat[i][loop] = paillierDP.encrypt_pack(start, start+flag, slot_bits);
        }
    }
    return largeMat;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::encryptTrans_largeMat(const Vec< Vec<Integer> > &plaintexts, size_t slot_bits){
    const size_t n_plaintexts = paillierDP.get_packCount(slot_bits);
    size_t flag = plaintexts.length() % n_plaintexts; // 剩余向量的数量
    size_t loop = plaintexts.length() / n_plaintexts; // 矩阵分块数量

    Vec<Vec<PackedCiphertext>> largeMat;
    largeMat.SetLength(plaintexts[0].length());

    // 转置矩阵
    auto transPlaint = transpose(plaintexts);
    // 开始加密并push到密文矩阵
    for (size_t i = 0; i < largeMat.length(); i++)
    {
        if (flag)
            largeMat[i].SetLength(loop+1);
        else
            largeMat[i].SetLength(loop);

        // 分片明文矩阵的cols
        for (size_t j = 0; j < loop; j++)
        {
            Integer* start = transPlaint[i].begin() + j*n_plaintexts;
            largeMat[i][j] = paillierDP.encrypt_pack(start, start+n_plaintexts, slot_bits);
        }
        if (flag)
        {
            Integer* start = transPlaint[i].begin() + loop*n_plaintexts;
            largeMat[i][loop] = paillierDP.encrypt_pack(start, start+flag, slot_bits);
        }
    }
    return largeMat;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMul_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    for (size_t i = 0; i < result.length(); i++)//结果密文的rows
    {
        for (size_t j = 0; j < result[i].length(); j++)//结果密文的cols = plaintexts的cols
        {   
            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k ==0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else
                    result[i][j].data += ciphertexts[k][i].data * plaintexts[k][j];
            }
        }
    }
    return result;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMulFast_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts, const size_t plaintextBits){
    
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    Vec<Vec<Vec<Ciphertext>>> preMat;
    preMat.SetLength(ciphertexts.length());
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    for (size_t i = 0; i < ciphertexts[0].length(); i++) //结果密文的rows
    {
        for (size_t j = 0; j < plaintexts[0].length(); j++) //结果密文的cols = plaintexts的cols
        {
            if (0 == i && 0 == j)
            {//预计算c,2c,4c,8c,...
                for (size_t t = 0; t < ciphertexts.length(); t++)
                {
                    preMat[t].SetLength(ciphertexts[t].length());
                    for (size_t tt = 0; tt < ciphertexts[t].length(); tt++)
                    {
                        preMat[t][tt].SetLength(plaintextBits);
                        for (size_t ttt = 0; ttt < plaintextBits; ttt++)
                        {
                            preMat[t][tt][ttt] = ciphertexts[t][tt].data * Integer(2).pow(ttt);
                        }
                    }
                }
            }
            
            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k ==0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else{
                    long tmp = plaintexts[k][j].to_long();
                    size_t bits = plaintexts[k][j].size_bits();
                    if (bits > plaintextBits)
                    {
                        error_exit("明文矩阵的值超出范围\n");
                    }
                    if (tmp >= 0)
                    {
                        for (size_t y = 0; y < bits; y++){ 
                            if (tmp & 0x1){
                                result[i][j].data += preMat[k][i][y];
                            }
                            tmp >>= 1;
                        }
                    }else
                    {
                        tmp *= -1;
                        for (size_t y = 0; y < bits; y++){ 
                            if (tmp & 0x1){
                                result[i][j].data -= preMat[k][i][y];
                            }
                            tmp >>= 1;
                        }
                    }
                }                
            }
        }
    }
    return result;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMulFast_signed_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts, const size_t plaintextBits){
    
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    Vec<Vec<Vec<Ciphertext>>> preMat;
    preMat.SetLength(ciphertexts.length());
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    for (size_t i = 0; i < ciphertexts[0].length(); i++) //结果密文的rows
    {
        for (size_t j = 0; j < plaintexts[0].length(); j++) //结果密文的cols = plaintexts的cols
        {
            if (0 == i && 0 == j)
            {//预计算c,2c,4c,8c,...
                for (size_t t = 0; t < ciphertexts.length(); t++)
                {
                    preMat[t].SetLength(ciphertexts[t].length());
                    for (size_t tt = 0; tt < ciphertexts[t].length(); tt++)
                    {
                        preMat[t][tt].SetLength(2*(plaintextBits));
                        for (size_t ttt = 0; ttt < plaintextBits; ttt++)
                        {
                            preMat[t][tt][ttt] = ciphertexts[t][tt].data * Integer(2).pow(ttt);
                            preMat[t][tt][ttt + plaintextBits] = ciphertexts[t][tt].data * Integer(2).pow(ttt * (-1));
                        }
                    }
                }
            }
            
            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k ==0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else{
                    long tmp = plaintexts[k][j].to_long();
                    size_t bits = plaintexts[k][j].size_bits();
                    if (bits > plaintextBits)
                    {
                        error_exit("明文矩阵的值超出范围\n");
                    }
                    if (tmp >= 0)
                    {
                        for (size_t y = 0; y < bits; y++){ 
                            if (tmp & 0x1){
                                result[i][j].data += preMat[k][i][y];
                            }
                            tmp >>= 1;
                        }
                    }else
                    {
                        tmp *= -1;
                        for (size_t y = 0; y < bits; y++){ 
                            if (tmp & 0x1){
                                result[i][j].data -= preMat[k][i][y+plaintextBits-1];
                            }
                            tmp >>= 1;
                        }
                    }
                }                
            }
        }
    }
    return result;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMulBinary_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<bool> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    
    for (size_t i = 0; i < result.length(); i++)//结果密文的rows
    {
        for (size_t j = 0; j < result[i].length(); j++)//结果密文的cols = plaintexts的cols
        {
            
            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k == 0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else if (plaintexts[k][j])
                    result[i][j].data += ciphertexts[k][i].data;
            }
        }
    }
    return result;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMulBinary_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<Integer> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    for (size_t i = 0; i < result.length(); i++)//结果密文的rows
    {
        for (size_t j = 0; j < result[i].length(); j++)//结果密文的cols = plaintexts的cols
        {
            
            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k == 0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else 
                {
                    if (plaintexts[k][j] == 1)
                    {
                        result[i][j].data += ciphertexts[k][i].data;
                    } 
                    else if (plaintexts[k][j] != 0)
                    {
                        error_exit("明文矩阵的值超出范围\n");
                    }
                }
            }
        }
    }
    return result;
}

Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMulAbsOne_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, Vec<Vec<short> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    Vec<Vec<Vec<Ciphertext>>> preMat;
    preMat.SetLength(ciphertexts.length());
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    for (size_t i = 0; i < result.length(); i++)//结果密文的rows
    {
        for (size_t j = 0; j < result[i].length(); j++)//结果密文的cols = plaintexts的cols
        {
            if (0 == i && 0 == j)
            {//预计算c,2c,4c,8c,...
                for (size_t t = 0; t < ciphertexts.length(); t++)
                {
                    preMat[t].SetLength(ciphertexts[t].length());
                    for (size_t tt = 0; tt < ciphertexts[t].length(); tt++)
                    {
                        preMat[t][tt].SetLength(2);
                        preMat[t][tt][0] = ciphertexts[t][tt].data * Integer(1);
                        preMat[t][tt][1] = ciphertexts[t][tt].data * Integer(-1);
                    }
                }
            }

            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k ==0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else{
                    if (plaintexts[k][j] ==1)
                        result[i][j].data += preMat[k][i][0];
                    else if (plaintexts[k][j] == -1)
                        result[i][j].data += preMat[k][i][1];
                    else
                        error_exit("明文矩阵的值超出范围\n");
                }    
            }
        }
    }
    return result;
}
Vec<Vec<PackedCiphertext>> opaillierlib::PaillierMat::matrixMulAbsOne_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts, const Vec<Vec<Integer> > plaintexts){
    if (plaintexts.length() < 1)
        error_exit("plaintextMatrix has no elements!\n");
    if (ciphertexts.length() < 1)
        error_exit("ciphertextMatrix has no elements!\n");
    Vec<Vec<Vec<Ciphertext>>> preMat;
    preMat.SetLength(ciphertexts.length());
    Vec<Vec<PackedCiphertext>> result;
    result.SetLength(ciphertexts[0].length());
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(plaintexts[0].length());
    }
    for (size_t i = 0; i < result.length(); i++)//结果密文的rows
    {
        for (size_t j = 0; j < result[i].length(); j++)//结果密文的cols = plaintexts的cols
        {
            if (0 == i && 0 == j)
            {//预计算c,2c,4c,8c,...
                for (size_t t = 0; t < ciphertexts.length(); t++)
                {
                    preMat[t].SetLength(ciphertexts[t].length());
                    for (size_t tt = 0; tt < ciphertexts[t].length(); tt++)
                    {
                        preMat[t][tt].SetLength(2);
                        preMat[t][tt][0] = ciphertexts[t][tt].data * Integer(1);
                        preMat[t][tt][1] = ciphertexts[t][tt].data * Integer(-1);
                    }
                }
            }

            for (size_t k = 0; k < ciphertexts.length(); k++) 
            {
                if (k ==0)
                    result[i][j] = PackedCiphertext(ciphertexts[k][i].data * plaintexts[k][j],ciphertexts[k][i].n_plaintexts, ciphertexts[k][i].plaintext_bits);
                else{
                    if (plaintexts[k][j] == 1)
                        result[i][j].data += preMat[k][i][0];
                    else if (plaintexts[k][j] == -1)
                        result[i][j].data += preMat[k][i][1];
                    else{
                        std::cout << "plaintexts[k][j]=" << plaintexts[k][j] << std::endl;
                        error_exit("明文矩阵的值超出范围\n");
                    }
                        
                }    
            }
        }
    }
    return result;
}

Vec< Vec<Integer> > opaillierlib::PaillierMat::decrypt_largeMat(const Vec<Vec<PackedCiphertext>> &ciphertexts){
    Vec< Vec<Integer> > results;
    for (size_t i = 0; i < ciphertexts.length(); i++)
    {
        auto tmp = this->decrypt(ciphertexts[i]);
        results.append(tmp);   
    }
    return results;
}