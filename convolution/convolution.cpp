/*
 作者：吴伟彬
 联系方式：wuweibin@nuaa.edu.cn
 */

#include "convolution.h"

using namespace opaillierlib;

PaillierNN::PaillierNN(const size_t key_size_bits):paillierMat(key_size_bits){

}

PaillierMat& PaillierNN::getPaillierMat(){
    return this->paillierMat;
}

PaillierNN::MatrixSize PaillierNN::getOutSize_conv(Image &image, const FilterSize filterSize, const Stride stride, PaddingType padType){
    MatrixSize ms;
    switch (padType)
    {
    case VALID: //"VALID" = without padding
        ms.width = ceil((image.width - filterSize.width + 1)/stride.width);
        ms.height = ceil((image.height - filterSize.height + 1)/stride.height);
        ms.channels = filterSize.filterNum;
        break;
    case SAME: // "SAME" = with zero padding
        error_exit("暂时未实现 `SAME` 模式的Padding!\n");
        break;
    default:
        error_exit("padding模式错误!\n");
        break;
    }
    return ms;
}

Vec<Vec<Integer>> PaillierNN::encode(Image &image, const FilterSize filterSize, const Stride stride, PaddingType padType, const size_t slotBits){
    Vec<Vec<Integer>> result;
    result.SetLength(filterSize.channels * filterSize.height * filterSize.width);
    size_t out_w,out_h;

    // 计算每个通道需要输出的行列大小
    switch (padType)
    {
    case VALID: // "VALID" = without padding
        out_w = ceil((image.height-filterSize.height + 1)/stride.height); //ceil()->向上取整
        out_h = ceil((image.width-filterSize.width + 1)/stride.width); //ceil()->向上取整
        break;
    case SAME: // "SAME" = with zero padding
        out_w = ceil((image.height-filterSize.height)/stride.height); //ceil()->向上取整
        out_h = ceil((image.width-filterSize.width)/stride.width); //ceil()->向上取整
        break;
    default:
        break;
    }
    size_t length = out_h * out_w;
    for (size_t i = 0; i < result.length(); i++)
    {
        result[i].SetLength(length);
    }
    
    for (size_t ch = 0; ch < filterSize.channels; ch++) //channel
    {
        for (size_t h = 0; h < filterSize.height; h++) //height
        {
            for (size_t w = 0; w < filterSize.width; w++) //width
            {
                switch (padType)
                {
                case VALID: //"VALID" = without padding
                    for (size_t i = 0; i < out_h; i++)
                    {
                        for (size_t j = 0; j < out_w; j++)
                        {
                            result[ch*filterSize.height*filterSize.width + h*filterSize.width + w][i*out_w + j] = image.data[ch][i*stride.height+h][j*stride.width+w];
                        }
                    }
                    
                    break;
                case SAME: // "SAME" = with zero padding
                    error_exit("暂时未实现 `SAME` 模式的Padding!\n");
                    break;
                default:
                    error_exit("padding模式错误!\n");
                    break;
                }
                
            }
        }
    }
    
    return result;
}

PaillierNN::IMGPackedCiphertext PaillierNN::encrypt(Vec<Vec<Integer>> image, size_t slotBits){
    return paillierMat.encrypt_largeMat(image, slotBits);
    // return paillierMat.encryptTrans_largeMat(paillierMat.transpose(image), slotBits);
}

Vec<Vec<Vec<Integer>>> PaillierNN::decrypt(PaillierNN::IMGPackedCiphertext ciphertexts, PaillierNN::MatrixSize outMatSize){
    auto res = paillierMat.decrypt_largeMat(ciphertexts);
    Vec<Vec<Vec<Integer>>> result;
    result.SetLength(outMatSize.channels);
    for (size_t ch = 0; ch < result.length(); ch++)
    {
        result[ch].SetLength(outMatSize.height);
        for (size_t h = 0; h < outMatSize.height; h++)
        {
            result[ch][h].SetLength(outMatSize.width);
            for (size_t w = 0; w < outMatSize.width; w++)
            {
                result[ch][h][w] = res[h*outMatSize.width+w][ch];
            }
        }
    }
    return result;
}

PaillierNN::IMGPackedCiphertext PaillierNN::convolution(IMGPackedCiphertext imgCiphertexts, Filters filters){
    // Vec<PaillierNN::IMGPackedCiphertext> result;
    // result.SetLength(filters.length());
    Vec<Vec<Integer> > plaintexts;
    plaintexts.SetLength(filters[0].channels * filters[0].height * filters[0].width);
    for (size_t i = 0; i < plaintexts.length(); i++)
    {
        plaintexts[i].SetLength(filters.length());   
    }
    size_t col =0;
    for (auto &&filter : filters)
    {
        for (size_t ch = 0; ch < filter.channels; ch++)
        {
            for (size_t h = 0; h < filter.height; h++)
            {
                for (size_t w = 0; w < filter.width; w++)
                {
                    plaintexts[ch*filter.height*filter.width + h*filter.width + w][col] = filter.data[ch][h][w];
                }
            }
        }
        col ++;
    }
    
    return paillierMat.matrixMul_largeMat(imgCiphertexts,plaintexts);
}

PaillierNN::IMGPackedCiphertext PaillierNN::convolutionFast(PaillierNN::IMGPackedCiphertext imgCiphertexts, Filters filters, size_t plaintextBits){
    Vec<Vec<Integer> > plaintexts;
    plaintexts.SetLength(filters[0].channels * filters[0].height * filters[0].width);
    for (size_t i = 0; i < plaintexts.length(); i++)
    {
        plaintexts[i].SetLength(filters.length());   
    }
    size_t col = 0;
    for (auto &&filter : filters)
    {
        for (size_t ch = 0; ch < filter.channels; ch++)
        {
            for (size_t h = 0; h < filter.height; h++)
            {
                for (size_t w = 0; w < filter.width; w++)
                {
                    plaintexts[ch*filter.height*filter.width + h*filter.width + w][col] = filter.data[ch][h][w];
                }
            }
        }
        col ++;
    }
    
    return paillierMat.matrixMulFast_largeMat(imgCiphertexts,plaintexts,plaintextBits);
}


PaillierNN::IMGPackedCiphertext PaillierNN::convolution_Binary(PaillierNN::IMGPackedCiphertext imgCiphertexts, PaillierNN::Filters filters){
    // Vec<PaillierNN::IMGPackedCiphertext> result;
    // result.SetLength(filters.length());
    Vec<Vec<bool> > plaintexts;
    plaintexts.SetLength(filters[0].channels * filters[0].height * filters[0].width);
    for (size_t i = 0; i < plaintexts.length(); i++)
    {
        plaintexts[i].SetLength(filters.length());   
    }
    size_t col =0;
    for (auto &&filter : filters)
    {
        for (size_t ch = 0; ch < filter.channels; ch++)
        {
            for (size_t h = 0; h < filter.height; h++)
            {
                for (size_t w = 0; w < filter.width; w++)
                {
                    if (filter.data[ch][h][w] == 0)
                        plaintexts[ch*filter.height*filter.width + h*filter.width + w][col] = false;
                    else if (filter.data[ch][h][w] == 1)
                        plaintexts[ch*filter.height*filter.width + h*filter.width + w][col] = true;
                    else
                        error_exit("plaintext value is not 1 or 0!\n");
                    // plaintexts[ch*filter.height*filter.width + h*filter.width + w][col] = filter.data[ch][h][w];
                }
            }
        }
        col ++;
    }
    
    return paillierMat.matrixMulBinary_largeMat(imgCiphertexts,plaintexts);
}
PaillierNN::IMGPackedCiphertext PaillierNN::convolution_AbsOne(PaillierNN::IMGPackedCiphertext imgCiphertexts, PaillierNN::Filters filters){
    // Vec<PaillierNN::IMGPackedCiphertext> result;
    // result.SetLength(filters.length());
    Vec<Vec<short> > plaintexts;
    plaintexts.SetLength(filters[0].channels * filters[0].height * filters[0].width);
    for (size_t i = 0; i < plaintexts.length(); i++)
    {
        plaintexts[i].SetLength(filters.length());   
    }
    size_t col =0;
    for (auto &&filter : filters)
    {
        for (size_t ch = 0; ch < filter.channels; ch++)
        {
            for (size_t h = 0; h < filter.height; h++)
            {
                for (size_t w = 0; w < filter.width; w++)
                {
                    plaintexts[ch*filter.height*filter.width + h*filter.width + w][col] = filter.data[ch][h][w].to_int();
                }
            }
        }
        col ++;
    }
    
    return paillierMat.matrixMulAbsOne_largeMat(imgCiphertexts,plaintexts);
}
