#include "benchmark.h"
#include "../matrixMul/matMul.h"
#include "../utils/myRandom.h"



using namespace opaillierlib;
using namespace std;

template <typename T>
void printVector(const Vec<T>& values){
    for (size_t i = 0; i < values.length(); i++)
    {
        std::cout << values[i] << " ";
    }
    std::cout << endl;
}

void benchmark_PaillierVM(size_t times, size_t key_size, size_t vecSize, size_t cipherBits, size_t plaintBits)
{
    NewPaillier pai(key_size);
    pai.generate_keys();
    Vec<Integer> A, B, R;
    Vec<Ciphertext> cts;
    cts.SetLength(vecSize);
    Ciphertext ct_result;

    std::cout << "keySize=" << key_size << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << endl;
    std::cout << "vector size:" << vecSize << std::endl;
    
    //init TIMER
    TIMER_INIT(encrypt,times)
    TIMER_INIT(dotProduct,times)
    TIMER_INIT(decrypt,times)

    for (size_t t = 0; t < times; t++)
    {
        genVector_signed(A,vecSize,cipherBits);
        genVector_signed(B,vecSize,cipherBits);
        TIMER_CONTINUE(encrypt)
        for (size_t i = 0; i < A.length(); i++)
        {
            cts[i] = pai.encrypt(A[i]);
        }
        TIMER_PAUSE(encrypt)

        TIMER_CONTINUE(dotProduct)
        ct_result = cts[0] * B[0];
        for (size_t i = 1; i < B.length(); i++)
        {
            ct_result += cts[i] * B[i];
        }
        TIMER_PAUSE(dotProduct)

        TIMER_CONTINUE(decrypt)
        auto res = pai.decrypt(ct_result);
        TIMER_PAUSE(decrypt)
    }
    TIMER_PRINT(encrypt)
    TIMER_PRINT(dotProduct)
    TIMER_PRINT(decrypt)
    TIMER_TOTAL(encrypt,dotProduct,decrypt)    
}

void benchamrk_VM(size_t times, size_t key_size, size_t vecSize, size_t cipherBits, size_t plaintBits){
    std::cout << "start benchmark_Vector Multiplication:" << std::endl;
    Vec<Integer> A, B, R;
    Vec<bool> B_bool;
    Vec<short> B_short;
    opaillierlib::PaillierDP paiDP(key_size);
    int slotBits = (plaintBits+1)+(cipherBits+1)+log2(vecSize)+1;
    // int slotBits = 31;
    // if (plaintBits == 0)
    // {
    //     slotBits += 1;
    // }
    // if (cipherBits == 0)
    // {
    //     slotBits += 1;
    // }

    // auto slotNum = paiDP.get_packCount(slotBits);

    // genVector_signed(R,slotNum,cipherBits);
    // auto R_cts = paiDP.getPaillier().encrypt_pack(R,slotBits);

    // cout << "encode:" << endl;
    std::cout << "keySize=" << key_size <<"; slotBits=" << slotBits << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << endl;

    std::cout << "vector size:" << vecSize << std::endl;
    
    std::cout << "check Slots: " << paiDP.checkSlotSize(slotBits,cipherBits,plaintBits,vecSize) << endl;
    
    //init TIMER
    TIMER_INIT(encrypt,times)
    TIMER_INIT(dotProduct,times)
    TIMER_INIT(decrypt,times)
    for (size_t i = 0; i < times; i++)
    {
        genVector_signed(A,vecSize,cipherBits);
        Integer sum = 0;
        switch (plaintBits)
        {
        case 0: //binary
            genVector_binary(B_bool,vecSize);
            for (size_t j = 0; j < vecSize; j++)
            {
                if (B_bool[j])
                {
                    sum += A[j];
                }
            }
            break;
        case 1://absOne
            genVector_absOne(B_short,vecSize);
            for (size_t j = 0; j < vecSize; j++)
            {
                if (B_short[j] == 1)
                {
                    sum += A[j];
                }
                else{
                    sum -= A[j];
                }
            }
            break;
        default:
            genVector_signed(B,vecSize,cipherBits);
            for (size_t j = 0; j < vecSize; j++)
            {
                sum += A[j]*B[j];
            }
            break;
        }

        // cout << "encrypt" << endl;
        TIMER_CONTINUE(encrypt)
        auto cts = paiDP.encrypt(A,slotBits);
        TIMER_PAUSE(encrypt)
        
        PackedCiphertext resultCt;
        switch (plaintBits)
        {
        case 0: //binary
            TIMER_CONTINUE(dotProduct)
            resultCt = paiDP.dotProductBinary(cts,B_bool);
            // resultCt.data += R_cts.data;
            TIMER_PAUSE(dotProduct)
            break;
        case 1://absOne
            TIMER_CONTINUE(dotProduct)
            resultCt = paiDP.dotProductAbsOne(cts,B_short);
            // resultCt.data += R_cts.data;
            TIMER_PAUSE(dotProduct)
            break;
        default:
            TIMER_CONTINUE(dotProduct)
            resultCt = paiDP.dotProduct(cts,B);
            // resultCt.data += R_cts.data;
            TIMER_PAUSE(dotProduct)
            break;
        }

        // cout<<"dec"<< endl;
        
        TIMER_CONTINUE(decrypt)
        auto results = paiDP.decrypt(resultCt,vecSize);
        TIMER_PAUSE(decrypt)

        if (results != sum)//R[slotNum/2]
        {
            std::cout << "failed"  << "; sum=" << sum << "; results=" << results << endl;
        }
        
    }
    TIMER_PRINT(encrypt)
    TIMER_PRINT(dotProduct)
    TIMER_PRINT(decrypt)
    TIMER_TOTAL(encrypt,dotProduct,decrypt)

}

//plaintBits: 0=binary 1=absOne x=x-bits
void benchmark_SMP(size_t times, size_t key_size, size_t mat_n1, size_t mat_n2, size_t mat_n3, size_t cipherBits, size_t plaintBits){
    std::cout << "start benchmark_SMP:" << std::endl;
    Vec<Vec<Integer>> A,B,R;
    opaillierlib::PaillierMat paillierMat(key_size);
    int slotBits = plaintBits+cipherBits+log2(mat_n2)+1;
    if (plaintBits == 0)
    {
        slotBits += 1;
    }
    if (cipherBits == 0)
    {
        slotBits += 1;
    }
    
    switch (cipherBits)
    {
    case 0: //binary
        genMatrix_binary(R,mat_n1,mat_n3);
        break;
    case 1://absOne
        genMatrix_absOne(R,mat_n1,mat_n3);
        break;
    default:
        genMatrix_signed(R,mat_n1,mat_n3,cipherBits);
        break;
    }

    auto R_cts = paillierMat.encryptTrans_largeMat(R,slotBits);

    // cout << "encode:" << endl;
    std::cout << "keySize=" << key_size <<"; slotBits=" << slotBits << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << endl;

    std::cout << "Matrix A:(" << mat_n1 << "," << mat_n2 \
            << "); Matrix B:(" << mat_n2 << "," << mat_n3 << ")" << endl;
    
    std::cout << "check Slots: " << paillierMat.checkSlotSize(slotBits,cipherBits,plaintBits,mat_n2) << endl;
    
    //init TIMER
    TIMER_INIT(encrypt,times)
    TIMER_INIT(matrixMulFast_largeMat,times)
    TIMER_INIT(decrypt,times)

    for (size_t i = 0; i < times; i++)
    {
        // gen random data
        switch (cipherBits)
        {
        case 0: //binary
            genMatrix_binary(A,mat_n1,mat_n2);
            break;
        case 1://absOne
            genMatrix_absOne(A,mat_n1,mat_n2);
            break;
        default:
            genMatrix_signed(A,mat_n1,mat_n2,cipherBits);
            break;
        }

        switch (plaintBits)
        {
        case 0: //binary
            genMatrix_binary(B,mat_n2,mat_n3);
            break;
        case 1://absOne
            genMatrix_absOne(B,mat_n2,mat_n3);
            break;
        default:
            genMatrix_signed(B,mat_n2,mat_n3,plaintBits);
            break;
        }

        // cout << "encrypt" << endl;
        TIMER_CONTINUE(encrypt)
        auto cts = paillierMat.encryptTrans_largeMat(A,slotBits);
        TIMER_PAUSE(encrypt)
        if (i == 0)
        {
            std::cout << "client sent " << cts.length() * cts[0].length() << " ciphertexts" << std::endl;
        }
        
        Vec<Vec<PackedCiphertext>> resultCt;
        switch (plaintBits)
        {
        case 0: //binary
            TIMER_CONTINUE(matrixMulFast_largeMat)
            resultCt = paillierMat.matrixMulBinary_largeMat(cts,B);
            for (size_t h = 0; h < resultCt.length(); h++)
            {
                for (size_t w = 0; w < resultCt[h].length(); w++)
                {
                    resultCt[h][w].data += R_cts[w][h].data;
                }
            }
            TIMER_PAUSE(matrixMulFast_largeMat)
            break;
        case 1://absOne
            TIMER_CONTINUE(matrixMulFast_largeMat)
            resultCt = paillierMat.matrixMulAbsOne_largeMat(cts,B);
            for (size_t h = 0; h < resultCt.length(); h++)
            {
                for (size_t w = 0; w < resultCt[h].length(); w++)
                {
                    resultCt[h][w].data += R_cts[w][h].data;
                }
            }
            TIMER_PAUSE(matrixMulFast_largeMat)
            break;
        default:
            TIMER_CONTINUE(matrixMulFast_largeMat)
            resultCt = paillierMat.matrixMulFast_signed_largeMat(cts,B,plaintBits);
            for (size_t h = 0; h < resultCt.length(); h++)
            {
                for (size_t w = 0; w < resultCt[h].length(); w++)
                {
                    resultCt[h][w].data += R_cts[w][h].data;
                }
            }
            TIMER_PAUSE(matrixMulFast_largeMat)
            break;
        }
        if (0 == i)
        {
            std::cout << "server sent " << resultCt.length() * resultCt[0].length() << " ciphertexts" << std::endl;
        }
        
        
        TIMER_CONTINUE(decrypt)
        auto results = paillierMat.decrypt_largeMat(resultCt);
        TIMER_PAUSE(decrypt)
    }
    TIMER_PRINT(encrypt)
    TIMER_PRINT(matrixMulFast_largeMat)
    TIMER_PRINT(decrypt)
    TIMER_TOTAL(encrypt,matrixMulFast_largeMat,decrypt)

}

void benchmark_conv(size_t times, size_t key_size,
                     size_t image_channels, size_t image_height, size_t image_width,
                     size_t filter_num, size_t filter_channels, size_t filter_height, size_t filter_width,
                     size_t stride_height, size_t stride_width, 
                     size_t cipherBits, size_t plaintBits){

    PaillierNN::Image image;
    PaillierNN::Filters filters;

    PaillierNN::Stride stride;
    stride.height = stride_height;
    stride.width = stride_width;

    PaillierNN::FilterSize filterSize;
    filterSize.channels = filter_channels;
    filterSize.height = filter_height;
    filterSize.width = filter_width;
    filterSize.filterNum = filter_num;

    PaillierNN::PaddingType padType = PaillierNN::PaddingType::VALID;

    PaillierNN paillierNN(key_size);

    size_t slotBits = plaintBits+cipherBits+log(filter_channels * filter_height * filter_width)+1;
    if (plaintBits == 0)
    {
        slotBits += 1;
    }
    if (cipherBits == 0)
    {
        slotBits += 1;
    }

    auto w = ceil((image_width - filterSize.width + 1)/stride.width);
    auto h = ceil((image_height - filterSize.height + 1)/stride.height);
    // auto co = filter_num;
    // PaillierNN::Image  R;
    // genImg_signed(R, co, h, w, cipherBits);

    // auto imgCt = paillierNN.encrypt(R.data,slotBits);

    // cout<<"image.width" << image.width << "; filterSize.width: "<< filterSize.width << "; stride.width:" << stride.width << endl;
    Vec<Vec<Integer>> R;
    switch (cipherBits)
    {
    case 0: //binary
        genMatrix_binary(R,w*h,filter_num);
        break;
    case 1://absOne
        genMatrix_absOne(R,w*h,filter_num);
        break;
    default:
        genMatrix_signed(R,w*h,filter_num,cipherBits);
        break;
    }

    auto R_cts = paillierNN.getPaillierMat().encryptTrans_largeMat(R,slotBits);
    // cout<<"R_cts.nums: " << R_cts.length() << "; R_cts.nums[0]: " << R_cts[0].length() << "; filter_num:" <<filter_num << endl;

    // cout << "encode:" << endl;
    std::cout << "keySize=" << key_size <<"; slotBits=" << slotBits << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << std::endl;
    cout << "Image("<< image_height << "*" << image_width <<"," << image_channels <<"); Filters("<< filter_height << "*"<< filter_width << ","  << filter_num <<"); Stride(" << stride_height << "," << stride_width << ")"<< endl;

    //init TIMER
    TIMER_INIT(encode,times)
    TIMER_INIT(encrypt,times)
    TIMER_INIT(convolution,times)
    TIMER_INIT(decrypt,times)

    for (size_t i = 0; i < times; i++)
    {
        // gen random data
        switch (cipherBits)
        {
        case 0: //binary
            genImg_binary(image, image_channels, image_height, image_width);
            break;
        case 1://absOne
            genImg_absOne(image, image_channels, image_height, image_width);
            break;
        default:
            genImg_signed(image, image_channels, image_height, image_width, cipherBits);
            break;
        }

        switch (plaintBits)
        {
        case 0: //binary
            genFilters_binary(filters, filter_num, filter_channels, filter_height, filter_width);
            break;
        case 1://absOne
            genFilters_absOne(filters, filter_num, filter_channels, filter_height, filter_width);
            break;
        default:
            genFilters_unsigned(filters, filter_num, filter_channels, filter_height, filter_width, plaintBits);
            break;
        }

        // cout << "encrypt" << endl;
        TIMER_CONTINUE(encode)
        auto encodeImg = paillierNN.encode(image,filterSize,stride,padType,slotBits);
        TIMER_PAUSE(encode)

        TIMER_CONTINUE(encrypt)
        auto imgCt = paillierNN.encrypt(encodeImg,slotBits);
        TIMER_PAUSE(encrypt)
        if (i == 0)
        {
            std::cout << "client sent " << imgCt.length() * imgCt[0].length() << " ciphertexts->size="<< (1.0 * imgCt.length() * imgCt[0].length() *key_size * 2)/ (1024*8) << "KB"<< std::endl;
        }
        
        PaillierNN::IMGPackedCiphertext  resultCt;
        switch (plaintBits)
        {
        case 0: //binary
            TIMER_CONTINUE(convolution)
            resultCt = paillierNN.convolution_Binary(imgCt,filters);
            for (size_t h = 0; h < resultCt.length(); h++)
            {
                for (size_t w = 0; w < resultCt[h].length(); w++)
                {
                    resultCt[h][w].data += R_cts[w][h].data;
                }
            }
            TIMER_PAUSE(convolution)
            break;
        case 1://absOne
            TIMER_CONTINUE(convolution)
            resultCt = paillierNN.convolution_AbsOne(imgCt,filters);
            for (size_t h = 0; h < resultCt[0].length(); h++)
            {
                for (size_t w = 0; w < resultCt.length(); w++)
                {
                    resultCt[h][w].data += R_cts[w][h].data;
                }
            }
            TIMER_PAUSE(convolution)
            break;
        default:
            TIMER_CONTINUE(convolution)
            resultCt = paillierNN.convolutionFast(imgCt,filters,plaintBits);
            // for (size_t h = 0; h < resultCt.length(); h++)
            // {
            //     for (size_t w = 0; w < resultCt[h].length(); w++)
            //     {
            //         resultCt[h][w].data += R_cts[w][h].data;
            //     }
            // }
            TIMER_PAUSE(convolution)
            break;
        }
        if (0 == i)
        {
            std::cout << "server sent " << resultCt.length() * resultCt[0].length() << " ciphertexts->size="<< (1.0 * resultCt.length() * resultCt[0].length() * key_size * 2) / (1024*8) << "KB"<< std::endl;
            // std::cout << "resultCt.length(): " << resultCt.length() << "; resultCt[0].length(): " << resultCt[0].length() << std::endl;
        }
        
        auto matSize = paillierNN.getOutSize_conv(image,filterSize,stride,padType);

        TIMER_CONTINUE(decrypt)
        auto result = paillierNN.decrypt(resultCt,matSize);
        TIMER_PAUSE(decrypt)
    }
    TIMER_PRINT(encode)
    TIMER_PRINT(encrypt)
    TIMER_PRINT(convolution)
    TIMER_PRINT(decrypt)
    TIMER_TOTAL_CC(encode,encrypt,convolution,decrypt)
}