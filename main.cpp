#include "ophelib/util.h"
#include <stdio.h>
#include <time.h>
#include <random>

#include "dotProduct/dotProduct.h"
#include "paillier/newPaillier.h"
#include "matrixMul/matMul.h"
#include "convolution/convolution.h"
#include "ophelib/paillier_fast.h"

#include "utils/myRandom.h"
#include "benchmarks/benchmark.h"

#include "ophelib/random.h"

using namespace std;
using namespace ophelib;
using namespace opaillierlib;

#define  TIMES 10
size_t imageChannels = 3;
size_t imageHeight = 224;
size_t imageWidth = 224;

size_t filterNum = 64;
size_t filterHeight = 7;
size_t filterWidth = 7;

size_t strideHeight = 2;
size_t strideWidth = 2;

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


template <typename T>
void printVector(const Vec<T>& values);
void printBool(const Vec<bool>& values);
void printShort(const Vec<short>& values);

void testNewPaillier();

void testDP_vec_vec();
void testDP_vec_vec_Binary();
void testDP_vec_vec_AbsOne();

void testDP_mat_vec();
void testDP_mat_vec_Binary();
void testDP_mat_vec_AbsOne();


void testMatMul();
void testMatMul_Binary();
void testMatMul_AbsOne();

void testLargeMatMul();
void testLargeMatMul_Binary();
void testLargeMatMul_AbsOne();
void matrixMulFast_largeMat();

void testConvolution();
void testConvolution_Binary();
void testConvolution_AbsOne();
void testConvolutionFast();

struct ConvParams
{
    // image
    size_t cipherBits = 0;
    size_t image_channels;
    size_t image_height;
    size_t image_width;

    // filters
    size_t plaintBits = 0;
    size_t filter_num;
    size_t filter_channels;
    size_t filter_height;
    size_t filter_width;

    // stride
    size_t stride_height;
    size_t stride_width;
};


int main(int argc, char const *argv[])
{
    initTime();
    /*-----------测试原Paillier性能---------*/
    // PaillierFast pail(2048);
    // NewPaillier pai(2048);
    // pai.generate_keys();
    // pail.generate_keys();
    // size_t n = 100;
    // TIMER_INIT(encrypt_Paillier,n)
    // TIMER_INIT(decrypt_Paillier,n)       
    // TIMER_INIT(encrypt_NewPaillier,n)
    // TIMER_INIT(decrypt_NewPaillier,n)       
    // Random& rand = Random::instance();
    // size_t bits[] = {16,32,64,500,1000,1500,2048} ;
    // for (size_t i = 0; i < n; i++)
    // {
    //     Integer a = rand.rand_int_bits(16);
    //     // cout << a << endl;
    //     TIMER_CONTINUE(encrypt_Paillier)
    //     auto ct = pail.encrypt(a);
    //     TIMER_PAUSE(encrypt_Paillier)

    //     TIMER_CONTINUE(encrypt_NewPaillier)
    //     auto ct1 = pai.encrypt(a);
    //     TIMER_PAUSE(encrypt_NewPaillier)

    //     TIMER_CONTINUE(decrypt_Paillier)
    //     auto res = pail.decrypt(ct);
    //     TIMER_PAUSE(decrypt_Paillier)

    //     TIMER_CONTINUE(decrypt_NewPaillier)
    //     auto res1 = pai.decrypt(ct1);
    //     TIMER_PAUSE(decrypt_NewPaillier)
    // }
    // TIMER_PRINT(encrypt_Paillier)
    // TIMER_PRINT(decrypt_Paillier)
    // TIMER_PRINT(encrypt_NewPaillier)
    // TIMER_PRINT(decrypt_NewPaillier)


    /*----------测试Paillier派生类------------*/
    // testNewPaillier();
    // ophelib::PaillierFast pai(1024);
    // pai.generate_keys();
    // Integer a = 50;
    // auto ct = pai.encrypt(a);
    

    // // // pai.get_pub().g.pow_mod_n(pai.get_pub().n, *pai.get_n2().get());
    // // Integer ret = ((Integer::L(pai.get_fast_mod().get()->pow_mod_n2(ct.data, pai.get_pub().n-1),pai.get_pub().n)) 
    // //             / Integer::L(pai.get_pub().g.pow_mod_n(pai.get_pub().n-1, *pai.get_n2().get()),pai.get_pub().n))% pai.get_pub().n;
    // // cout << "result = " << ret << endl;
    
    // Integer ret = ((Integer::L(pai.get_fast_mod().get()->pow_mod_n2(ct.data, pai.get_pub().n-1),pai.get_pub().n)) 
    //             / (Integer::L(pai.get_fast_mod().get()->pow_mod_n2(pai.get_pub().g, pai.get_pub().n-1),pai.get_pub().n))) 
    //             % pai.get_pub().n;
    // cout << "result = " << ret << endl;

    // opaillierlib::NewPaillier paillier(2048);
    // paillier.generate_keys();
    // Integer a = 50;
    // auto ct = paillier.encrypt(a);
    // // ct *= 5; 
    // // ct += ct;
    // auto result = paillier.decrypt_test(ct);
    // cout << "result = " << result << endl;

    

    // opaillierlib::NewPaillier paillier(2048);
    // // ophelib::PaillierFast pai(2048);
    // // pai.generate_keys();
    // paillier.generate_keys();
    // int itor = 100;
    // TIMER_INIT(encrypt_pack,itor)
    // // TIMER_INIT(encrypt,itor)
    // TIMER_INIT(mul,itor)
    // TIMER_INIT(add_ctct,itor)
    // TIMER_INIT(add_ctpt,itor)
    // TIMER_INIT(decrypt_pack,itor)
    // Vec<Integer> r;
    // Vec<Integer> a,b;
    // genVector_one(b,93);
    
    // for (size_t i = 0; i < itor; i++)
    // {
    //     auto aCt = paillier.encrypt_pack(b.begin(),b.end(),22);
    //     genVector_unsigned(r,2,2000);
    //     genVector_unsigned(a,93,8);

    //     TIMER_CONTINUE(encrypt_pack)
    //     // auto aCt = paillier.encrypt_pack(a.begin(),a.end(),22);
    //     //编码
    //     Integer sum = a[0];
    //     for (auto iter = a.begin()+1; iter < a.end(); iter++) {
    //         if(iter->size_bits() > 22)
    //             error_exit("plaintext size too large!\n");
    //         sum <<= 22;
    //         sum += *iter;
    //     }
    //     // Integer c = paillier.get_pub().g.pow_mod_n(sum,*paillier.get_n2());
    //     // auto aCt = PackedCiphertext(ct0.data *sum,ct0.n_plaintexts,22);
    //     aCt.data *= sum;
    //     // const PackedCiphertext packed = Vector::encrypt_pack(a, 21, pai);
    //     TIMER_PAUSE(encrypt_pack)
    //     // cout << "a:" << endl;
    //     // printVector(a);
    //     // auto res_tmp = paillier.decrypt_pack(aCt);
    //     // cout << "dec:" << endl;
    //     // printVector(res_tmp);

    //     // // auto ct = pai.encrypt(r[0]);
    //     // auto ct = paillier.encrypt(r[0]);
    //     // cout <<"r[0] = " << r[0] << "; r[1] = " << r[1] << endl;
    //     // auto tmp = ct + paillier.get_pub().g.pow_mod_n(r[1],*paillier.get_n2());;
    //     // cout << paillier.decrypt(tmp) << endl;


    //     TIMER_CONTINUE(mul)
    //     auto result = PackedCiphertext(aCt.data * r[1],aCt.n_plaintexts, aCt.plaintext_bits);
    //     // auto result = PackedCiphertext(packed.data * r[1],packed.n_plaintexts, packed.plaintext_bits);
    //     TIMER_PAUSE(mul)

    //     TIMER_CONTINUE(add_ctct)
    //     // result.data = packed.data + result.data;
    //     result.data = aCt.data + result.data;
    //     TIMER_PAUSE(add_ctct)

    //     TIMER_CONTINUE(add_ctpt)
    //     // result.data = packed.data + result.data;
    //     result.data = aCt.data + r[1];
    //     TIMER_PAUSE(add_ctpt)

    //     TIMER_CONTINUE(decrypt_pack)
    //     auto res = paillier.decrypt_pack(result);
    //     // auto res = Vector::decrypt_pack(result, pai);
    //     TIMER_PAUSE(decrypt_pack)
    // }
    // TIMER_PRINT(encrypt_pack)
    // TIMER_PRINT(mul)
    // TIMER_PRINT(add_ctct)
    // TIMER_PRINT(add_ctpt)
    // TIMER_PRINT(decrypt_pack)
    
    /*----------测试`向量*向量`点积------------*/
    // testDP_vec_vec();
    // testDP_vec_vec_Binary();
    // testDP_vec_vec_AbsOne();

    /*----------测试`矩阵*向量`点积------------*/
    // testDP_mat_vec();
    // testDP_mat_vec_Binary();
    // testDP_mat_vec_AbsOne();

    /*----------测试矩阵乘法------------*/
    // testMatMul();
    // testMatMul_Binary();
    // testMatMul_AbsOne();

    /*----------测试大矩阵乘法------------*/
    // testLargeMatMul();
    // testLargeMatMul_Binary();
    // testLargeMatMul_AbsOne();
    // matrixMulFast_largeMat();

    /*----------测试卷积------------*/
    // testConvolution();
    // testConvolutionFast();
    // std::cout << std::endl;
    // testConvolution_Binary();
    // std::cout << std::endl;
    // testConvolution_AbsOne();
    // std::cout << std::endl;

    /*---------- benchmark VM ------------*/
    // size_t times = 20;
    // size_t key_size = 2048;
    // size_t vecSize[] = {8,16,32};// :max:44
    // size_t cipherBits = 8;
    // size_t plaintBits = 8;
    // for (size_t i = 0; i < 3; i++)
    // {
    //     benchamrk_VM(times,key_size,vecSize[i],cipherBits,plaintBits);
    //     // benchmark_PaillierVM(times,key_size,vecSize[i],cipherBits,plaintBits);
    //     cout <<  endl;
    // }
    
    

    /*---------- benchmark SMP ------------*/
    size_t times = 20;
    size_t key_size = 1024;
    size_t mat_n1[] = {128,256,512};
    size_t mat_n2 = 128;
    size_t mat_n3[] = {128,256,512};
    size_t cipherBits[6] = {0,1,2,4,8,16};
    size_t plaintBits[5] = {0,1,4,8,16};
    benchmark_SMP(times,key_size,mat_n1[0],mat_n2,mat_n3[0],cipherBits[4],0);std::cout << endl;
    benchmark_SMP(times,key_size,mat_n1[1],mat_n2,mat_n3[1],cipherBits[4],0);std::cout << endl;
    benchmark_SMP(times,key_size,mat_n1[2],mat_n2,mat_n3[2],cipherBits[4],0);std::cout << endl;
    // for (size_t n = 0; n < 3; n++)
    // {
    //     for (size_t pbit = 0; pbit < 1; pbit++)
    //     {
    //         for (size_t cbit = 4; cbit < 5; cbit++)
    //         {
    //             benchmark_SMP(times,key_size,mat_n1[n],mat_n2,mat_n3[n],cipherBits[cbit],plaintBits[pbit]);
    //             std::cout << endl;
    //         }
    //     }
    // }

    /*------------------ benchmark Convolution ------------------------*/
    // size_t times = 5;
    // size_t key_size = 2048;

    // size_t cipherBits[6] = {0,1,2,4,8,16};
    // size_t plaintBits[5] = {0,1,4,8,16};
    // Vec<ConvParams> convParams;
    // convParams.SetLength(8);

    // convParams[0].image_channels = 1;
    // convParams[0].image_height = 28;
    // convParams[0].image_width = convParams[0].image_height;
    // convParams[0].filter_num = 5;
    // convParams[0].filter_channels = convParams[0].image_channels;
    // convParams[0].filter_height = 5;
    // convParams[0].filter_width = convParams[0].filter_height;
    // convParams[0].stride_height = 2;
    // convParams[0].stride_width = convParams[0].stride_height;

    // convParams[1].image_channels = 32;
    // convParams[1].image_height = 32;
    // convParams[1].image_width = convParams[1].image_height;
    // convParams[1].filter_num = 32;
    // convParams[1].filter_channels = convParams[1].image_channels;
    // convParams[1].filter_height = 3;
    // convParams[1].filter_width = convParams[1].filter_height;
    // convParams[1].stride_height = 2;
    // convParams[1].stride_width = convParams[1].stride_height;

    // convParams[2].image_channels = 128;
    // convParams[2].image_height = 16;
    // convParams[2].image_width = convParams[2].image_height;
    // convParams[2].filter_num = 128;
    // convParams[2].filter_channels = convParams[2].image_channels;
    // convParams[2].filter_height = 3;
    // convParams[2].filter_width = convParams[2].filter_height;
    // convParams[2].stride_height = 2;
    // convParams[2].stride_width = convParams[2].stride_height;

    // convParams[3].image_channels = 3;
    // convParams[3].image_height = 224;
    // convParams[3].image_width = convParams[3].image_height;
    // convParams[3].filter_num = 64;
    // convParams[3].filter_channels = convParams[3].image_channels;
    // convParams[3].filter_height = 7;
    // convParams[3].filter_width = convParams[3].filter_height;
    // convParams[3].stride_height = 2;
    // convParams[3].stride_width = convParams[3].stride_height;

    // convParams[4].image_channels = 3;
    // convParams[4].image_height = 224;
    // convParams[4].image_width = convParams[4].image_height;
    // convParams[4].filter_num = 64;
    // convParams[4].filter_channels = convParams[4].image_channels;
    // convParams[4].filter_height = 3;
    // convParams[4].filter_width = convParams[4].filter_height;
    // convParams[4].stride_height = 2;
    // convParams[4].stride_width = convParams[4].stride_height;

    // convParams[5].image_channels = 256;
    // convParams[5].image_height = 56;
    // convParams[5].image_width = convParams[5].image_height;
    // convParams[5].filter_num = 64;
    // convParams[5].filter_channels = convParams[5].image_channels;
    // convParams[5].filter_height = 1;
    // convParams[5].filter_width = convParams[5].filter_height;
    // convParams[5].stride_height = 1;
    // convParams[5].stride_width = convParams[5].stride_height;

    // convParams[6].image_channels = 64;
    // convParams[6].image_height = 56;
    // convParams[6].image_width = convParams[6].image_height;
    // convParams[6].filter_num = 256;
    // convParams[6].filter_channels = convParams[6].image_channels;
    // convParams[6].filter_height = 1;
    // convParams[6].filter_width = convParams[6].filter_height;
    // convParams[6].stride_height = 1;
    // convParams[6].stride_width = convParams[6].stride_height;

    // convParams[7].image_channels = 3;
    // convParams[7].image_height = 480;
    // convParams[7].image_width = convParams[7].image_height;
    // convParams[7].filter_num = 64;
    // convParams[7].filter_channels = convParams[7].image_channels;
    // convParams[7].filter_height = 3;
    // convParams[7].filter_width = convParams[7].filter_height;
    // convParams[7].stride_height = 2;
    // convParams[7].stride_width = convParams[7].stride_height;

    
    // for (size_t cbit = 4; cbit < 5; cbit++)
    // {
    //     for (size_t pbit = 0; pbit < 1; pbit++)
    //     {
    //         for (size_t i = 0; i < convParams.length(); i++)
    //         {
    //             benchmark_conv(times,key_size,
    //                             convParams[i].image_channels, convParams[i].image_height, convParams[i].image_width,
    //                             convParams[i].filter_num, convParams[i].filter_channels, convParams[i].filter_height, convParams[i].filter_width,
    //                             convParams[i].stride_height, convParams[i].stride_width,
    //                             cipherBits[cbit], 0);//plaintBits[pbit]
    //             std:;cout << endl;
    //         }
    //     }
    // }



    // benchmark_conv(times,key_size,
    //                             convParams[5].image_channels, convParams[5].image_height, convParams[5].image_width,
    //                             convParams[5].filter_num, convParams[5].filter_channels, convParams[5].filter_height, convParams[5].filter_width,
    //                             convParams[5].stride_height, convParams[5].stride_width,
    //                             cipherBits[4], plaintBits[3]);

    cout << "Running on ophelib version " << ophelib::ophelib_version() << ", git ref " << ophelib::ophelib_ref() << endl;
    return 0;
}


/*------------- test Convolution ------------------*/
void testConvolutionFast(){
    std::cout << "start verify testConvolutionFast:" << endl;
    const size_t key_size = 2048;

    PaillierNN::Image image;
    size_t cipherBits = 8;
    size_t image_channels = imageChannels;
    size_t image_height = imageHeight;
    size_t image_width = imageWidth;

    size_t plaintBits = 8;
    PaillierNN::Filters filters;
    size_t filter_num = filterNum;
    size_t filter_channels = image_channels;
    size_t filter_height = filterHeight;
    size_t filter_width = filterWidth;
    
    PaillierNN::Stride stride;
    stride.height = strideHeight;
    stride.width = strideWidth;

    PaillierNN::FilterSize filterSize;
    filterSize.channels = filter_channels;
    filterSize.height = filter_height;
    filterSize.width = filter_width;
    filterSize.filterNum = filter_num;
    
    PaillierNN::PaddingType padType = PaillierNN::PaddingType::VALID;
    
    size_t slotBits = plaintBits+cipherBits+log(filter_channels * filter_height * filter_width)+1;
    
    PaillierNN paillierNN(key_size);

    // cout << "encode:" << endl;
    std::cout << "keySize=" << key_size <<"; slotBits=" << slotBits << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << endl;
    cout << "Image("<< image_height << "*" << image_width <<"," << image_channels <<"); Filters("<< filter_height << "*"<< filter_width << ","  << filter_num <<")" << endl;
    //init TIMER
    TIMER_INIT(encode,TIMES)
    TIMER_INIT(encrypt,TIMES)
    TIMER_INIT(convolutionFast,TIMES)
    TIMER_INIT(decrypt,TIMES)

    for (size_t i = 0; i < TIMES; i++)
    {
        // gen random data
        genImg_signed(image, image_channels, image_height, image_width, cipherBits);
        genFilters_signed(filters, filter_num, filter_channels, filter_height, filter_width, plaintBits);

        TIMER_CONTINUE(encode)
        auto encodeImg = paillierNN.encode(image,filterSize,stride,padType,slotBits);
        TIMER_PAUSE(encode)

        // cout << "encrypt" << endl;
        TIMER_CONTINUE(encrypt)
        auto imgCt = paillierNN.encrypt(encodeImg,slotBits);
        TIMER_PAUSE(encrypt)

        TIMER_CONTINUE(convolutionFast)
        auto resultCt = paillierNN.convolutionFast(imgCt,filters,plaintBits);
        TIMER_PAUSE(convolutionFast)

        auto matSize = paillierNN.getOutSize_conv(image,filterSize,stride,padType);

        TIMER_CONTINUE(decrypt)
        auto result = paillierNN.decrypt(resultCt,matSize);
        TIMER_PAUSE(decrypt)
    }
    TIMER_PRINT(encode)
    TIMER_PRINT(encrypt)
    TIMER_PRINT(convolutionFast)
    TIMER_PRINT(decrypt)
    // cout << "output" << endl;
    // for (size_t i = 0; i < result.length(); i++)
    // {
    //     for (size_t j = 0; j < result[i].length(); j++)
    //     {
    //         printVector(result[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
    // cout<<endl;
}

void testConvolution(){
    std::cout << "start verify testConvolution:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    PaillierNN::Image image;
    image.channels = 128;
    image.height = 16;
    image.width = 16;
    size_t cipherBits = 8;
    image.data.SetLength(image.channels);
    for (size_t ch = 0; ch < image.channels; ch++)
    {
        image.data[ch].SetLength(image.height);
        for (size_t h = 0; h < image.height; h++)
        {
            image.data[ch][h].SetLength(image.width);
            for (size_t w = 0; w < image.width; w++)
            {
                image.data[ch][h][w] = rand()&((1<<cipherBits)-1);
            }
            // printVector(image.data[ch][h]);
        }
        // std::cout << endl;
    }
    

    PaillierNN::Filters filters;
    filters.SetLength(128);
    size_t plaintBits = 8;
    for (size_t i = 0; i < filters.length(); i++)
    {
        filters[i].channels = image.channels;
        filters[i].height = 3;
        filters[i].width = 3;
        filters[i].data.SetLength(filters[i].channels);
        for (size_t j = 0; j < filters[i].data.length(); j++)
        {
            filters[i].data[j].SetLength(filters[i].height);
            for (size_t k = 0; k < filters[i].data[j].length(); k++)
            {
                filters[i].data[j][k].SetLength(filters[i].width);
                for (size_t t = 0; t < filters[i].data[j][k].length(); t++)
                {
                    filters[i].data[j][k][t] = rand() & ((1<<plaintBits)-1);
                }
                // printVector(filters[i].data[j][k]);
            }
            // cout << endl;
        }
    }
    PaillierNN::FilterSize filterSize;
    filterSize.channels = filters[0].channels;
    filterSize.height = filters[0].height;
    filterSize.width = filters[0].width;
    filterSize.filterNum = filters.length();
    PaillierNN::Stride stride;
    stride.height = 2;
    stride.width = 2;
    PaillierNN::PaddingType padType = PaillierNN::PaddingType::VALID;
    
    size_t slotBits = 32;
    PaillierNN paillierNN(1024);
    // cout << "encode" << endl;

    cout << "Image("<< image.height << "*" << image.width <<"," << image.channels <<"); Filters("<< filters[0].height << "*"<< filters[0].width<< ","  << filters.length()<<")" << endl;
    auto encodeImg = paillierNN.encode(image,filterSize,stride,padType,slotBits);

    // for (size_t i = 0; i < encodeImg.length(); i++)
    // {
    //     printVector(encodeImg[i]);
    // }
    // std::cout << endl;

    // cout << "encrypt" << endl;
    TIMER_START(encrypt)
    auto imgCt = paillierNN.encrypt(encodeImg,slotBits);
    TIMER_STOP(encrypt)

    TIMER_START(convolution)
    auto resultCt = paillierNN.convolution(imgCt,filters);
    TIMER_STOP(convolution)

    auto matSize = paillierNN.getOutSize_conv(image,filterSize,stride,padType);

    TIMER_START(decrypt)
    auto result = paillierNN.decrypt(resultCt,matSize);
    TIMER_STOP(decrypt)

    // cout << "output" << endl;
    // for (size_t i = 0; i < result.length(); i++)
    // {
    //     for (size_t j = 0; j < result[i].length(); j++)
    //     {
    //         printVector(result[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
    cout<<endl;
}

void testConvolution_Binary(){
    std::cout << "start verify testConvolution_Binary:" << endl;
    const size_t key_size = 2048;

    PaillierNN::Image image;
    size_t cipherBits = 8;
    size_t image_channels = imageChannels;
    size_t image_height = imageHeight;
    size_t image_width = imageWidth;

    size_t plaintBits = 1;
    PaillierNN::Filters filters;
    size_t filter_num = filterNum;
    size_t filter_channels = image_channels;
    size_t filter_height = filterHeight;
    size_t filter_width = filterWidth;
    
    PaillierNN::Stride stride;
    stride.height = strideHeight;
    stride.width = strideWidth;

    PaillierNN::FilterSize filterSize;
    filterSize.channels = filter_channels;
    filterSize.height = filter_height;
    filterSize.width = filter_width;
    filterSize.filterNum = filter_num;
    
    PaillierNN::PaddingType padType = PaillierNN::PaddingType::VALID;
    
    size_t slotBits = plaintBits+cipherBits+log(filter_channels * filter_height * filter_width)+1;
    
    PaillierNN paillierNN(key_size);

    // cout << "encode:" << endl;
    std::cout << "keySize=" << key_size <<"; slotBits=" << slotBits << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << endl;
    cout << "Image("<< image_height << "*" << image_width <<"," << image_channels <<"); Filters("<< filter_height << "*"<< filter_width << ","  << filter_num <<")" << endl;
    
    //init TIMER
    TIMER_INIT(encode,TIMES)
    TIMER_INIT(encrypt,TIMES)
    TIMER_INIT(convolution_Binary,TIMES)
    TIMER_INIT(decrypt,TIMES)

    for (size_t i = 0; i < TIMES; i++)
    {
        // gen random data
        genImg_unsigned(image, image_channels, image_height, image_width, cipherBits);
        genFilters_binary(filters, filter_num, filter_channels, filter_height, filter_width);

        TIMER_CONTINUE(encode)
        auto encodeImg = paillierNN.encode(image,filterSize,stride,padType,slotBits);
        TIMER_PAUSE(encode)

        // cout << "encrypt" << endl;
        TIMER_CONTINUE(encrypt)
        auto imgCt = paillierNN.encrypt(encodeImg,slotBits);
        TIMER_PAUSE(encrypt)

        TIMER_CONTINUE(convolution_Binary)
        auto resultCt = paillierNN.convolution_Binary(imgCt,filters);
        TIMER_PAUSE(convolution_Binary)

        auto matSize = paillierNN.getOutSize_conv(image,filterSize,stride,padType);

        TIMER_CONTINUE(decrypt)
        auto result = paillierNN.decrypt(resultCt,matSize);
        TIMER_PAUSE(decrypt)
    }
    TIMER_PRINT(encode)
    TIMER_PRINT(encrypt)
    TIMER_PRINT(convolution_Binary)
    TIMER_PRINT(decrypt)
    // cout << "output" << endl;
    // for (size_t i = 0; i < result.length(); i++)
    // {
    //     for (size_t j = 0; j < result[i].length(); j++)
    //     {
    //         printVector(result[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
    
}

void testConvolution_AbsOne(){
    std::cout << "start verify testConvolution_AbsOne:" << endl;
    const size_t key_size = 2048;

    PaillierNN::Image image;
    size_t cipherBits = 8;
    size_t image_channels = imageChannels;
    size_t image_height = imageHeight;
    size_t image_width = imageWidth;

    size_t plaintBits = 1;
    PaillierNN::Filters filters;
    size_t filter_num = filterNum;
    size_t filter_channels = image_channels;
    size_t filter_height = filterHeight;
    size_t filter_width = filterWidth;
    
    PaillierNN::Stride stride;
    stride.height = strideHeight;
    stride.width = strideWidth;

    PaillierNN::FilterSize filterSize;
    filterSize.channels = filter_channels;
    filterSize.height = filter_height;
    filterSize.width = filter_width;
    filterSize.filterNum = filter_num;
    
    PaillierNN::PaddingType padType = PaillierNN::PaddingType::VALID;
    
    size_t slotBits = plaintBits+cipherBits+log(filter_channels * filter_height * filter_width)+1;
    
    PaillierNN paillierNN(key_size);

    // cout << "encode:" << endl;
    std::cout << "keySize=" << key_size <<"; slotBits=" << slotBits << endl;
    std::cout << "cipherBits=" << cipherBits <<"; plaintBits=" << plaintBits << endl;
    cout << "Image("<< image_height << "*" << image_width <<"," << image_channels <<"); Filters("<< filter_height << "*"<< filter_width << ","  << filter_num <<")" << endl;
    
    //init TIMER
    TIMER_INIT(encode,TIMES)
    TIMER_INIT(encrypt,TIMES)
    TIMER_INIT(convolution_AbsOne,TIMES)
    TIMER_INIT(decrypt,TIMES)

    for (size_t i = 0; i < TIMES; i++)
    {
        // gen random data
        genImg_unsigned(image, image_channels, image_height, image_width, cipherBits);
        genFilters_absOne(filters, filter_num, filter_channels, filter_height, filter_width);

        TIMER_CONTINUE(encode)
        auto encodeImg = paillierNN.encode(image,filterSize,stride,padType,slotBits);
        TIMER_PAUSE(encode)

        // cout << "encrypt" << endl;
        TIMER_CONTINUE(encrypt)
        auto imgCt = paillierNN.encrypt(encodeImg,slotBits);
        TIMER_PAUSE(encrypt)

        TIMER_CONTINUE(convolution_AbsOne)
        auto resultCt = paillierNN.convolution_AbsOne(imgCt,filters);
        TIMER_PAUSE(convolution_AbsOne)

        auto matSize = paillierNN.getOutSize_conv(image,filterSize,stride,padType);

        TIMER_CONTINUE(decrypt)
        auto result = paillierNN.decrypt(resultCt,matSize);
        TIMER_PAUSE(decrypt)
    }
    TIMER_PRINT(encode)
    TIMER_PRINT(encrypt)
    TIMER_PRINT(convolution_AbsOne)
    TIMER_PRINT(decrypt)

    // cout << "output" << endl;
    // for (size_t i = 0; i < result.length(); i++)
    // {
    //     for (size_t j = 0; j < result[i].length(); j++)
    //     {
    //         printVector(result[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
    
}


/*------------- test Large MatMul ------------------*/
void testLargeMatMul_AbsOne(){
    std::cout << "start verify testLargeMatMul_AbsOne:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A;
    Vec<Vec<short>> B;
    int slot_bits = 16;
    int A_rows = 256;
    int A_cols = 256;
    int B_cols = 256;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            A[i][j] = Integer(rand()%256);
            if (A[i][j] == 0)
            {
                A[i][j] = -1;
            }
        }
        printVector(A[i]);
    }
    std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            auto r = rand()%2;
            if (r)
            {
                B[i][j] = 1;
                for (size_t k = 0; k < A_rows; k++)
                {
                    v_result[k][j] += A[k][i];
                }
            }   
            else
            {
                B[i][j] = -1;
                for (size_t k = 0; k < A_rows; k++)
                {
                    v_result[k][j] -= A[k][i];
                }
            }      
        }
        printShort(B[i]);
    }
    std::cout << "verify result:" << endl;
    for (size_t i = 0; i < v_result.length(); i++)
        printVector(v_result[i]);

    std::cout << "check Slots: " << paillierMat.checkSlotSize(slot_bits,8,1,A_cols) << endl;
    TIMER_START(encrypt_largeMat);
    auto cts = paillierMat.encryptTrans_largeMat(A,slot_bits);
    TIMER_STOP(encrypt_largeMat);
    TIMER_START(matrixMulAbsOne_largeMat);
    auto resultCt = paillierMat.matrixMulAbsOne_largeMat(cts,B);
    TIMER_STOP(matrixMulAbsOne_largeMat);
    TIMER_START(decrypt_largeMat);
    auto results = paillierMat.decrypt_largeMat(resultCt);
    TIMER_STOP(decrypt_largeMat);
    // 开始验证
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
                flag = false;
        }
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testLargeMatMul_Binary(){
    std::cout << "start verify testLargeMatMul:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A;
    Vec<Vec<bool>> B;
    int slot_bits = 32;
    int A_rows = 128;
    int A_cols = 128;
    int B_cols = 128;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            A[i][j] = Integer(rand()%256);
        }
        printVector(A[i]);
    }
    std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            B[i][j] = rand()%2;
            if (B[i][j])
            {
                for (size_t k = 0; k < A_rows; k++)
                {
                    v_result[k][j] += A[k][i];
                }
            }         
        }
        printBool(B[i]);
    }
    std::cout << "verify result:" << endl;
    for (size_t i = 0; i < v_result.length(); i++)
        printVector(v_result[i]);

    std::cout << "check Slots: " << paillierMat.checkSlotSize(slot_bits,8,1,A_cols) << endl;
    TIMER_START(encrypt_largeMat);
    auto cts = paillierMat.encryptTrans_largeMat(A,slot_bits);
    TIMER_STOP(encrypt_largeMat);
    TIMER_START(matrixMulBinary_largeMat);
    auto resultCt = paillierMat.matrixMulBinary_largeMat(cts,B);
    TIMER_STOP(matrixMulBinary_largeMat);
    TIMER_START(decrypt_largeMat);
    auto results = paillierMat.decrypt_largeMat(resultCt);
    TIMER_STOP(decrypt_largeMat);
    // 开始验证
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
                flag = false;
        }
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testLargeMatMul(){
    std::cout << "start verify testLargeMatMul:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A,B;
    int A_rows = 3;
    int A_cols = 3;
    int B_cols = 3;
    size_t plaintBits = 4,cipherBits = 2;
    int slot_bits = plaintBits+cipherBits+log(A_cols)+1;//这里+2是因为log(A_cols)需要向上取整,但是log()函数并没有向上取整,因此这里是plaintBits+cipherBits+(log(A_cols)+1)+1;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    // std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    // std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            A[i][j] = Integer(rand()%(1<<cipherBits) - (1<<(cipherBits-1)));
        }
        printVector(A[i]);
    }
    // std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            B[i][j] = Integer(rand()%(1<<plaintBits) - (1<<(plaintBits-1)));
            for (size_t k = 0; k < A_rows; k++)
            {
                v_result[k][j] += A[k][i]*B[i][j];
            }
        }
        printVector(B[i]);
    }
    // std::cout << "verify result:" << endl;
    // for (size_t i = 0; i < v_result.length(); i++)
    //     printVector(v_result[i]);

    std::cout << "keySize=" << key_size <<"; slotBits=" << slot_bits << endl;
    std::cout << "Matrix A:(" << A.length() << "," << A[0].length() \
            << "); Matrix B:(" << B.length() << "," << B[0].length() << ")" << endl;
    std::cout << "check Slots: " << paillierMat.checkSlotSize(slot_bits,cipherBits,plaintBits,A_cols) << endl;
    TIMER_START(encrypt_largeMat);
    auto cts = paillierMat.encryptTrans_largeMat(A,slot_bits);
    TIMER_STOP(encrypt_largeMat);
    TIMER_START(matrixMul_largeMat);
    auto resultCt = paillierMat.matrixMul_largeMat(cts,B);
    TIMER_STOP(matrixMul_largeMat);
    TIMER_START(decrypt_largeMat);
    auto results = paillierMat.decrypt_largeMat(resultCt);
    TIMER_STOP(decrypt_largeMat);
    // 开始验证
    cout<<"开始验证" << endl;
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
                flag = false;
        }
        //  printVector(results[i]);
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
    std::cout << endl;
}

void matrixMulFast_largeMat(){
    std::cout << "start verify matrixMulFast_largeMat:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A,B;
    int A_rows = 3;
    int A_cols = 3;
    int B_cols = 3;
    size_t plaintBits = 1,cipherBits = 2;
    int slot_bits = plaintBits+cipherBits+log2(A_cols)+1;//这里+2是因为log(A_cols)需要向上取整,但是log()函数并没有向上取整,因此这里是plaintBits+cipherBits+(log(A_cols)+1)+1;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    // std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            // A[i][j] = Integer(rand()%(1<<cipherBits) - (1<<(cipherBits-1)));
            A[i][j] = Integer(rand()%2);
            if (A[i][j]==0)
            {
                A[i][j] = -1;
            }
            
        }
        printVector(A[i]);
    }
    std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            B[i][j] = Integer(rand()%(1<<plaintBits) - (1<<(plaintBits-1)));
            // B[i][j] = Integer(rand()%(1<<plaintBits) - (1<<(plaintBits-1)));
            for (size_t k = 0; k < A_rows; k++)
            {
                v_result[k][j] += A[k][i]*B[i][j];
            }
        }
        printVector(B[i]);
    }
    std::cout << "verify result:" << endl;
    for (size_t i = 0; i < v_result.length(); i++)
        printVector(v_result[i]);

    std::cout << "keySize=" << key_size <<"; slotBits=" << slot_bits << endl;
    std::cout << "Matrix A:(" << A.length() << "," << A[0].length() \
            << "); Matrix B:(" << B.length() << "," << B[0].length() << ")" << endl;
    std::cout << "check Slots: " << paillierMat.checkSlotSize(slot_bits,cipherBits,plaintBits,A_cols) << endl;
    TIMER_START(encrypt_largeMat);
    auto cts = paillierMat.encryptTrans_largeMat(A,slot_bits);
    TIMER_STOP(encrypt_largeMat);
    TIMER_START(matrixMulFast_largeMat);
    auto resultCt = paillierMat.matrixMulFast_largeMat(cts,B,plaintBits);
    TIMER_STOP(matrixMulFast_largeMat);
    TIMER_START(decrypt_largeMat);
    auto results = paillierMat.decrypt_largeMat(resultCt);
    TIMER_STOP(decrypt_largeMat);
    // 开始验证
    cout<<"开始验证" << endl;
    // std::cout << "output size: " << results.length() << "*" << results[0].length()<<std::endl;
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
                flag = false;
        }
         printVector(results[i]);
    }

    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

/*------------- test MatMul ------------------*/
void testMatMul_AbsOne(){
    std::cout << "start verify testMatMul_AbsOne :" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A;
    Vec<Vec<short>> B;
    int A_rows = 2;
    int A_cols = 3;
    int B_cols = 3;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            A[i][j] = Integer(rand()%20);
        }
        printVector(A[i]);
    }
    std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            B[i][j] = rand()%2;
            if (B[i][j] == 1)
            {
                B[i][j] = 1;
                for (size_t k = 0; k < A_rows; k++)
                {
                    v_result[k][j] += A[k][i];
                }
            }else
            {
                B[i][j] = -1;
                for (size_t k = 0; k < A_rows; k++)
                {
                    v_result[k][j] -= A[k][i];
                }
            }     
        }
        printShort(B[i]);
    }
    auto cts = paillierMat.encrypt(A,64);
    auto resultCt = paillierMat.matrixMulAbsOne(cts,B);
    auto results = paillierMat.decrypt(resultCt);
    for (size_t i = 0; i < results.length(); i++)
    {
        printVector(results[i]);
    }
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
                flag = false;
        }
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testMatMul_Binary(){
    std::cout << "start verify testMatMul_Binary:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A;
    Vec<Vec<bool>> B;
    int A_rows = 2;
    int A_cols = 3;
    int B_cols = 3;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            A[i][j] = Integer(rand()%20);
        }
        printVector(A[i]);
    }
    std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            B[i][j] = rand()%2;
            if (B[i][j])
            {
                for (size_t k = 0; k < A_rows; k++)
                {
                    v_result[k][j] += A[k][i];
                }
            }     
        }
        printBool(B[i]);
    }
    auto cts = paillierMat.encrypt(A,64);
    auto resultCt = paillierMat.matrixMulBinary(cts,B);
    auto results = paillierMat.decrypt(resultCt);
    for (size_t i = 0; i < results.length(); i++)
    {
        printVector(results[i]);
    }
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
                flag = false;
        }
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testMatMul(){
    std::cout << "start verify matrix Mul `A*B`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierMat paillierMat(key_size);
    srand(time(NULL));
    Vec<Vec<Integer>> A,B;
    int A_rows = 2;
    int A_cols = 3;
    int B_cols = 3;
    A.SetLength(A_rows);
    B.SetLength(A_cols);
    Vec<Vec<Integer>> v_result;
    v_result.SetLength(A_rows);
    for (size_t i = 0; i < A_rows; i++)
    {
        A[i].SetLength(A_cols);
        v_result[i].SetLength(B_cols);
    }
    std::cout << "testMatMul:1" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {
        B[i].SetLength(B_cols);
    }
    std::cout << "matrix A:" << endl;
    for (size_t i = 0; i < A_rows; i++)
    {    
        for (size_t j = 0; j < A_cols; j++)
        {
            A[i][j] = Integer(rand()%20);
        }
        printVector(A[i]);
    }
    std::cout << "matrix B:" << endl;
    for (size_t i = 0; i < A_cols; i++)
    {    
        for (size_t j = 0; j < B_cols; j++)
        {
            B[i][j] = Integer(rand()%20);
            for (size_t k = 0; k < A_rows; k++)
            {
                v_result[k][j] += A[k][i]*B[i][j];
            }
        }
        printVector(B[i]);
    }
    auto cts = paillierMat.encrypt(A,64);
    auto resultCt = paillierMat.matrixMul(cts,B);
    std::cout << "testMatMul:6" << endl;
    auto results = paillierMat.decrypt(resultCt);
    std::cout << "testMatMul:7" << endl;
    for (size_t i = 0; i < results.length(); i++)
    {
        printVector(results[i]);
    }
    bool flag = true;
    for (size_t i = 0; i < results.length(); i++)
    {
        for (size_t j = 0; j < results[i].length(); j++)
        {
            if (results[i][j] != v_result[i][j])
            {
                std::cout << "verify failed." << endl;
                flag = false;
            }
        }
    }
    if (flag)
    {
        std::cout << "verify success." << endl;
    }
}

/*------------- test mat_vec_Mul ------------------*/
void testDP_mat_vec_AbsOne(){
    std::cout << "start testDP_mat_vec_AbsOne `A*b`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierDP paillierDP(key_size);
    srand(time(NULL));
    Vec<short> b;
    Vec<Vec<Integer>> A;
    int rows = 10;
    int cols = 10;
    size_t plaintBits = 8,cipherBits = 8;
    A.SetLength(rows);
    b.SetLength(cols);
    for (size_t i = 0; i < rows; i++)
    {
        A[i].SetLength(cols);
    }
    
    Vec<Integer> v_result;
    v_result.SetLength(rows);
    for (size_t i = 0; i < cols; i++)
    {
        auto r = rand()%2;
        for (size_t j = 0; j < rows; j++)
        {
            A[j][i] = Integer(rand()&((2<<plaintBits)-1));
            if (r == 1)
            {
                b[i] = 1;
                v_result[j] += A[j][i];
            }
            else
            {
                b[i] = -1;
                v_result[j] -= A[j][i];
            }
        }
    }
    std::cout << "random matrix A result:" << endl;
    for (size_t i = 0; i < rows; i++)
    {
        printVector(A[i]);
    }
    std::cout << "random vector B result:" << endl;
    printShort(b);
    auto ct = paillierDP.encrypt_pack(A,64);
    auto dpResult = paillierDP.dotProductAbsOne_pack(ct,b);
    auto result = paillierDP.decrypt_pack(dpResult);
    std::cout << "decrypt result: " << endl;
    printVector(result);
    bool flag = true;
    for (size_t i = 0; i < rows; i++)
    {
        if (v_result[i] != result[i])
        {
            flag = false;
        }
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testDP_mat_vec_Binary(){
    std::cout << "start testDP_mat_vec_Binary `A*b`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierDP paillierDP(key_size);
    srand(time(NULL));
    Vec<bool> b;
    Vec<Vec<Integer>> A;
    int rows = 10;
    int cols = 10;
    size_t plaintBits = 8,cipherBits = 8;
    A.SetLength(rows);
    b.SetLength(cols);
    for (size_t i = 0; i < rows; i++)
    {
        A[i].SetLength(cols);
    }
    
    Vec<Integer> v_result;
    v_result.SetLength(rows);
    for (size_t i = 0; i < cols; i++)
    {
        b[i] = rand()%2;
        for (size_t j = 0; j < rows; j++)
        {
            A[j][i] = Integer(rand()&((2<<plaintBits)-1));
            if (b[i])
            {
                v_result[j] += A[j][i];
            }
        }
    }
    std::cout << "random matrix A result:" << endl;
    for (size_t i = 0; i < rows; i++)
    {
        printVector(A[i]);
    }
    std::cout << "random vector B result:" << endl;
    printBool(b);
    auto ct = paillierDP.encrypt_pack(A,64);
    auto dpResult = paillierDP.dotProductBinary_pack(ct,b);
    auto result = paillierDP.decrypt_pack(dpResult);
    std::cout << "decrypt result: " << endl;
    printVector(result);
    bool flag = true;
    for (size_t i = 0; i < rows; i++)
    {
        if (v_result[i] != result[i])
        {
            flag = false;
        }
    }
    if (flag)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testDP_mat_vec(){
    std::cout << "start testDP_mat_vec `A*b`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierDP paillierDP(key_size);
    srand(time(NULL));
    Vec<Integer> b;
    Vec<Vec<Integer>> A;
    int rows = 10;
    int cols = 10;
    A.SetLength(rows);
    b.SetLength(cols);
    for (size_t i = 0; i < rows; i++)
    {
        A[i].SetLength(cols);
    }
    
    Vec<Integer> v_result;
    v_result.SetLength(rows);
    for (size_t i = 0; i < cols; i++)
    {
        b[i] = Integer(rand()%20);
        
        for (size_t j = 0; j < rows; j++)
        {
            A[j][i] = Integer(rand()%20);
            v_result[j] += A[j][i]* b[i];
        }
    }
    std::cout << "random matrix A result:" << endl;
    for (size_t i = 0; i < rows; i++)
    {
        printVector(A[i]);
    }
    std::cout << "random vector B result:" << endl;
    printVector(b);
    auto ct = paillierDP.encrypt_pack(A,64);
    auto dpResult = paillierDP.dotProduct_pack(ct,b);
    auto result = paillierDP.decrypt_pack(dpResult);
    std::cout << "decrypt result: " << endl;
    printVector(result);
    for (size_t i = 0; i < rows; i++)
    {
        if (v_result[i] != result[i])
        {
            std::cout << "verify failed." << endl;
        }
    }
    std::cout << "verify success." << endl;
    
}

/*------------- test vec_vec_Mul ------------------*/
void testDP_vec_vec_AbsOne(){
    std::cout << "start verify testDP_vec_vec_AbsOne `a*b`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierDP paillierDP(key_size);
    srand(time(NULL));
    Vec<Integer> a;
    Vec<short> b;
    int length = 5;
    a.SetLength(length);
    b.SetLength(length);
    Integer v_result = 0;
    for (size_t i = 0; i < length; i++)
    {
        a[i] = Integer(rand()%20);
        auto r = rand()%2;
        if (r)
        {
            b[i] = 1;
            v_result += a[i];
        }
        else
        {
            b[i] = -1;
            v_result -= a[i];
        }
    }
    std::cout << "random A result:" << endl;
    printVector(a);
    std::cout << "random B result:" << endl;
    printShort(b);
    auto ct = paillierDP.encrypt(a,64);
    auto dpResult = paillierDP.dotProductAbsOne(ct,b);
    auto result = paillierDP.decrypt(dpResult,b.length());
    std::cout << "decrypt result: " << result << endl;
    if (v_result == result)
        std::cout << "verify success." << endl;
    else
        std::cout << "verify failed." << endl;
}

void testDP_vec_vec_Binary(){
    std::cout << "start verify testDP_vec_vec_Binary `a*b`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierDP paillierDP(key_size);
    srand(time(NULL));
    Vec<Integer> a;
    Vec<bool> b;
    int length = 5;
    a.SetLength(length);
    b.SetLength(length);
    Integer v_result = 0;
    for (size_t i = 0; i < length; i++)
    {
        a[i] = Integer(rand()%20);
        b[i] = rand()%2;
        if (b[i])
        {
            v_result += a[i];
        }
    }
    std::cout << "random A result:" << endl;
    printVector(a);
    std::cout << "random B result:" << endl;
    printBool(b);
    auto ct = paillierDP.encrypt(a,64);
    auto dpResult = paillierDP.dotProductBinary(ct,b);
    auto result = paillierDP.decrypt(dpResult,b.length());
    std::cout << "decrypt result: " << result << endl;
    if (v_result == result)
    {
        std::cout << "verify success." << endl;
    }
}

void testDP_vec_vec(){
    std::cout << "start verify vector `a*b`:" << endl;
    const size_t key_size = 1024;
    opaillierlib::PaillierDP paillierDP(key_size);
    srand(time(NULL));
    Vec<Integer> a,b;
    int length = 5;
    a.SetLength(length);
    b.SetLength(length);
    Integer v_result = 0;
    for (size_t i = 0; i < length; i++)
    {
        a[i] = Integer(rand()%20);
        b[i] = Integer(rand()%20);
        v_result += a[i]* b[i];
    }
    std::cout << "random A result:" << endl;
    printVector(a);
    std::cout << "random B result:" << endl;
    printVector(b);
    auto ct = paillierDP.encrypt(a,64);
    auto dpResult = paillierDP.dotProduct(ct,b);
    auto result = paillierDP.decrypt(dpResult,0);
    std::cout << "decrypt result: " << result << "; v_result=" << v_result << endl;
    if (v_result == result)
    {
        std::cout << "verify success." << endl;
    }
}

/*------------- test NewPaillier ------------------*/
void testNewPaillier(){
    const size_t key_size = 2048;
    opaillierlib::NewPaillier paillier(key_size);
    paillier.generate_keys();
    Vec<Integer> a;
    size_t slotBit = 16;
    size_t length = key_size/slotBit;
    genVector_signed(a, length, 8);

    TIMER_START(encrypt_pack)
    paillier.encrypt_pack(a.begin(),a.end(),slotBit);
    TIMER_STOP(encrypt_pack)
    TIMER_START(encrypt)
    paillier.encrypt(a[0]);
    TIMER_STOP(encrypt)
    // srand(time(NULL));
    // Vec<Integer> a,b;
    // int length = 2;
    // a.SetLength(length);
    // b.SetLength(length);
    // for (size_t i = 0; i < length; i++)
    // {
    //     a[i] = Integer(rand()%20);
    //     b[i] = Integer(rand()%20);
    // }
    // std::cout << "random A result:" << endl;
    // printVector(a);
    // std::cout << "random B result:" << endl;
    // printVector(b);
    // auto ct = paillier.encrypt_pack(a,64);
    // PackedCiphertext tmp(ct);
    // tmp.data *= Integer(1)<<(tmp.plaintext_bits*1);//移位;

    // Integer base =1;
    // size_t shift = 1;
    // Ciphertext tmp_mul;
    // ophelib::PackedCiphertext result = ophelib::PackedCiphertext(ct.data * b[0], ct.n_plaintexts, ct.plaintext_bits);
    // for (size_t i = 1; i < b.length(); i++)
    // {
    //     tmp_mul = ct.data * b[i];
    //     tmp_mul *= (base<<(ct.plaintext_bits*shift));//移位
    //     result.data += tmp_mul;
    //     shift++;
    // }
    // std::cout << "mul result:" << endl;
    // auto result_tmp = paillier.decrypt_pack(result);
    // printVector(result_tmp);
}

/*------------- print ------------------*/
template <typename T>
void printVector(const Vec<T>& values){
    for (size_t i = 0; i < values.length(); i++)
    {
        std::cout << values[i] << " ";
    }
    std::cout << endl;
}

void printBool(const Vec<bool>& values){
    for (size_t i = 0; i < values.length(); i++)
    {
        std::cout << values[i] << " ";
    }
    std::cout << endl;
}

void printShort(const Vec<short>& values){
    for (size_t i = 0; i < values.length(); i++)
    {
        std::cout << values[i] << " ";
    }
    std::cout << endl;
}


