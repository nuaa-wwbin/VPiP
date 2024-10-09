
#include "myRandom.h"
#include <time.h>
#include <random>

void initTime()
{
    srand(time(NULL));
}

/*********** genVector APIs ***************/

void genVector_zero(Vec<Integer> &vec, size_t size){
    // srand(time(NULL));
    vec.SetLength(size);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = 0;
    }
}

void genVector_one(Vec<Integer> &vec, size_t size){
    // srand(time(NULL));
    vec.SetLength(size);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = 1;
    }
}

void genVector_binary(Vec<bool> &vec, size_t size){
    // srand(time(NULL));
    vec.SetLength(size);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = rand()%2;
    }
}

void genVector_absOne(Vec<short> &vec, size_t size){
    // srand(time(NULL));
    vec.SetLength(size);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = rand()%2;
        if (vec[i] == 0)
        {
            vec[i] = -1;
        }
    }
}

void genVector_signed(Vec<Integer> &vec, size_t size, size_t bits){
    // srand(time(NULL));
    vec.SetLength(size);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = (rand() % (1<<bits)) - (1<<(bits-1));
    }
}

void genVector_unsigned(Vec<Integer> &vec, size_t size, size_t bits){
    // srand(time(NULL));
    vec.SetLength(size);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = rand()%(1<<bits);
    }
}

/***************** gen Matrix APIs ***********************/

void genMatrix_binary(Vec<Vec<Integer>> &mat, size_t height, size_t width){
    // srand(time(NULL));
    mat.SetLength(height);
    for (size_t h = 0; h < height; h++)
    {
        mat[h].SetLength(width);
        for (size_t w = 0; w < width; w++)
        {
            mat[h][w] = rand()%2;
        }
    }
}

void genMatrix_absOne(Vec<Vec<Integer>> &mat, size_t height, size_t width){
    // srand(time(NULL));
    mat.SetLength(height);
    for (size_t h = 0; h < height; h++)
    {
        mat[h].SetLength(width);
        for (size_t w = 0; w < width; w++)
        {
            mat[h][w] = rand()%2;
            if (mat[h][w] == 0)
            {
                mat[h][w] = -1;
            }
        }
    }
}

void genMatrix_signed(Vec<Vec<Integer>> &mat, size_t height, size_t width, size_t bits){
    // srand(time(NULL));
    mat.SetLength(height);
    for (size_t h = 0; h < height; h++)
    {
        mat[h].SetLength(width);
        for (size_t w = 0; w < width; w++)
        {
            mat[h][w] = (rand()%(1<<bits)) - (1<<(bits-1));
        }
    }
}

void genMatrix_unsigned(Vec<Vec<Integer>> &mat, size_t height, size_t width, size_t bits){
    // srand(time(NULL));
    mat.SetLength(height);
    for (size_t h = 0; h < height; h++)
    {
        mat[h].SetLength(width);
        for (size_t w = 0; w < width; w++)
        {
            mat[h][w] = rand()%(1<<bits);
        }
    }
}

/*********** genImg APIs ***************/
void genImg_binary(PaillierNN::Image &image, size_t channels, size_t height, size_t width){
    // srand(time(NULL));
    image.channels = channels;
    image.height = height;
    image.width = width;
    image.data.SetLength(image.channels);
    for (size_t ch = 0; ch < image.channels; ch++)
    {
        image.data[ch].SetLength(image.height);
        for (size_t h = 0; h < image.height; h++)
        {
            image.data[ch][h].SetLength(image.width);
            for (size_t w = 0; w < image.width; w++)
            {
                image.data[ch][h][w] = rand()%2;
            }
        }
    }
}

void genImg_absOne(PaillierNN::Image &image, size_t channels, size_t height, size_t width){
    // srand(time(NULL));
    image.channels = channels;
    image.height = height;
    image.width = width;
    image.data.SetLength(image.channels);
    for (size_t ch = 0; ch < image.channels; ch++)
    {
        image.data[ch].SetLength(image.height);
        for (size_t h = 0; h < image.height; h++)
        {
            image.data[ch][h].SetLength(image.width);
            for (size_t w = 0; w < image.width; w++)
            {
                image.data[ch][h][w] = rand()%2;
                if (image.data[ch][h][w] == 0)
                {
                    image.data[ch][h][w] = -1;
                }
            }
        }
    }
}

void genImg_unsigned(PaillierNN::Image &image, size_t channels, size_t height, size_t width, size_t bits){
    // srand(time(NULL));
    image.channels = channels;
    image.height = height;
    image.width = width;
    image.data.SetLength(image.channels);
    for (size_t ch = 0; ch < image.channels; ch++)
    {
        image.data[ch].SetLength(image.height);
        for (size_t h = 0; h < image.height; h++)
        {
            image.data[ch][h].SetLength(image.width);
            for (size_t w = 0; w < image.width; w++)
            {
                image.data[ch][h][w] = rand()%(1<<bits);
            }
        }
    }
}

void genImg_signed(PaillierNN::Image &image, size_t channels, size_t height, size_t width, size_t bits){
    // srand(time(NULL));
    image.channels = channels;
    image.height = height;
    image.width = width;
    image.data.SetLength(image.channels);
    for (size_t ch = 0; ch < image.channels; ch++)
    {
        image.data[ch].SetLength(image.height);
        for (size_t h = 0; h < image.height; h++)
        {
            image.data[ch][h].SetLength(image.width);
            for (size_t w = 0; w < image.width; w++)
            {
                image.data[ch][h][w] = (rand()%(1<<bits)) - (1<<(bits-1));
            }
        }
    }
}


void genFilters_binary(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width){
    filters.SetLength(filterNum);
    // srand(time(NULL));
    for (size_t i = 0; i < filterNum; i++)
    {
        filters[i].channels = channels;
        filters[i].height = height;
        filters[i].width = width;
        filters[i].data.SetLength(filters[i].channels);
        for (size_t j = 0; j < filters[i].data.length(); j++)
        {
            filters[i].data[j].SetLength(filters[i].height);
            for (size_t k = 0; k < filters[i].data[j].length(); k++)
            {
                filters[i].data[j][k].SetLength(filters[i].width);
                for (size_t t = 0; t < filters[i].data[j][k].length(); t++)
                {
                    filters[i].data[j][k][t] = rand()%2;
                }
            }
        }
    }
}

void genFilters_absOne(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width){
    filters.SetLength(filterNum);
    // srand(time(NULL));
    for (size_t i = 0; i < filterNum; i++)
    {
        filters[i].channels = channels;
        filters[i].height = height;
        filters[i].width = width;
        filters[i].data.SetLength(filters[i].channels);
        for (size_t j = 0; j < filters[i].data.length(); j++)
        {
            filters[i].data[j].SetLength(filters[i].height);
            for (size_t k = 0; k < filters[i].data[j].length(); k++)
            {
                filters[i].data[j][k].SetLength(filters[i].width);
                for (size_t t = 0; t < filters[i].data[j][k].length(); t++)
                {
                    filters[i].data[j][k][t] = rand()%2;
                    if (filters[i].data[j][k][t] == 0)
                    {
                        filters[i].data[j][k][t] = -1;
                    }
                }
            }
        }
    }
}

void genFilters_unsigned(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width, size_t bits){
    filters.SetLength(filterNum);
    // srand(time(NULL));
    for (size_t i = 0; i < filterNum; i++)
    {
        filters[i].channels = channels;
        filters[i].height = height;
        filters[i].width = width;
        filters[i].data.SetLength(filters[i].channels);
        for (size_t j = 0; j < filters[i].data.length(); j++)
        {
            filters[i].data[j].SetLength(filters[i].height);
            for (size_t k = 0; k < filters[i].data[j].length(); k++)
            {
                filters[i].data[j][k].SetLength(filters[i].width);
                for (size_t t = 0; t < filters[i].data[j][k].length(); t++)
                {
                    filters[i].data[j][k][t] = rand()%(1<<bits);
                }
            }
        }
    }
}

void genFilters_signed(PaillierNN::Filters &filters, size_t filterNum,size_t channels, size_t height, size_t width, size_t bits){
    filters.SetLength(filterNum);
    // srand(time(NULL));
    for (size_t i = 0; i < filterNum; i++)
    {
        filters[i].channels = channels;
        filters[i].height = height;
        filters[i].width = width;
        filters[i].data.SetLength(filters[i].channels);
        for (size_t j = 0; j < filters[i].data.length(); j++)
        {
            filters[i].data[j].SetLength(filters[i].height);
            for (size_t k = 0; k < filters[i].data[j].length(); k++)
            {
                filters[i].data[j][k].SetLength(filters[i].width);
                for (size_t t = 0; t < filters[i].data[j][k].length(); t++)
                {
                    filters[i].data[j][k][t] = (rand()%(1<<bits)) - (1<<(bits-1));
                }
            }
        }
    }
}