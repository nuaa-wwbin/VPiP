# 项目简介(Introduction)
本项目为基于Paillier库(ophelib)设计的paillier打包下点积,矩阵乘法,以及卷积操作的实现.
This project is an implementation of dot product, matrix multiplication, and convolution operations under the paillier package designed based on ophelib. Please refer to the paper “VPiP: Values Packing in Paillier for Communication Efficient Oblivious Linear Computations” for details.

本项目分为4部分:
- paillier: 一个ophelib::paillierFast的派生类,是所有计算操作的基础
- dotProduct: 用于计算点积操作的类:包含`向量*向量`和`矩阵*向量`;
- matrixMul: 用于计算矩阵乘法的类: 包含`一般矩阵乘法`和`大矩阵乘法`;
- convolution: 用于计算卷积操作的类.

Note: `点积`和`矩阵乘法`都包含3种针对不同明文数据类型的优化实现:
    - 明文为:`实数`,明文时正常的实数,这是一个通用型实现
    - 明文为:`Binary`,明文值为0/1的bool型,这是一个最优实现,只有大整数的加法
    - 明文为:`AbsOne`,明文绝对值为1(即值为-1或1)的short型,这是一个最优实现,只有大整数的加减法

## 依赖
- g++ version: 9.4 (原则上,更早的版本应该也没问题,但是本项目源码只在G++ V9.4测试过)
- cmake version: V3.13
- OPHELib version: 0.3.4 (git url: https://github.com/abb-iss/ophelib.git)


## 准备工作
以下为linux环境下的指导

### 安装g++ 和Cmake
编译工具的安装请直接google/百度,或者参考博客:https://blog.csdn.net/yvhqbat/article/details/50853196

### 安装OPHELib
OPHELib库的按照github上OPHELib项目的build文档直接安装. 

## 运行
在完成准备工作后,执行以下命令

```
git clone http://code.oppoer.me/S9047741/paillier_atomic.git
cd paillier_atomic
mkdir build && cd build
cmake ..
make -j8
./bin/paillier_atomic
```

# 联系方式
如有问题,请直接联系作者:吴伟彬,邮箱(wuweibin@nuaa.edu.cn或者nuaawuweibin@gmail.com)
