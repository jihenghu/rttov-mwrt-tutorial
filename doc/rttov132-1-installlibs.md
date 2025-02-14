
# RTTOV模式笔记：(一) 依赖安装

本节介绍RTTOV的服务器安装，除必要的软件支持，都使用普通用户权限来完成，做个记录；
感谢[大气快速辐射传输模型RTTOV12.2安装教程及心得体会](https://blog.csdn.net/weixin_43471242/article/details/103248318)的详细说明，这里记录一下自己的安装过程。

### 安装环境和安装包准备 

环境：
- Ubuntu 20.04.5 LTS (GNU/Linux 5.15.0-75-generic x86_64)
- gcc version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.1)
	- gfortran 9.4.0 
	- g++ 9.4.0
- GNU Make 4.2.1

安装包：
- zlib-1.2.11.tar.gz
- hdf5-1.8.21.tar.gz
- netcdf-c-4.9.2.tar.gz
- netcdf-fortran-4.6.1.tar.gz
- RTTOV 13.2

除了RTTOV, 其他的安装包你可以在[Compiling WRF](https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php) 或者官方网站下载

安装顺序就按照以上顺序，否则会出现依赖问题。
将上述安装包放置到/home/hjh/下，并分别新建目录，结构如下：
```bash
root@node05:/home/hjh$ tree -L 1
.
├── hdf5
├── netcdf
├── zlib
├── RTTOV13
├── hdf5-1.8.21.tar.gz
├── netcdf-c-4.9.2.tar.gz
├── netcdf-fortran-4.6.1.tar.gz
└── zlib-1.2.11.tar.gz
```
#### zlib
```bash
~$ tar -zvxf zlib-1.2.11.tar.gz
~$ cd zlib-1.2.11
~/zlib-1.2.11$ ./configure --prefix=/home/hjh/zlib
~/zlib-1.2.11$ make
~/zlib-1.2.11$ make check
~/zlib-1.2.11$ make install
$ cd
$ rm -r zlib-1.2.11
```
从make开始一般比较顺利，没有报错信息；
netcdf的安装需要依赖zlib,先在/home/hjh/.bashrc定义zlib路径，以配置netcdf的编译信息：

```bash
~$ vim ~/.bashrc 
## 追加已下两行
## zlib
export LD_LIBRARY_PATH=/home/hjh/zlib/lib:$LD_LIBRARY_PATH

## 保存退出
~$ source ~/.bashrc
~$ echo $LD_LIBRARY_PATH
```


安装zlib*，这里需要联系管理员（权限狗:P）
安装时避免目录下有zlib开头的文件，否则会导致正则匹配并提示找不到安装包；

```bash
(一个不含有zlib*的文件夹)# sudo apt install zlib*
```

#### HDF5

```bash
$ cd
~ $ tar zxvf hdf5-1.8.21.tar.gz
~ $ cd hdf5-1.8.21
~/hdf5-1.8.21 $ ./configure --with-zlib=/home/hjh/zlib --prefix=/home/hjh/hdf5 FC=gfortran CC=gcc --enable-fortran --enable-cxx
```
配置完以后会出现配置总结：
```
	     Installation point: /home/hjh/hdf5
					......
         C Compiler: /usr/bin/gcc
					......
         Fortran: yes
			Fortran Compiler: /usr/bin/gfortran
					......
         C++: yes
            C++ Compiler: /usr/bin/g++
					......
         Parallel HDF5: no
					......
```

编译安装，make 过程非常长，warning可以忽略
```bash
$ make 
$ make check
$ make install
```
可以在安装目录下/home/hjh/hdf5/bin 发现可执行文件；/home/hjh/hdf5/lib下发现库文件

```bash
~$ vim ~/.bashrc 
## 追加以下内容
## HDF5 1.8.21
export PATH=/home/hjh/hdf5/bin:$PATH
export LD_LIBRARY_PATH=/home/hjh/hdf5/lib:$LD_LIBRARY_PATH
export CPPFLAGS=-I/home/hjh/hdf5/include  ## comment this after compiling
export LDFLAGS=-L/home/hjh/hdf5/lib       #
## 保存退出
```
```bash
~$ source ~/.bashrc
~$ echo $LD_LIBRARY_PATH
~$ echo $CPPFLAGS
~$ echo $LDFLAGS

$which h5fc 
/home/hjh/hdf5/bin/h5fc

```

#### netcdf 

```bash
$ cd
~$ tar zxvf netcdf-c-4.9.2.tar.gz
~$ cd netcdf-c-4.9.2
netcdf-c-4.9.2 $ 
```

联系管理员安装
```bash
#apt install m4
#apt install libcurl4-openssl-dev
```
编译安装
```bash
netcdf-c-4.9.2 $ ./configure --prefix=/home/hjh/netcdf --enable-netcdf-4 --disable-libxml2
netcdf-c-4.9.2 $ make
netcdf-c-4.9.2 $ make check
netcdf-c-4.9.2 $ make install
```

安装完成发现/home/hjh/netcdf/bin下有很多可执行文件，但是/home/hjh/netcdf/lib 没有fortran相关的库文件，所以后面需要安装netcdf-fortran

```bash
~$ vim ~/.bashrc 

## HDF5 1.8.21
export PATH=/home/hjh/hdf5/bin:$PATH
export LD_LIBRARY_PATH=/home/hjh/hdf5/lib:$LD_LIBRARY_PATH
## 注释以下内容
#export CPPFLAGS=-I/home/hjh/hdf5/include  
#export LDFLAGS=-L/home/hjh/hdf5/lib       

##追加以下内容

## NETCDF
export PATH=/home/hjh/netcdf/bin:$PATH

export LD_LIBRARY_PATH=/home/hjh/netcdf/lib:$LD_LIBRARY_PATH
export CPPFLAGS=-I/home/hjh/netcdf/include  ## comment this after compiling
export LDFLAGS=-L/home/hjh/netcdf/lib       ## 编译完成后注释掉
## 保存退出
```
```bash
~$ source ~/.bashrc
~$ echo $LD_LIBRARY_PATH
~$ echo $CPPFLAGS
~$ echo $LDFLAGS

$which nccopy
/home/hjh/netcdf/bin/nccopy
```

#### NETCDF-Fortran

```bash
 $ cd 
~ $ tar xzvf netcdf-fortran-4.6.1.tar.gz
~ $ cd netcdf-fortran-4.6.1
```
编译安装
```bash
netcdf-fortran-4.6.1 $ ./configure --prefix=/home/hjh/netcdf FC=gfortran
netcdf-fortran-4.6.1 $ make
netcdf-fortran-4.6.1 $ make check
netcdf-fortran-4.6.1 $ make install
```
安装完成发现/home/hjh/netcdf/lib下新增了fortran相关的库文件

至此，安装RTTOV的依赖已经安装好。

### 参考
[大气快速辐射传输模型RTTOV12.2安装教程及心得体会](https://blog.csdn.net/weixin_43471242/article/details/103248318)
[Compiling WRF](https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php)
