# Tutorial of microwave radiative transfer in RTTOV v13.2
Tutorial proj for microwave radiative simulation based on RTTOV v13.2

## Note
<p style="color:red"> 今天RTTOV更新到了V14版本了，相比V13.2最大的改动是deplicate了coefficient table，在V15版本将完全去除（预计2027年）， 但是会更加好用--- 2025-2-12</p>

## Preface
本文原发布于我的主页，转载于此，欢迎大家前来玩耍 🐱 -> [jihenghu.github.io](https://jihenghu.github.io)

## Chapters
- [(一) 依赖安装](./doc/rttov132-1-installlibs.md) -- 介绍RTTOV的服务器安装所需要的依赖，包括NETCDF,HDF5,zlib等库的安装;
- [(二) RTTOV V13.2安装](./doc/rttov132-2-install.md) -- 介绍RTTOV的服务器安装和消光散射系数的下载，及模式介绍;
- [(三) 约定和特性](./doc/rttov132-3-conventions.md) -- 粗浅解读一下RTTOV的一些基础设定。内容包括廓线插值等;
- [(四) RTTOV 变量和结构体](./doc/rttov132-4-variables.md) -- 解读一下RTTOV的源码，介绍一下RTTOV变量和结构体;
- [(五) 基于Direct Forward的晴空模拟](./doc/rttov132-5-direct-fw.md) -- 晴空正向微波辐射传输的源码结构及模拟实践;
- [(六) 基于MW-SCAT的水凝物模拟](./doc/rttov132-6-mw-scat.md) -- RTTOV MW-SCATT模块对于水凝物散射情形的模拟;
- [(七) 全天气地表微波比辐射反演方案](./doc/rttov132-7-emissivity-retrieve.md) -- 介绍( emissivity反演算法Baordo and Geer, 2016)和相关实践;
- [番外 TELSEM Emissivity Atlas的离线使用](./doc/rttov132-8-telsem2-atlas.md)。


## RTTOV简介 (个人理解)
RTTOV (Radiative Transfer for TOVS)具备对被动可见光、红外和微波卫星辐射计、光谱仪和干涉仪的卫星观测进行快速辐射传输模拟的能力。
也就是说，给定大气温度、湿度以及可选的微量气体、气溶胶和水凝物的大气廓线，再加上表面参数和星载仪器的观测几何参数，RTTOV可以实现对多种卫星任务搭载的传感器进行天顶辐亮度的计算。
其快速辐射传输计算能力得益于内置的多种类型的大气介质在不同波长下的吸收和散射系数的基于机器学习方法的参数表。
不仅如此，RTTOV还提供了对有雅可比矩阵的计算，其描述了天顶亮温对任一环境因子的敏感性，可以应用于多种变分反演算法和数值天气预报模型中；
RTTOV内置了诸多地表发射率、反射率atlas。
- 引用RTTOV
> Saunders, R., Hocking, J., Turner, E., Rayer, P., Rundle, D., Brunel, P., Vidot, J., Roquet, P., Matricardi, M., Geer, A., Bormann, N., and Lupu, C., 2018: An update on the RTTOV fast radiative transfer model (currently at version 12), Geosci. Model Dev., 11, 2717-2737, https://doi.org/10.5194/gmd-11-2717-2018.
- RTTOV主页地址
	https://www.nwpsaf.eu/site/software/rttov/ 。


## 更新进度
- 2023-08-31 
完成模式的安装，更新了前两章，学习源码；
- 2023-08-06 
目前已经初步完成了前四章节的内容，后续会更出剩下的章节，取决于科研进度和这个项目的推进速度，因为不是主线任务，存在更新较慢的可能；
第三节基于user guide的解读还会继续更新；
- 2023-08-10
更新第五章，部分有待补充；
- 2023-08-18
更新第六章，部分有待补充；
- 2023-08-24
更新第七章，部分有待补充；
后续更新随缘了。

## Todo list
- RTTOV Jacobian matrix and variational inversion methods.
- update ERA5 CDSAPI portal.
- update RTTOV to V14.0


## License
This project is licensed under the [BSD-3-Clause license](./LICENSE).