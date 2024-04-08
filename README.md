# iPSM
source codes for employing iPSM (individual Predictive Soil Mapping)
## 1.简介
iPSM是基于地理学第三定律的土壤属性空间推测方法
## 2. 文件说明
- doc下面包含了Google Code style 中文帮助、一张图读懂Google Code style，之前组里培训编写的MPI学习文档

- cmake：自动创建C++项目的

- code：实现iPSM功能的C++代码

- solim 下包含了重构代码，包括两套基础代码（parallel-base, block-base)，一套功能实现代码(apps)，以及代码实现过程中用到的第三方代码 (3rdparty)
    基础代码主要负责数据的读写功能，根据读写要求的不同，分为**并行基础代码（parallel-base）**，和**分块读取基础代码（block-base）**两组代码，现对基础代码进行修改，使其为功能型代码提供统一接口，方便代码的修改和维护。
        - parallel-base 用于实现MPI并行读写，用于解决计算效率低下的问题
        - block-base 用于实现分块读写，用于解决处理较大文件时内存不足的问题
    功能实现代码包含在 'apps' 文件夹中，该代码包含了实现iPSM预测制图的核心功能，

## 3. 工程配置与编译
### 3.1 新建iPSM工程

**（1） 安装CMake并配置环境变量**

我们下载和测试的CMake版本为```cmake-3.10.0-win64-x64.msi```。在安装过程中除了需要选择相应的安装路径，还需选择自动配置系统环境变量，即选择```Add CMake to the system PATH for all users```。

CMake的安装具体步骤可详见：https://jingyan.baidu.com/article/da1091fb645ab4027849d6bc.html

**（2）安装MPI并配置环境变量**

在Github中，我们已经提供了MPI相应的下载，在【library】文件夹中，安装【V8】文件夹下的两个SDK程序```msmpisdk.msi```和```MSMpiSetup.exe```。安装完成后配置MPI的环境变量，在MPI的安装目录下，复制【bin】文件夹目录字符串，如```C:\Program Files (x86)\MPI\Bin```，将其复制到系统环境变量中，完成MPI配置。

**（3）配置GDAL**

推荐使用已编译的C++版本，各个编译版本可从以下网站获取：http://www.gisinternals.com/release.php 。 我们选择的GDAL版本为```GDAL 3.4 64位```。配置环境变量,在已编译后的GDAL文件夹中复制【bin】文件夹路径至系统的环境变量中。

**（4） 从Github上Clone或Download项目**

**（5）在iPSM项目文件夹根目录新建【build】文件夹**

**（6）利用VS开发人员命令提示工具，定位到build文件夹**

**（7）进行项目工程生成** 输入```cmake ..```。

**（8）生成工程** 生成之后可在【build】文件夹中打开iPSM的VS工程

**特别注意：**
1. CMake、MPI和GDAL所用版本的位数需要一致，例如都采用64位。
2. 在执行CMake反向工程时，默认采用32位，和之前三个类库保持一致，如果之前的三个类库是64位，同时还需考虑到Visual Studio的版本问题，例如```cmake -G "Visual Studio 12 2013 Win64" ..```
3. 根据开发目的不同，可以在生成工程（步骤7）前对基础代码进行选择：
    运行功能实现代码（即文件夹apps下的代码）时只需要一套基础代码，可通过修改codes文件夹下的CMakeLists.txt 中第二行来选择你所需要的基础代码：

    - 选择parallel-base作为基础代码：

    ```cmake
    SET(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/parallel-base)
    ```

    - 选择block-base作为基础代码：

    ```cmake
    SET(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/block-base)
    ```


### 3.2 在已有工程中新建项目
由于我们采用CMake对iPSM进行反向工程，因此当我们需要在已有项目新建项目或类库扩展功能时，需要特别注意工作流程。

1. 假设此处以Mapping为例，我们需新建Mapping项目，找到需要新建项目的根目录，例如在【apps】下新建文件夹【Mapping】。（**特别注意**，这里的文件夹根目录不是build中的，而是原始代码文件）。
2. 在【Mapping】文件夹中创建```*.cpp```文件和```Readme.txt```文件。（**特别注意**，由于我们之前很多人同学都采用中文写注释，这里新建的```*.cpp```文件需要在创建的时候将编码改为```UTF-8```）
3. 修改CMake配置文件。打开【apps】文件夹下的【CMakeLists.txt】
4. 模仿其他类库和项目，修改相关的配置信息，一下代码是我们已经存在的项目，只需复制并修改代码中mapping或MAPPING中的部分即可，代码如下：
```  
# <mapping>
SET(<MAPPING>_BASED_DIR .</apps/mapping>)
file(GLOB_RECURSE FREQUENCY_BASED_SOURCES ${<MAPPING>_BASED_DIR}/*.cpp)
file(GLOB_RECURSE <MAPPING>_BASED_HEADERS ${<MAPPING>_BASED_DIR}/*.h)
source_group("Header Files" FILES ${<MAPPING>_BASED_HEADERS})
ADD_EXECUTABLE(<mapping> ${<MAPPING>_BASED_SOURCES} ${<MAPPING>_BASED_HEADERS})
```
5. 添加在以下代码段中的<>中添加配置信息：
```
SET(iPSM_TARGETS
				  #demo
				  ipsm
				  ipsmNeighbor
                  mapping
                  )
```
6. 在【build】文件中重新生成iPSM工程文件。


### 3.3 编译程序
1. 使用visual studio打开【build】文件夹中新生成的```*.sln```工程文件

2. 在【解决方案资源管理器中】展开apps文件夹，选择需要编译的项目，如ipsm

3. 右键目标项目（如ipsm），选择【属性】，打开属性界面，将【配置】修改为```Release```，在【配置属性->常规】中，将将输出目录修改为【bin】文件夹的路径

4. 在主界面中【解决方案配置】修改为```Release```，右键目标项目，选择【生成】，生成成功后即可在bin文件夹中找到编译好的ipsm.exe文件

### 3.4 iPSM测试
#### 3.4.1 iPSM参数设置
iPSM项目在进行利用各类算法进行计算时，由于需要使用大量环境协变量，涉及到很多图层，因此，参数的设置采用特殊符号分段的方法。



+ **ipsm**类库，共**7**个参数

| 参数符号 | 解释 |
| :- | :- |
| -inlayers |  输入环境因子图层文件名（以#分割） |
| -datatypes |  每个环境因子图层的数据类型，包含类别型(categorical)和连续型(continuous)，（以#分割） |
| -sample |  样点数据文件名(csv格式)，文件中包含样点地理坐标(名为:x,y,大小写均可,并且与环境因子图层的空间参考对应)和要推测的属性值 |
| -target |  要推测的属性名(存在于样点文件的表头中) |
| -simithred |  环境相似度阈值，默认值0.5（可选） |
| -predmap |  输出的推测图层文件名 |
| -uncmap |  输出的不确定性图层文件名 |

    参数设置示例：
    ```-inlayers D:/data/xc/geo.asc#D:/data/xc/planc.asc#D:/data/xc/preci.asc#D:/data/xc/profc.asc#D:/data/xc/slope.asc#D:/data/xc/tempr.asc#D:/data/xc/twi.asc -datatypes categorical#continuous#continuous#continuous#continuous#continuous#continuous -sample D:/data/xc/samples_xc.csv -target SOMB -simithred 0.5 -predmap D:/data/test/pred.tif -uncmap D:/data/test/unc.tif```

#### 3.4.2 调用示例
    在生成编译程序（ipsm.exe等）的bin文件夹，打开命令提示符（cmd.exe），依次输入目标执行程序名和调用参数，执行程序。调用示例：
    ```ipsm -inlayers D:/data/xc/geo.asc#D:/data/xc/planc.asc#D:/data/xc/preci.asc#D:/data/xc/profc.asc#D:/data/xc/slope.asc#D:/data/xc/tempr.asc#D:/data/xc/twi.asc -datatypes categorical#continuous#continuous#continuous#continuous#continuous#continuous -sample D:/data/xc/samples_xc.csv -target SOMB -simithred 0.5 -predmap D:/data/test/pred.tif -uncmap D:/data/test/unc.tif```

## 4. 文献

### 4.1 推测理论
- <b> Geographic similarity </b>

    - Zhu, A-X., Lv, G. N., Liu, J., Qin, C.-Z., Zhou, C. H., 2018. Spatial prediction based on Third Law of Geography. Annals of GIS 24(4), 225-240. https://doi.org/10.1080/19475683.2018.1534890 

    - Zhu, A-X., Turner, M., 2022. How is the Third Law of Geography different?. Annals of GIS 28(1), 57-67. https://doi.org/10.1080/19475683.2022.2026467

### 4.2 iPSM 方法
- <b> iPSM </b>
    - Zhu, A-X., Liu, J., Du, F., Zhang, S. J., Qin, C.-Z., Burt, J. E., Behrens, T., Scholten, T., 2015. Predictive soil mapping with limited sample data. European Journal of Soil Science. 66(3), 535-547. https://doi.org/10.1111/ejss.12244

- <b> iPSM considering data availability (iPSM-FilterNA)</b>
    - Fan, N. Q., Zhu, A-X., Qin, C. Z., Liang, P. 2020. Digital soil mapping over large areas with invalid environmental covariate data. ISPRS International Journal of Geo-Information. 9(2), 102. http://dx.doi.org/10.3390/ijgi9020102 

- <b> iPSM considering adaptive covariate applicability (iPSM-WCovar)</b>
    - Fan, N.-Q., Zhao, F.-H., Zhu, L.-J., Qin, C.-Z., Zhu, A-X. 2022. Digital soil mapping with adaptive consideration of the applicability of environmental covariates over large areas. International Journal of Applied Earth Observation and Geoinformation 113, 102986. https://doi.org/10.1016/j.jag.2022.102986

- <b> iPSM considering spatial distances (iPSM-IDW) </b>
    - Fan, X., Fan, N., Qin, C.-Z., Zhao, F.-H., Zhu, L.-J., Zhu, A-X. 2023. Large-area soil mapping based on environmental similarity with adaptive consideration of spatial distance to samples. Geoderma 439, 116683. https://doi.org/10.1016/j.geoderma.2023.116683

    - Qin, C.-Z., An Y., Liang, P., Zhu, A-X., Yang, L. 2021. Soil property mapping by combining spatial distance information into the Soil Land Inference Model (SoLIM). Pedosphere, 31(4), 638-644. https://doi.org/10.1016/S1002-0160(20)60016-9 


### 4.3 软件
- <b><a href=https://lreis2415.github.io/SoLIMSolutions/software.html> Software: SoLIM Solutions</a></b>
    - Zhu, A-X., Qin, C.-Z., Liang, P., Du, F., 2018. Digital Soil Mapping for Smart Agriculture: the SoLIM method and Software Platforms. RUDN Journal of Agronomy and Animal Industries, 13(4), 317-335. http://dx.doi.org/10.22363/2312-797X-2018-13-4-317-335

- <b> <a href=https://lreis2415.github.io/SoLIMSolutions/software.html> Software: iSoLIM</a></b>
    - Zhao, F.-H., Zhu, A-X., Zhu, L.-J., & Qin, C.-Z. 2024. iSoLIM: a similarity-based spatial prediction software for the big data era. Annals of GIS, 1-15. https://doi.org/10.1080/19475683.2024.2324381 

- <b> <a href=http://www.easygeoc.net:8090/>Webservices: EasyGC</a></b>
    - Jiang J., Zhu, A-X., Qin, C.-Z., Zhu, T., Liu, J., Du, F., Liu, J., Zhang, G., An, Y. 2016. CyberSoLIM: A Cyber Platform for Digital Soil Mapping. Geoderma, 263, 234-243. https://doi.org/10.1016/j.geoderma.2015.04.018 


