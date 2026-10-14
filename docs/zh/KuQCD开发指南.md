# KuQCD开发指南<a name="ZH-CN_TOPIC_0000002424460372"></a>

- **[总述](#ZH-CN_TOPIC_0000002458154593)**

- **[KuQCD函数说明](#ZH-CN_TOPIC_0000002458074433)**

## 总述<a name="ZH-CN_TOPIC_0000002458154593" id="ZH-CN_TOPIC_0000002458154593"></a>

KuQCD是基于鲲鹏平台的格点量子色动力学领域加速库，采用C/C++实现，针对求解器等功能进行优化，并提供相应接口，以便在应用中使用。目前KuQCD仅支持在鲲鹏920新型号上使用，依赖KUPL，支持毕昇编译器和gcc编译器。

## KuQCD函数说明<a name="ZH-CN_TOPIC_0000002458074433" id="ZH-CN_TOPIC_0000002458074433"></a>

- **[数据结构](#ZH-CN_TOPIC_0000002458171977)**

- **[函数定义](#ZH-CN_TOPIC_0000002424493252)**

### 数据结构<a name="ZH-CN_TOPIC_0000002458171977" id="ZH-CN_TOPIC_0000002458171977"></a>

- **[Lattice](#ZH-CN_TOPIC_0000002461238893)**

- **[LattBase](#ZH-CN_TOPIC_0000002424507534)**

- **[CommBase](#ZH-CN_TOPIC_0000002458106157)**

- **[SpinorField](#ZH-CN_TOPIC_0000002424667406)**

- **[GaugeField](#ZH-CN_TOPIC_0000002458186257)**

#### Lattice<a name="ZH-CN_TOPIC_0000002461238893" id="ZH-CN_TOPIC_0000002461238893"></a>

lattice为格点，结构体里参数如[表1](#table-1)所示。

**表1** Lattice结构体<a name="table-1" id="table-1"></a>

|变量名|类型|描述|
|--|--|--|
|nx|unsigned short|格点x方向大小|
|ny|unsigned short|格点y方向大小|
|nz|unsigned short|格点z方向大小|
|nt|unsigned short|格点t方向大小|
|V|size_t|格点整体大小，值为nx*ny*nz*nt|

#### LattBase<a name="ZH-CN_TOPIC_0000002424507534" id="ZH-CN_TOPIC_0000002424507534"></a>

LattBase结构体参数如[表2](#table-2)所示。

**表2** LattBase结构体<a name="table-2" id="table-2"></a>

|变量名|类型|描述|
|--|--|--|
|latSupBulkFL|Lattice|全进程，无halo，全格点|
|latSupBulkCB|Lattice|全进程，无halo，奇或偶格点|
|latSubFullFL|Lattice|子进程，有halo，全格点|
|latSubFullCB|Lattice|子进程，有halo，奇或偶格点|
|latSubBulkFL|Lattice|子进程，无halo，全格点|
|latSubBulkCB|Lattice|子进程，无halo，奇或偶格点|
|latFull|Lattice *|根据site_subset，为子进程有halo格点|
|latBulk|Lattice *|根据site_subset，为子进程无halo格点|
|lenFull|size_t|根据site_subset，为子进程有halo总格点大小|
|lenBulk|size_t|根据site_subset，为子进程无halo总格点大小|
|gridDim|array<unsigned short, 4>|各维度进程切分策略|
|haloDepth|unsigned short|halo的深度，向量时设置为0，矩阵时为3|
|Nbatch|unsigned int|x方向批量大小|
|myRank|int|当前进程所属进程号|
|np|int|进程数|
|prec|KuQCDPrecision|精度选择，有如下值：<ul><li>KUQCD_HALF</li><li>KUQCD_SINGLE</li><li>KUQCD_DOUBLE</li></ul>|
|gauge_order|KuQCDGaugeOrder|是否做x方向批量，有如下值：<ul><li>KUQCD_NORMAL_ORDER</li><li>KUQCD_BATCH_ORDER</li></ul>|
|rankOrder|KuQCDRankOrder|格点排布顺序，有如下值：<ul><li>XYZT_ORDER</li><li>TZYX_ORDER</li></ul>|
|site_subset|KuQCDSiteSubset|子格点，有如下值：<ul><li>KUQCD_EVEN</li><li>KUQCD_ODD</li><li>KUQCD_ALL</li></ul>|
|myIdx|Index|格点所在的index|

#### CommBase<a name="ZH-CN_TOPIC_0000002458106157" id="ZH-CN_TOPIC_0000002458106157"></a>

CommBase结构体如[表3](#table-3)所示，仅规范场用到。

**表3** CommBase结构体<a name="table-3" id="table-3"></a>

|变量名|类型|描述|
|--|--|--|
|latBase|LattBase &|所用到的格点|
|rankHyperPlane|array<array<int, 2>, 4>|一维通信进程号|
|lenHyperPlane|array<size_t, 4>|一维通信包长度|
|rankPlane|array<array<int, 4>, 6>|二维通信进程号|
|lenPlane|array<size_t, 6>|二维通信包长度|
|rankStripe|array<array<int, 8>, 4>|三维通信进程号|
|lenStripe|array<size_t, 4>|三维通信包长度|
|rankCorner|array<int, 16>|四维通信进程号|
|lenCorner|size_t|四维通信包长度|
|nbIdx|Index|邻居格点序号|
|haloIdxSend1|array<array<Index, 2>, 4>|一维通信需要发送的位置序号|
|haloIdxRecv1|array<array<Index, 2>, 4>|一维通信需要接收的位置序号|
|haloIdxSend2|array<array<Index, 4>, 6>|二维通信需要发送的位置序号|
|haloIdxRecv2|array<array<Index, 4>, 6>|二维通信需要接收的位置序号|
|haloIdxSend3|array<array<Index, 8>, 4>|三维通信需要发送的位置序号|
|haloIdxRecv3|array<array<Index, 8>, 4>|三维通信需要接收的位置序号|
|haloIdxSend4|array<Index, 16>|四维通信需要发送的位置序号|
|haloIdxRecv4|array<Index, 16>|四维通信需要接收的位置序号|

#### SpinorField<a name="ZH-CN_TOPIC_0000002424667406" id="ZH-CN_TOPIC_0000002424667406"></a>

SpinorField为旋量场，属于模板类，模板的具体内容如[表4](#table-4)所示。

**表4** SpinorField模板<a name="table-4" id="table-4"></a>

|模板名|类型|描述|
|--|--|--|
|T|typename|旋量场的数据类型，取值可为float,double|
|Nbatch|int|x方向批量大小，对于float可为1或16，对于double为1|
|Ns|int|旋量场分量数，对于Wilson型费米子为4，对于Stagger型费米子为1|


SpinorField的参数的具体内容如[表5](#table-5)所示。

**表5** SpinorField参数<a name="table-5" id="table-5"></a>

|变量名|类型|描述|
|--|--|--|
|latBase|LattBase &|所用到的格点|
|data|ColorSpinor<T,Nbatch>|数据|
|locSize|int|每个时空格点储存浮点数个数|
|length|size_t|总时空格点数，考虑奇偶格点|

#### GaugeField<a name="ZH-CN_TOPIC_0000002458186257" id="ZH-CN_TOPIC_0000002458186257"></a>

GaugeField为规范场，属于模板类，模板的具体内容如[表6](#table-6)所示。

**表6** GaugeField模板<a name="table-6" id="table-6"></a>

|模板名|类型|描述|
|--|--|--|
|T|typename|规范场的数据类型，取值可为float,double|
|Nbatch|int|x方向批量大小，对于float可为1或16，对于double为1|
|Recon|KuQCDGaugeReconstruct|SU(3)矩阵压缩方案，目前仅可取值KUQCD_RECON_NO|


GaugeField的参数的具体内容如[表7](#table-7)所示。

**表7** GaugeField参数<a name="table-7" id="table-7"></a>

|变量名|类型|描述|
|--|--|--|
|even|SU3<T, Nbatch, Recon>|指向规范场的偶部分|
|odd|SU3<T, Nbatch, Recon>|指向规范场的奇部分|
|comm|CommBase &|规范场通信|
|commHyperPlane|array<array<Communicator<GaugeT>, 2>, 4>|一维通信子|
|commPlane|array<array<Communicator<GaugeT>, 4>, 6>|二维通信子|
|commStripe|array<array<Communicator<GaugeT>, 8>, 4>|三维通信子|
|commCorner|array<Communicator<GaugeT>, 16>|四维通信子|
|latBase|LattBase &|所用到的时空格点|
|data|SU3<T, Nbatch, Recon>*|数据|
|locSize|int|每个时空格点储存浮点数个数|
|length|size_t|总时空格点数，全格点|

### 函数定义<a name="ZH-CN_TOPIC_0000002424493252" id="ZH-CN_TOPIC_0000002424493252"></a>

- **[kuqcd_init](#ZH-CN_TOPIC_0000002424667410)**

- **[kuqcd_load_gaugefield](#ZH-CN_TOPIC_0000002424507542)**

- **[kuqcd_solve_multishift](#ZH-CN_TOPIC_0000002458106165)**

- **[kuqcd_gauge_free](#ZH-CN_TOPIC_0000002424667414)**

- **[kuqcd_end](#ZH-CN_TOPIC_0000002458186269)**

- **[示例](#ZH-CN_TOPIC_0000002461420081)**

#### kuqcd_init<a name="ZH-CN_TOPIC_0000002424667410" id="ZH-CN_TOPIC_0000002424667410"></a>

初始化

>![](public_sys-resources/icon-notice.gif) **须知：** 
>程序需要在使用其他函数前调用该函数一次，且只允许调用一次，运行时需添加输入文件，若可执行文件为run.exe，输入文件为input，需用以下命令运行：mpirun <xxx mpi参数> ./run.exe < input

**接口定义<a name="section18456193311376"></a>**

```cpp
KuQCDStatus kuqcd_init();
```

**返回值<a name="section65022020475"></a>**

|返回值|类型|描述|
|--|--|--|
|KUQCD_SUCCESS|int|执行成功|
|KUQCD_ERR_INVALID_PARAM|int|无效参数|
|KUQCD_ERR_UNSUPPORTED|int|不支持的选项|
|KUQCD_ERR_NULL_POINTER|int|空指针|

#### kuqcd_load_gaugefield<a name="ZH-CN_TOPIC_0000002424507542" id="ZH-CN_TOPIC_0000002424507542"></a>

加载规范场对象

**接口定义<a name="section18456193311376"></a>**

```cpp
KuQCDStatus kuqcd_load_gaugefield(void *h_gauge, KuQCDGaugeParam *param);
```

**参数<a name="section18911913153917"></a>**

|参数名|类型|描述|输入/输出|
|--|--|--|--|
|h_gauge|void *|规范场的数据|输入|
|param|KuQCDGaugeParam *|规范场的配置选项，具体内容如[表8](#table-8)所示|输入|

**表8** KuQCDGaugeParam结构体字段描述<a name="table-8" id="table-8"></a>

|字段名|类型|描述|
|--|--|--|
|lat|unsigned int [4]|x,y,z,t方向格点数|
|halo_depth|int|halo区大小，对于旋量场取0；对于规范场，若做rhmc（Rational Hybrid Monte Carlo），取3，只做inverter，取0。|
|od|KuQCDRankOrder|进程排布顺序，可取KUQCD_XYZT_ORDER（x方向快变）或KUQCD_TZYX_ORDER（t方向快变）|
|site_subset|KuQCDSiteSubset|子格点，可取值KUQCD_EVEN, KUQCD_ODD, KUQCD_ALL|
|gauge_order|KuQCDGaugeOrder|数据排布顺序，可取：<ul><li>KUQCD_NORMAL_ORDER（对于旋量场(even-odd)-y-z-y-x-Ns-3-(real-imag)，对于规范场(even-odd)-t-z-y-x-direction-col-row-(real-imag)）</li><li>KUQCD_BATCH_ORDER（对于旋量场(even-odd)-y-z-y-(x/Nbatch)-Ns-3-(real*Nbatch-imag*Nbatch)，对于规范场(even-odd)-t-z-y-(x/Nbatch)-direction-col-row-(real*Nbatch-imag*Nbatch)）</li></ul>|
|prec|KuQCDPrecision|精度，可取KUQCD_SINGLE和KUQCD_DOUBLE|
|link_type|KuQCDLinkType|链接类型，取值有：<ul><li>KUQCD_WILSON_LINK（暂不支持）</li><li>KUQCD_FAT_LINK</li><li>KUQCD_LONG_LINK</li></ul>|
|recon|KuQCDReconstructType|SU(3)矩阵压缩方式，取值有：<ul><li>KUQCD_GAUGE_RECON_NO</li><li>KUQCD_GAUGE_RECON_U14（暂不支持）</li><li>KUQCD_GAUGE_RECON_12（暂不支持）</li><li>KUQCD_GAUGE_RECON_8（暂不支持）</li></ul>|

**返回值<a name="section65022020475"></a>**

|返回值|类型|描述|
|--|--|--|
|KUQCD_SUCCESS|int|执行成功|
|KUQCD_ERR_INVALID_PARAM|int|无效参数|
|KUQCD_ERR_UNSUPPORTED|int|不支持的选项|
|KUQCD_ERR_NULL_POINTER|int|空指针|

#### kuqcd_solve_multishift<a name="ZH-CN_TOPIC_0000002458106165" id="ZH-CN_TOPIC_0000002458106165"></a>

进行求解

**接口定义<a name="section18456193311376"></a>**

```cpp
KuQCDStatus kuqcd_solve_multishift(void **_hp_x, void *_hp_b, KuQCDInvertParam *param);
```

**参数<a name="section18911913153917"></a>**

|参数名|类型|描述|输入/输出|
|--|--|--|--|
|_hp_x|void **|指向多偏移的解数据|输入/输出|
|_hp_b|void *|指向源数据|输入|
|param|KuQCDInvertParam *|求解进行配置选项，具体内容如[表9](#table-9)所示|输入|

**表9** KuQCDInvertParam结构体字段描述<a name="table-9" id="table-9"></a>

|字段名|类型|描述|
|--|--|--|
|dslash_type|KuQCDDslashType|狄拉克算子类型，取值有：<ul><li>KUQCD_WILSON_DIRAC（暂不支持）</li><li>KUQCD_WILSONPC_DIRAC（暂不支持）</li><li>KUQCD_HISQ_DIRAC</li><li>KUQCD_HISQPC_DIRAC</li></ul>|
|inverter_type|KuQCDInverterType|求解器类型，取值有：<ul><li>KUQCD_CG（暂不支持）</li><li>KUQCD_CGM</li><li>KUQCD_BICGSTAB（暂不支持）</li><li>KUQCD_BICGSTABM（暂不支持）</li></ul>|
|matpc_type|KuQCDMatPCType|预条件子类型，取值有：<ul><li>KUQCD_EVEN_EVEN</li><li>KUQCD_ODD_ODD</li></ul>|
|prec|KuQCDPrecision|求解精度，需要和KuQCDGaugeParam参数的值保持一致，取值有：<ul><li>KUQCD_HALF（暂不支持）</li><li>KUQCD_SINGLE</li><li>KUQCD_DOUBLE</li></ul>|
|iter|int|迭代次数，求解完之后读取该值可获取求解迭代的次数|
|num_shift|int|解的个数|
|max_iter|int|最大迭代次数|
|tol|double|收敛精度|
|shifts|double *|偏移量数组|

#### kuqcd_gauge_free<a name="ZH-CN_TOPIC_0000002424667414" id="ZH-CN_TOPIC_0000002424667414"></a>

释放规范场对象

**接口定义<a name="section18456193311376"></a>**

```cpp
void kuqcd_gauge_free();
```

#### kuqcd_end<a name="ZH-CN_TOPIC_0000002458186269" id="ZH-CN_TOPIC_0000002458186269"></a>

求解结束

**接口定义<a name="section18456193311376"></a>**

```cpp
void kuqcd_end();
```

#### 示例<a name="ZH-CN_TOPIC_0000002461420081" id="ZH-CN_TOPIC_0000002461420081"></a>

```cpp
#include <omp.h>
#include <mpi.h>
#include <cstdio>
#include "kuqcd.h"

using namespace kuqcd;

template<typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void TestCGMPC(const KuQCDGaugeOrder gauge_order);


int main (int argc, char** argv)
{
    int provided = 0;
    int required = MPI_THREAD_FUNNELED;
    int flag = MPI_Init_thread(&argc, &argv, required, &provided);

    TestCGMPC<float, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>(KUQCD_NORMAL_ORDER);
    MPI_Finalize();
    return 0;
}


template<typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void TestCGMPC(const KuQCDGaugeOrder gauge_order)
{
    int repos = 1;
    int nx = 4, ny = 4, nz = 4, nt = 1, seed = 1;
    const int num_shift = 1;
    std::array<unsigned short, 4> gridDim{1, 1, 1, 1};
    KuQCDRankOrder od = KUQCD_TZYX_ORDER;
    KuQCDSiteSubset site_subset = KUQCD_ODD;
    int haloDepthSpin = 0;
    int haloDepthGauge = 3;

    Lattice latSup(nx, ny, nz, nt);

    KuQCDPrecision prec = (sizeof(T)==sizeof(double)) ? KUQCD_DOUBLE : KUQCD_SINGLE;

    rootLogger.Info("initializing spinor field lattice base");
    LattBase *latBaseOdd = nullptr;
    KuQCDGaugeParam Gagugeoptions;
    Gagugeoptions.lat[0] = nx;
    Gagugeoptions.lat[1] = ny;
    Gagugeoptions.lat[2] = nz;
    Gagugeoptions.lat[3] = nt;
    Gagugeoptions.gauge_order = gauge_order;
    Gagugeoptions.od = od;
    Gagugeoptions.prec = prec;
    site_subset = KUQCD_ALL;
    Gagugeoptions.site_subset = site_subset;
    Gagugeoptions.halo_depth = haloDepthGauge;


    double shifts[num_shift];
    shifts[0] = 3.67167631338816e-05;
    FILE* fp;
    float* origin;
    int np = GetNp();
    int rk = GetRank();

    size_t vol = (size_t)nx *ny * nt *nz / Nbatch / (np);
    size_t volh = vol / 2;
    size_t spinorLen = volh* Nbatch * 6;
    size_t gaugeLen = vol* Nbatch * 4 * ReconS;
    size_t gaugeLenNaik = vol* Nbatch * 4 * ReconN;

    printf("vol is %zu, spinorLen %zu, gaugeLen %zu\n",vol,spinorLen,gaugeLen);
    origin = (float*)malloc(((num_shift + 1)*spinorLen+gaugeLen+gaugeLenNaik)*sizeof(float));
    char name[100];

    for(int i = 0; i<(num_shift + 1)*spinorLen+gaugeLen+gaugeLenNaik;i++) {
        origin[i] = 1.0;
    }
    kuqcd_init();
    float *originSrc = origin;
    float *origingaugeSmeared = &(origin[(num_shift+1)*spinorLen]);
    float *origingaugeNaik = &(origin[(num_shift+1)*spinorLen+gaugeLen]);

    rootLogger.Info("SpinorFieldSetValues1 done");
    Gagugeoptions.link_type = KUQCD_FAT_LINK;
    Gagugeoptions.recon = ReconS;
    kuqcd_load_gaugefield((void*)origingaugeSmeared, &Gagugeoptions);
    Gagugeoptions.link_type = KUQCD_LONG_LINK;
    Gagugeoptions.recon = ReconN;
    kuqcd_load_gaugefield((void*)origingaugeNaik, &Gagugeoptions);

    rootLogger.Info("initializing SetValues done");

    KuQCDInvertParam invertOpts;
    invertOpts.dslash_type = KUQCD_HISQPC_DIRAC;
    invertOpts.inverter_type = KUQCD_CGM;
    invertOpts.matpc_type = KUQCD_ODD_ODD;
    int iter = 0;
    invertOpts.num_shift = num_shift;
    invertOpts.max_iter = 6500;
    invertOpts.tol = 1e-12;
    invertOpts.shifts = shifts;
    invertOpts.prec =prec;
    void *dst[num_shift] = {nullptr};
    float *dst_value[num_shift];
    for(int i = 0;i<num_shift;i++) {
        dst_value[i] = (float*)malloc(spinorLen * sizeof(float));
        dst[i] = (void *)(dst_value[i]);
    }

    double tStart = MPI_Wtime();

    kuqcd_solve_multishift(dst, (void *)originSrc, &invertOpts);
    double tEnd = MPI_Wtime();
    Rank0Printf("done\n");
    fflush(0);
    Rank0Printf("total time =  %e [s]\n", tEnd - tStart);
    Rank0Printf("time/iter  =  %e [s]\n", (tEnd - tStart) / invertOpts.iter);
    fflush(stdout);



}
```

