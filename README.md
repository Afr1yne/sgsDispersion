# 亚格子扩散模型

OpenCFD OpenFOAM **v2606** 用户库。运行时名称 `sgsDispersion`，类名 `SGSDispersion`，库名 `libsgsDispersion.so`。接入经典 `basicKinematicCloud`，本机原生测试使用 `kinematicParcelFoam`。

模型保存每个粒子所见流体的 SGS 速度 `UTurb`，以 `Uc + UTurb` 参与原有颗粒受力计算。生产率、解析速度方差和滑动时间窗口在模型内计算。无需配置 `controlDict/functions`。

论文式（4）—（8）的 SGS 随机模型已逐项核对，并通过独立解析基准及原生云测试；同库提供球形颗粒的论文式（11）阻力。具体测试数据、修正记录和实现约定见 [论文一致性与验证](论文一致性与验证.md)。该项目实现 SGS 粒子扩散库；载气求解与粒子运动积分采用下述 OpenFOAM 设置。

## 目录和编译

把仓库放在 `$WM_PROJECT_USER_DIR/src/sgsDispersion`，算例生成到 `$FOAM_RUN/sgsDispersion`。源码集中在一个开发目录。

```text
src/SGSDispersion/   模型、统计、随机更新与阻力
src/Make/           原生用户库编译配置
tests/              独立基准及原生缓存检查
tools/              算例生成、统计和点云绘制
validation/         本次实际核验结果和小型记录
```

```bash
# 先加载自己的 OpenFOAM v2606 环境
source /path/to/OpenFOAM-v2606/etc/bashrc
cd "$WM_PROJECT_USER_DIR/src/sgsDispersion"
bash Allwmake
bash Allvalidate
```

`Allwmake` 编译库、运行内部统计及随机过程测试，并编译原生缓存测试。`Allvalidate` 生成三个时间步的规定剪切流小惯性粒子测试，检查原生 SGS 场和 Taylor 扩散方差，再运行变时间步缓存检查。该组使用标准 kEqn，入口为固定值，可独立于本机数字滤波入口修改运行。首次验证建议使用新的 `validation/uniform` 算例目录；已有目录会复用其中的输入字典。

需要 C++17 编译器、OpenFOAM v2606 及 Python 3；基础检查与算例统计只使用 Python 标准库。点云图片需要 ParaView 的 pvpython，拼图另需 Pillow。未核验其他 OpenFOAM 分支或版本。

## 在已有 LES 算例中选择模型

`system/controlDict`：

```foam
libs ("libsgsDispersion.so");
```

`constant/kinematicCloudProperties` 的 `subModels`：

```foam
dispersionModel sgsDispersion;
sgsDispersionCoeffs
{
    averagingTime 400;
    C0 4;
    beta 2;
    cEpsilon 0.93;
    verticalDirection (0 0 1);
    stationaryInitial false;
}
```

`averagingTime` 是必填的内部时间平均窗口，单位 s。本例按域高和摩擦速度选择 H/u*=120/0.3=400 s。首步使用瞬时生产率、零解析方差；随后使用已积累窗口，完整 ABL 对比在窗口填满后释放颗粒。`stationaryInitial false` 从零 SGS 速度积累噪声；`true` 对初始 SGS 速度采样局部高斯分布，用于平衡态基准。

扩散库可与现有受力模型组合。若要使用本项目按论文式（11）实现的球形颗粒阻力，将原阻力条目替换为同库注册的 `schillerNaumannDrag`，例如：

```foam
particleForces
{
    schillerNaumannDrag;
    gravity;
}
```

三组 ABL 对比统一采用上述阻力和重力，单向耦合。阻力修正严格使用 `1 + 0.15 Re^0.687`，未另设高 Reynolds 数分支；新工况应检查该经验关系的适用范围。

## 公式和离散约定

依据 Comola et al. (2019), *Preferential Deposition of Snow and Dust Over Hills: Governing Processes and Relevant Scales*, [DOI](https://doi.org/10.1029/2018JD029614)、[作者 PDF](https://weather.ou.edu/~salesky/papers/comola19jgr.pdf)，第 7955—7956 页。

```text
du_s,i = [-alpha*u_s,i/Tp
          + 0.5*((D sigma²/Dt)*u_s,i/sigma² + ∂sigma²/∂x_i)]dt
          + sqrt(2*alpha*sigma²/Tp)dW_i

epsilon = <2*nut*S:S>, S = symm(grad(U))
e = (epsilon*delta/cEpsilon)^(2/3)
kResolved = 0.5*<|U-<U>|²>
alpha = e/(e+kResolved), sigma² = 2*e/3
Tf = 2*sigma²/(C0*epsilon)
Tp = Tf/sqrt(1+(beta*Up_vertical)²/sigma²)
```

`Up_vertical` 为粒子绝对竖向速度。`delta` 读取 LES 注册的实际滤波宽度场。生产率表达式用于不可压缩涡黏性 SGS 闭合；`e` 按平均生产率重建，载气输运的 LES `k` 用于其自身闭合。

时间统计按每个已完成载气步的时长加权，最旧区间可部分移除，速度方差采用中心矩更新。`kinematicParcelFoam` 先演化云、再解载气，因此可用场的时间为当前时间减当前步长；变化的场时间间隔用于 SGS 方差差分，已用原生不等步长测试核对。改接其他求解器时需要核查这一顺序。

重颗粒的轨迹导数采用 `∂t sigma² + Up·grad(sigma²)`；空间系数取所在单元，梯度使用 `fvc::grad`。这是本库的有限体积实现约定。冻结局部线性系数的 SDE 采用精确转移，内部噪声子步保持当地系数。粒子位置和速度由原生跟踪器及 `analytical` 速度积分更新。论文的载气为静态 Smagorinsky、谱方法/有限差分和 IBM，粒子动量采用 Verlet；本库的示例采用 kEqn 和经典云。论文未提供这些移植细节或可供逐行核对的作者源码。

## ABL 三组对比

```bash
mkdir -p "$FOAM_RUN/sgsDispersion"
cp tools/Allrun "$FOAM_RUN/sgsDispersion/Allrun"
bash "$FOAM_RUN/sgsDispersion/Allrun" smoke
bash "$FOAM_RUN/sgsDispersion/Allrun" abl
```

```text
abl/RANS_DRW       固定 U/k/epsilon 剖面 + stochasticDispersionRAS
abl/LES_resolved   LES + dispersionModel none
abl/LES_SGS        LES + sgsDispersion
smoke/            短时接入检查
results/abl/      三组统计、原生点云、图片和来源记录
results/validation/ 原生解析基准统计
```

生成器拒绝覆盖已有字典；运行入口保留已有正常结束并含终点粒子文件的算例。修改输入后请在新目录重算。

本机两组 LES 使用已有的 seed/逻辑索引确定性数字滤波入口，seed=20261011。入口修改不属于本库；本次载气一致性在该环境下验证。复现两组载气场逐文件相同的结果需要匹配本机入口实现；入口源码校验值在本机 `results/abl/provenance.json`。`Allvalidate` 可作为不依赖该修改的安装核验入口。

## 结果和适用范围

模型输出以云名开头的 `SGSVariance`、`SGSFraction`、`SGSDissipation`、`SGSResolvedK`、`SGSWindowDuration` 诊断场。原生点云在各组 `VTK/lagrangian/kinematicCloud`，统计为 `summary.json` 和 `dispersion.csv`。

已核验固定网格、串行、从零开始、单向耦合的经典云。非均匀方差下的流体示踪极限（beta=0）well-mixed 检查及删除梯度项的反例检查已通过；该结果限于独立随机/位置更新基准，尚未完成全原生非均匀 FV 网格、重颗粒沉积实验或 DNS 验证。ABL 三组点云表明给定设置下 SGS 产生额外扩散。点云宽度不能单独证明物理精度。

历史只在内存中保存，重启会重新积累统计。连续重启、MPI、多随机种子、ABL 网格收敛尚未核验；不支持拓扑改变。窗口样本数值存储约为 `32*nCells*nSamples` 字节，另有容器及临时场开销。3456 单元、800 样本约 88.5 MB。

独立 ParaView 绘图时可使用：

```bash
env -u LD_LIBRARY_PATH VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow \
 /path/to/ParaView/bin/pvpython --force-offscreen-rendering \
 "$WM_PROJECT_USER_DIR/src/sgsDispersion/tools/renderClouds.py"
```

图片每组显示 6000 个原生点，蓝色 20 μm、红色 80 μm，固定相机和显示尺度；标记大小用于显示。图片来源、校验值和相机保存在 `figure_provenance.json`。

许可证：GPL-3.0-or-later，见 [LICENSE.md](LICENSE.md)。
