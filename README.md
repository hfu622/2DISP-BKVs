# 2DISP-BKVs

本项目整理 Nesting 25 个样例在已审核文献中的非 learning、非精确 2DISP 候选 BKV 归属，并比较论文报告与作者代码或本项目实现得到的值。BKV 表示已知最好值，不表示数学最优；输入或配置尚未确认的来源单独标明。

## 1. 数据集来源

数据来自 [ESICUP 官方 2d_irregular 数据集](https://github.com/ESICUP/datasets/tree/main/2d_irregular)。本仓库 [25 个官方 XML](official_esicup/xml) 保留轮廓、数量、初始朝向与合法角度；`shapes2.xml` 对应上游 `blaz/blaz.xml`，文件内容不变。

## 2. 逐样例算法归属与引用

| 当前算法 | 对应样例 | 代码来源 |
|---|---|---|
| Sparrow [1] | Albano、Dagli、Fu、Jakobs1/2、Mao、Marques、Shapes0/1/2、Shirts、Swim、Trousers | 作者 MIT Rust 源码，固定提交 |
| GCS [2] | Dighe1/2 | 本项目 Rust/Jagua 专用实现 |
| F057 附表 SVGNest 固定序矢量对照 [3] | Han、Poly3b/4b | 作者 JS 包装与固定版本 SVGnest；提供按哈希获取源码的工具 |
| SNCPH [4] | Poly2a/3a/4a/5a | 本项目 Rust/Jagua 专用实现 |
| 3PM [5] | Poly1a | 本项目 Rust/Jagua、精确原生 DBM 与 Z3 有理数连续后端 |
| CHS2 [6] | Poly2b | 本项目 Rust/Jagua 与 CBC |
| SemiDiscreteBLF [7] | Poly5b | 本项目 Rust/Jagua 专用实现 |

1. Gardeyn、Vanden Berghe、Wauters (2025). [An open-source heuristic to reboot 2D nesting research](https://arxiv.org/abs/2509.13329v4)。[作者 Sparrow 代码](https://github.com/JeroenGar/sparrow)。
2. Elkeran (2013). [A new approach for sheet nesting problem using guided cuckoo search and pairwise clustering](https://doi.org/10.1016/j.ejor.2013.06.020)。
3. (2026). [Raster-Convolution Search for Two-Dimensional Irregular Packing with Curved and Holed Parts](https://doi.org/10.1109/AANN71145.2026.11709270)。本表采用其[附表固定序矢量对照](https://github.com/Hermione-HKX/Raster-Convolution-Search-and-CTM-Bench/blob/c7aaa7e8c31db15f4f7df7673b2cecfc362b21cd/results/svgnest_public_xml.csv)，不是论文提出的栅格方法。[作者包装源码](https://github.com/Hermione-HKX/Raster-Convolution-Search-and-CTM-Bench)与 [SVGnest](https://github.com/Jack000/SVGnest)。
4. Baldacci 等 (2014). [Algorithms for nesting with defects](https://doi.org/10.1016/j.dam.2012.03.026)。采用其中条带包装启发式分支。
5. (2016). [A model-based heuristic for the irregular strip packing problem](https://doi.org/10.1590/0101-7438.2016.036.03.0447)。
6. Antonio Martinez Sykora (2013). [Nesting problems: exact and heuristic algorithms](https://www.uv.es/marsyan/docs/thesis.pdf)，学位论文第 6 章构造启发式 CHS2。
7. Chehrazad、Roose、Wauters (2022). [A fast and scalable bottom-left-fill algorithm to solve nesting problems using a semi-discrete representation](https://doi.org/10.1016/j.ejor.2021.10.043)。

## 3. 复现方法

闭源获选方法使用 Rust 实现，几何基础为 [Jagua](https://github.com/JeroenGar/jagua-rs) 1.1.0。按论文定义实现 NFP/IFR、接触候选、价格更新、插入与局部压紧；保留凹槽和退化点、线接触。紧凑放置依据各论文的目标、邻域或离散包络判断：SNCPH 检查 phi 局部改善，半离散 BLF 搜索最早列和最低合法位置，3PM/CHS2 保持模型内最小长度目标，并区分限时可行与子问题已证明最优。几何合法性本身不能证明算法忠实或全局紧凑。

3PM 使用原生 DBM 离散搜索和 [Z3](https://github.com/Z3Prover/z3) 有理数连续优化；不能识别的模型回退 [CBC](https://github.com/coin-or/Cbc)。CHS2 使用 CBC 替代 CPLEX。启发式整体中包含精确优化子问题，与将完整 2DISP 作为精确算法求解有区别。布局按官方原轮廓独立用 Shapely/GEOS 复核；固定 0° 的退化情况另有有理数几何校验。不得通过缩小零件或显著改变条带宽度取得利用率。

[复现源码](reproductions/src)、[功能测试](reproductions/tests)、[作者源码及许可证](vendor)、[实际配置与未报告参数的处理](config/reproduction_decisions.json) 分别保留。源码许可证按各来源适用；Sparrow 为 MIT，Jagua 为 MPL-2.0。F057 包装许可未确认，故不转载包装本体，仅保留固定提交及 SHA256 获取清单。

## 4. 论文值与实测值

实测利用率按官方名义宽度计算，越大越好；长度已换算到论文对照单位，越小越好，不直接混用 XML 原始缩放单位。表中为当前披露配置下的实测最佳值。论文只报告长度时保留原长度，不臆造论文利用率。差异为“实测减论文”；利用率用百分点（pp），长度用 L。附表 F057 同时列出其舍入后的利用率和长度，极小差异主要受显示精度影响。

† 作者原输入或部分配置无法获取，论文值仅为条件参考；‡ 已接受的微小数值变体，详见限制。

| 样例 | 算法 | 论文报告 | 实测利用率 % | 实测长度 L | 差异 |
|---|---|---|---:|---:|---|
| albano | Sparrow | 89.82% | 89.785369 | 9695.863000 | -0.034631 pp |
| dagli | Sparrow | 90.17% | 89.825012 | 56.303917 | -0.344988 pp |
| dighe1† | GCS | 100.00% | 96.863871 | 103.237667 | -3.136129 pp |
| dighe2† | GCS | 100.00% | 99.737820 | 100.262870 | -0.262180 pp |
| fu‡ | Sparrow | 92.41% | 92.400787 | 30.843893 | -0.009213 pp |
| han | SVGNest-F057-fixed-order | 63.60%；L=53.131 | 63.603421 | 53.130952 | +0.003421 pp；-0.000048 L |
| jakobs1‡ | Sparrow | 89.26% | 89.090488 | 11.000052 | -0.169512 pp |
| jakobs2‡ | Sparrow | 87.73% | 86.071262 | 22.423280 | -1.658738 pp |
| mao | Sparrow | 86.87% | 86.829851 | 1697.535400 | -0.040149 pp |
| marques | Sparrow | 92.02% | 91.212726 | 75.837090 | -0.807274 pp |
| poly1a† | 3PM | L=15.8 | 66.386209 | 15.439954 | -0.360046 L |
| poly2a† | SNCPH | L=27.75 | 67.368523 | 30.396251 | +2.646251 L |
| poly2b† | CHS2 | L=32.6 | 63.342667 | 35.698686 | +3.098686 L |
| poly3a† | SNCPH | L=41.65 | 72.701201 | 42.250000 | +0.600000 L |
| poly3b | SVGNest-F057-fixed-order | 64.64%；L=47.182 | 64.643931 | 47.181536 | +0.003931 pp；-0.000464 L |
| poly4a† | SNCPH | L=55.5 | 71.026150 | 57.661876 | +2.161876 L |
| poly4b | SVGNest-F057-fixed-order | 66.49%；L=58.166 | 66.491027 | 58.165743 | +0.001027 pp；-0.000257 L |
| poly5a† | SNCPH | L=69.25 | 70.063058 | 73.068124 | +3.818124 L |
| poly5b† | SemiDiscreteBLF | L=69.5 | 64.647887 | 71.000000 | +1.500000 L |
| shapes0‡ | Sparrow | 69.98% | 69.366562 | 57.520510 | -0.613438 pp |
| shapes1‡ | Sparrow | 76.73% | 76.706854 | 52.016212 | -0.023146 pp |
| shapes2‡ | Sparrow | 86.23% | 86.009366 | 25.113544 | -0.220634 pp |
| shirts | Sparrow | 90.92% | 89.991151 | 60.005900 | -0.928849 pp |
| swim‡ | Sparrow | 79.83% | 78.915916 | 5604.745600 | -0.914084 pp |
| trousers | Sparrow | 92.62% | 92.333556 | 235.888020 | -0.286444 pp |

[机器可读的逐样例汇总](data/results_summary.json) 另保留实测均值、最差值、样本标准差（适用时）及来源链接；这些是科学结果摘要，不包含运行过程记录。

## 5. 限制、缺失信息与采用配置

获选归属是已审核来源中的候选 BKV，不能据此证明完整文献的全局 BKV，也不能把项目实测改善直接视为同实例论文纪录更新。16 例符合本项目输入比较规则，其中 9 例严格匹配，7 例为已量化的微小变体；Dighe1/2、Poly1a、Poly2a/b、Poly3a、Poly4a、Poly5a/b 共 9 例的作者原输入或部分配置无法获取，保留条件归属并使用官方 XML。作者原输入、预旋转或设置无法确认时，不宣称严格原作者等价；Poly 的合法角度以官方 XML 为准。

官方 Shirts 为 99 件，原本地副本为 85 件；此处采用官方版本。官方 Jakobs2 包含完整端点；Jakobs2/Swim 的微小轮廓变体保留可比性限制。Han 官方 XML 总面积为 1960，另一 XLS/PDF 版本为 1970，两者不混用。Sparrow 的 Fu、Jakobs1/2、Shapes0/1/2 六例采用了 0.01% 的微小条带扩宽；本表仍以官方名义宽度报告利用率，因此不应把它们当作严格原宽度零误差纪录。

| 算法 | 未说明、矛盾或不可获取的信息 | 本项目采用的配置 |
|---|---|---|
| SNCPH | K、MaxIt、模式数及逐样例栅格/面板设置未明确；价格增减仅规定为正 | K=20、MaxIt=2、每配置 10 模式，price=penalty=0.1、初始价格 1；有限面板长度 4H，栅格 512×128、1024×256、2048×512；SNH1/2 及二者初始化的 SNH3 各配 PH1/2/3，共 12 配置。此为项目约定，不能等同于恢复了作者未公开协议 |
| GCS | Lévy 指数、越界、替换取整未说明；正文与伪代码扩张基准矛盾 | Dighe 零聚类专用分支；15 巢、30 代、弃置率 0.25、停滞 200；Mantegna β=1.5、clip、每代替换最差 4 巢；按正文从当前长度扩张 1%、收缩 4%，独立保留最优合法布局；每次预算 600 秒含 NFP 准备 |
| 3PM | 作者参考顶点、Figure 6 的 s 未定义；批量规则文字与 Algorithm 4 冲突；CPLEX 不可用 | XML 第一个实际顶点；s=gmin，初邻域边长 4gmin；采用 Algorithm 4；g=2→1→0.5、md=3000；原顶点整数支持面，DBM/Z3 精确后端及 CBC 回退。子问题最优不代表全局装箱最优 |
| CHS2 | CPLEX 环境、CBC 推荐参数及热启动细节不可获取 | n_ini=12、逐 MIP 100 秒；串行 CBC2.10.13，整数预处理及 LP 预求解关闭，随机种子固定为 1；同模型可行 MIP start，旧切片继承但坐标自由，完整模型仍校验。求解器替代会影响限时轨迹 |
| SemiDiscreteBLF | 并列顺序及部分几何实现细节未规定 | R=0.5；单侧零长度顶点区间扩展并保留 L/R 标签，内部桥接关闭；等包围盒面积按输入索引稳定排序；每次上移重新检查，舍入停滞取下一个可表示 y，原轮廓否定候选后继续同列接触搜索 |
| Sparrow | 未固定 nightly 日期及历史传递依赖 | 作者固定源码，nightly-2026-01-15、SIMD、release、target-cpu=native、only_final_svg、锁定依赖；3 workers，960 秒探索及 240 秒压缩。硬件、系统与依赖差异仍会影响定时搜索量 |
| F057 附表矢量对照 | 包装许可未确认；不能把附表对照当成论文提出的方法 | 固定作者包装和 SVGnest 提交、固定顺序、官方角度、单次 placement worker；不添加遗传搜索或事后通用压紧 |

几何复核保留原轮廓和坐标：线性容差为 `1e-6 × max(条带宽度, 1)`，面积容差为 `1e-8 × max(零件总面积, 1)`；原 MIP 使用按行缩放的 `1e-7` 残差检查。允许的微小误差与严格零越界、零重叠不同，不能据容差内通过宣称严格精确可行。显著宽度改变、零件缩放或超容差几何不合法不能作为有效 BKV。

CHS2 当前配置实测利用率范围为 25.158157%–63.342667%，均值 55.857504%，样本标准差 9.084133 个百分点；与条件论文长度仍有差距，不能用“可行”推导为“复现性能等价”。所有无法获取的信息及实际选择集中记录于[配置文件](config/reproduction_decisions.json)。
