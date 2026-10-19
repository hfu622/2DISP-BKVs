# 复现源码

本目录包含 Rust/Jagua 算法实现、几何原语和优化后端。获选方法为 GCS 的 Dighe 零聚类分支、SNCPH、3PM、CHS2 和 SemiDiscreteBLF；其他模块为研究实现，不构成额外的 BKV 归属。

论文、逐样例对照值及限制见[项目说明](../README.md)。[配置文件](../config/reproduction_decisions.json)记录未报告参数和实际选择，并被编译嵌入；[测试](tests)用于检查几何和算法原语。限时子问题返回可行解与已证明最优分开判断；几何合法不证明算法忠实或全局紧凑。
