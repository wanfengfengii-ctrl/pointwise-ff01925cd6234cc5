# 校准包交接与求解正确性说明

本说明配合 CALIBRATION.md，描述两项联动能力的实现位置、固定修订包的交接
方式，以及有界等式校准二次最优权重是如何求得并被独立核验的。

## 模块与入口

- `survey/qp.py`：仅依赖标准库的有界等式凸二次规划求解器。
  - `bounded_minimum(D, q, A, b, lo, hi, denominators)` 返回唯一最优权重，
    目标为 `0.5 z'Dz - q'z`，约束 `Az=b`、`lo<=z<=hi`。
  - 校准目标 `sum (w-d)^2/d` 取 `D=1/d`、`q=1`，其梯度为 `Dw-1`。
- `survey/calibration.py`
  - `Calibrator.fit(snapshot, controls, targets, lower, upper)`：入口参数
    校验后，对总体和每个复制分别调用求解器；任一不可行即整体抛
    `CalibrationInfeasible`（`ValueError` 子类），不返回部分结果。
  - `CalibrationBundle`：冻结数据包，固定源快照的行身份、Row、修订、
    总体校准权重、全部复制校准权重与原系数。
- `survey/report.py`
  - `Reporter.ratios(bundle, domains)`：消费数据包，按域输出总体比例、
    每个复制比例、规定方差和缺失复制编号；`raw_ratios` 保持原未校准语义。
- `survey/__main__.py`：保留 `demo`，新增 `calibration-demo`，`all` 依次
  调用两者。运行清单 `.pointwise-runtime.json` 的 `demo` 指向 `all`，
  `verify` 仍用 unittest discover，自动同时发现原测试与新增测试。

## 固定修订包如何交接行身份、结果值与全部权重

`Calibrator.fit` 不复制或重排行，而是直接持有不可变的源
`ReplicateSnapshot`。包暴露的 `rows`、`row_ids`、`revision` 是源对象的只读
投影，`weights` 与 `replicate_weights` 是新求出的冻结元组，位置始终由源
`row_ids` 解释。随后的行账 `apply` 发布新修订，但旧快照对象及其数据包不被
改动，因此旧报告数值保持不变（由 `test_bundle_survives_ledger_edit_and_reorder`
和 demo 的 `old_bundle_after_edit` 覆盖）。零支持行（复制基础权重为 0）
不进入求解器，输出位置保持 0。

## 求解算法

1. **Phase I 可行性**：为每条等式加人工松弛变量，把原变量先固定在下界并令
   松弛量吸收差额，得到一个严格可行起点；最小化松弛平方（原变量仅带
   `1e-10` 的严格凸正则项）。松弛无法收敛到零时，再用枚举界配置的组合可行性
   兜底处理退化顶点；两者都失败才判定不可行。矩阵奇异或迭代停滞本身从不作为
   不可行依据。
2. **Phase II 原始积极集**：自由变量通过 KKT 块求等式约束下的目标点，按
   最小步长触界并加入工作集；步长为零时用归约乘子判断是否释放某个界约束，
   每次释放都重算投影方向确认可行下降，冗余等式产生的零残差退化也能释放。
3. **冗余与降秩**：每次构造 KKT 前用带行列标签的高斯消元选取极大线性无关
   等式行，允许重复控制信息、零列和删组后进一步降秩；恰在上下界的解照常接纳。
4. **容差**：求解器内部使用收紧的工作容差，`Calibrator` 再按 CALIBRATION.md
   的归一化残差 `1e-8`、界违反 `1e-8` 复核，失败按不可行处理。

## 独立最优性核验（不以生产求解器自证）

`tests/test_calibration.py` 的 `exact_reference` 完全用 `Fraction` 精确算术
独立实现：枚举每个正支持变量处于下界、上界或自由的全部划分，自由变量经精确
KKT 求该面上的驻点，校验全部原始等式（含被归约丢弃的冗余行）与界，取目标
最小者为参照最优值；无可行划分才判不可行。生产权重与参照逐值比较，并比较
目标 `Q`。该文件还覆盖：冗余/零控制列、上下界最优、零支持、总体可行但复制
不可行、域复制分母为零、更正与重排不变性，以及 60 行 4 控制 8 复制的全量
调用和独立方差公式复算。随机小例与多轮随机压力测试均与精确参照一致。
