# 有界校准复制分析

## 共同业务结果与现状

人口控制总量用于校准调查权重。点估计和删组复制方差必须对应同一校准方法、同一固定行账修订；不能仅校准总体，再拿未校准复制权重计算方差。

当前 README.md 中的 Ledger、LedgerSnapshot、Row、DeleteGroups、ReplicateSnapshot、Reporter.raw_ratios 及原有行为全部保留。现有实现只计算未校准总体比例，尚无本文件要求的校准或复制方差路径。两项新增能力为：Calibrator.fit 生成完整校准复制包；Reporter.ratios 消费该包计算分域比例与复制方差。包的内部类和求解组织由实现者决定。

## 调用入口和输入

在 survey 包顶层导出 Calibrator 和 CalibrationInfeasible。Calibrator.fit(snapshot, controls, targets, lower, upper) 可直接按此形式调用，无需构造实例。

- snapshot 必须为现有 ReplicateSnapshot，行数为 0 至 60；复制数沿用其 1 至 8 个限制。
- controls 是 1 至 4 个互异控制名称的序列，每个名称都须存在于 snapshot.control_names；targets 是按 controls 顺序对应的等长目标总量序列。名称通过原快照的列名解释，允许选择子集及重排，不能把所选名称的位置当成原控制向量列号。
- 总体基础权重沿用 snapshot.weights；复制基础权重沿用 snapshot.replicate_weights[r]。本次 fit 的正基础或复制权重须在 [0.1,100] 内，零复制权重合法。现有删组设计可能生成超过 100 的权重，其原行为不变，但该输入不属于本次 fit 的接纳范围。
- targets、lower、upper 只接纳有限 int/float，排除 bool；0 < lower <= 1 <= upper <= 10。目标可为负数。行控制和结果值仍遵守 Row 的原范围；复制系数沿用原快照的非负有限值，允许为零。
- 非法类型、维度、名称、数值范围或界倍数抛 ValueError；普通 Python 非可迭代参数等类型误用可沿用 README.md 的 TypeError 约定。合法但不能满足校准等式和权重界的输入抛 CalibrationInfeasible，该异常继承 ValueError。输入迭代器自身异常原样传播。
- 所有入口均不得更改传入快照、域集合或其他输入。fit 不修改行账修订。

## A：总体和全部复制分别校准

对总体及每个复制，分别令 d_i 为该样本的原权重，x_ij 为 controls[j] 对应的第 i 行控制值，t_j 为 targets[j]。求权重 w：

1. d_i=0 时 w_i=0；d_i>0 时 lower*d_i <= w_i <= upper*d_i。
2. 对每个所选控制量 j，sum_i(x_ij*w_i)=t_j。
3. 在以上约束下最小化 Q=sum_{i:d_i>0}((w_i-d_i)^2/d_i)。

复制使用自己的 d_i 和与总体相同的 targets、lower、upper，不能由总体校准权重缩放或删行代替独立校准。Q 在正支持上严格凸，所求权重唯一；零支持行仍保留原位置。控制允许重复线性信息、零列及冗余约束，删组后可能进一步降秩。恰在上下界的可行解必须接纳，冗余不等于不可行；空支持且全部目标为零时唯一权重全零，目标不全零则不可行。

总体或任一复制不可行时，整次 fit 抛 CalibrationInfeasible，不返回部分包。成功包须完整固定原快照的行身份、结果值、修订及全部校准权重，不受后续行账更正、删除或重排影响。

返回包至少公开下列只读属性，所有序列均不可变且不保留可变输入别名：

- revision、rows、row_ids：原快照的修订、原有 Row 序列及对应身份。
- weights：总体校准权重，按 rows 顺序。
- replicate_ids、replicate_weights、coefficients：复制名称及原顺序、全部校准复制权重矩阵、原系数。replicate_weights[r][i] 的两维顺序均与原 snapshot 对应。

无需指定包的具体类名，也无需跨进程保存；不可把修订号相同但来自不同行账的行和权重混用。允许复用原本已不可变的源对象。

### 数值判据

验收对每个总体或复制独立核对：

- 每个控制等式的归一化残差不超过 1e-8。归一化分母为 max(1, abs(t_j), sum_i(abs(x_ij)*d_i))。
- 每行权重界的绝对违反量不超过 1e-8；零支持权重保持为零。
- 相对于独立参照最优值 Q*，abs(Q-Q*) <= 1e-7*max(1,abs(Q*))。

验收含边界可行、约束冗余及明确不可行输入，不以小于上述容差的可行性差异判成败，也不依靠近乎不可分辨的秩或最优值扰动。不能把普通迭代停滞或矩阵奇异本身当成不可行证明。允许逐复制全量求解，不要求增量、并行或特定算法；没有额外的数据持久化目标。

## B：从校准包计算分域比例及复制方差

Reporter.ratios(bundle, domains) 作为现有 Reporter 的新增入口，可按此形式直接调用；bundle 接受 fit 生成的完整包，其他对象抛 ValueError。domains 沿用 raw_ratios 的域名到行 ID 可迭代对象映射、集合语义、名称规则及映射顺序。域可以重叠或为空，未知行 ID 使整个调用抛 ValueError。

返回只读报告，公开 revision 和按域映射顺序排列的不可变 estimates。每项公开：

- name、row_ids、count：沿用 DomainRatio 的语义，row_ids 按 ID 排序。
- numerator、denominator、ratio：使用包中的总体校准权重及对应 Row.value；前两项为域内 sum(w_i*y_i)、sum(w_i)，ratio 为两者之比。
- replicate_ratios：按 replicate_ids 原顺序保留全部复制的域比例，各自使用对应校准复制权重；该复制域分母为零时该项为 None。
- variance：theta=ratio，theta_r=replicate_ratios[r]，值为 sum_r(coefficients[r]*(theta_r-theta)^2)。
- missing_replicate_ids：域分母为零的复制 ID 元组，顺序与包一致。

总体分母为零时 ratio 和 variance 均为 None；replicate_ratios 与 missing_replicate_ids 仍按全部复制实际分母填写。仅复制分母为零时保留总体比例，variance 为 None，并保留全部缺失复制编号；即使对应系数为零，也不得静默删除该复制或重新归一化其他系数。无缺失复制且总体比例存在时使用规定公式，不另做自由度修正或开平方。空域 count、numerator、denominator 均为零。

域分子、分母、有效比例和方差相对独立参照的误差须不超过 1e-8*max(1,abs(参照值))。报告不应随调用方随后修改域映射而改变。允许逐域扫描整个包。raw_ratios 必须继续使用原基础权重，保持其未校准语义及原返回结构。

## 完整调用样例

以下示例使用现有行账及复制设计，期望新增入口完成后可运行：

```python
from survey import Ledger, Row, DeleteGroups, Calibrator, Reporter

ledger = Ledger(("population",), [
    Row("a", 1, (1,), 0, "A"),
    Row("b", 1, (1,), 1, "B"),
])
snapshot = DeleteGroups(("A", "B")).build(ledger.snapshot())
bundle = Calibrator.fit(snapshot, ("population",), (3,), 0.5, 2)
report = Reporter.ratios(bundle, {"all": {"a", "b"}, "a-only": {"a"}, "empty": set()})
```

总体权重应为 (1.5,1.5)，两个复制依次为 (0,3)、(3,0)，复制系数均为 0.5。all 的比例为 0.5，replicate_ratios 为 (1,0)，variance 为 0.25，missing_replicate_ids 为空元组。a-only 的总体比例为 0，replicate_ratios 为 (None,0)，variance 为 None，缺失编号为 ("A",)。empty 的总体比例与方差均为 None，全部复制比例为 None，缺失编号为 ("A","B")。

随后 ledger.apply([Row("a",1,(1,),1,"A")], order=["b","a"])，以上旧包的所有数值和报告保持原样。由新快照重新 fit 后全域比例为 1，复制比例均为 1，方差为 0。按稳定行 ID 比较，单纯重排输入不改变对应校准权重或域结果。

### 能区分“只校准总体”的联合例

四行按 ID a、b、c、d 排列，基础权重均为 1，控制名为 ("n","x")，控制向量依次为 (1,0)、(1,1)、(1,2)、(1,3)，结果值为 (0,0,0,1)，组为 (A,B,A,B)。使用 DeleteGroups(("A","B"))，目标为 (4,6)，lower=0.5、upper=2。

总体校准权重为 (1,1,1,1)，A、B 复制的校准权重分别为 (0,3,0,1)、(1,0,3,0)。全域比例为 1/4，复制比例为 (1/4,0)，方差为 1/32。若仅校准总体而沿用未校准复制权重，将错误得到方差 1/16。

### 边界和不可行例

两行基础权重均为 1、单控制分别为 1、2，组分别为 A、B，目标为 5，lower=0.5、upper=2。总体可行，但删除 B 后只剩第一行，复制基础权重为 (2,0)，最大控制总量为 4，因此整次 fit 必须抛 CalibrationInfeasible。

一个外部 ReplicateSnapshot 的源行控制为 (1,2)、基础权重为 1，只有一个复制且复制基础权重为 2；目标 (2,4)、lower=0.5、upper=2 时总体与复制均可行。总体达到自身上界 2，复制也得到权重 2；两条控制约束线性相关，不能因冗余或触界拒绝。目标改为 (2,5) 则不可行。

## 交付与核验

只完成校准包而继续调用 raw_ratios，没有接通校准后方差；只新增方差循环而沿用原复制权重，也没有完成共同业务结果。交付两项联动实现、简短的包交接及求解正确性说明，并提供独立小规模最优性核验，覆盖控制冗余、零支持、上下界最优、总体可行但复制不可行、域复制分母为零以及行账更正与重排。

完整核验应包含 60 行、4 个控制量和 8 个复制样本的调用；小例独立参照不能以生产求解器自证。无需增加其他估计量、标准误输出、置信区间、抽样设计、界面或缓存。

保留 tests/test_survey.py 的全部字节及原公开接口。原测试入口为：
```sh
python3.11 -B -m unittest discover -s tests -p 'test_*.py' -v
python3.11 -B -m survey demo
```

沿现有 .pointwise-runtime.json 接通完整 verify 和真实 demo：verify 应执行全部原回归及新增独立核验，demo 应展示校准后的复制方差、复制域缺失和更正后旧包不变。保留原 demo 命令及其原有输出语义，可新增命令用于校准演示；运行清单可以指向同时调用二者的包装入口。清单的 version、kind、setup、verify、demo 字段结构继续沿用现有文件。新增核验及演示须真实调用公开业务入口，不能只打印预写结果。
