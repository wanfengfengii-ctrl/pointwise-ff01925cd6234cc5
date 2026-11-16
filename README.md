# Survey Ledger

Python 3.11 标准库调查分析组件。行账接纳原子更正批次，分析固定在不可变修订上；删组设计公开复制权重及系数，报告按稳定行 ID 计算未校准的分域加权比例。

## 运行

在项目根目录运行，无安装步骤、无网络访问，也不写入业务文件：

```sh
python3.11 -B -m unittest discover -s tests -p 'test_*.py' -v
python3.11 -B -m survey demo
```

运行清单为 .pointwise-runtime.json。demo 使用 survey/sample.py 中四户有限样例，输出更正前报告、更正后旧快照报告、新报告和新快照的全部复制权重。原回归包含独立 Fraction 核账、固定随机种子以及 60 行、4 控制量、8 复制样本的分析。

## 行账

`Row(row_id, weight, controls, value, group)` 是不可变记录：

- row_id 与 group 为非空字符串；row_id 在一个修订中唯一。
- weight 是 [0.1,100] 内的基础权重。
- controls 是 1 至 4 项控制向量，value 是调查结果值，两者元素均在 [-100,100]。
- 数值只接受有限 int/float，排除 bool。输入控制向量被复制为不可变序列。

`Ledger(control_names, rows=())` 创建修订 0；control_names 是与每行控制向量逐列对应的 1 至 4 个唯一非空名称。空行账合法。

`ledger.apply(upserts=(), deletes=(), order=None)` 完整验证后发布一个新修订，并返回该 LedgerSnapshot。每次成功调用都加一，包括空批次。已有行替换后保持位置；新行按 upserts 顺序追加；删除只接受已有 ID。order 如提供，必须是最终全部 ID 的完整排列。重复更新、重复删除、同 ID 同时更新及删除、控制维度不符或排列不完整，均以 ValueError 拒绝整个批次。输入迭代器自身异常原样传播；所有失败均保持原快照对象及修订不变。

`ledger.snapshot()` 返回当前不可变 LedgerSnapshot，公开 revision、control_names、rows、row_ids。后续行账更正、删除或重排不改变已有快照。修订号仅在所属行账内有意义。

## 复制设计

`DeleteGroups(groups).build(ledger_snapshot)` 返回不可变 ReplicateSnapshot。groups 为 2 至 8 个唯一组名，给出公开复制顺序。快照必须恰好包含这些组且每组至少一行；空行账或缺组快照不能生成此设计。

设组数为 G，复制 r 删除 groups[r] 中的行：被删除行权重为 0，其余行为原基础权重乘 G/(G-1)。每个复制的系数为 (G-1)/G。这里采用非分层、等组删组设计；调用方负责选择符合调查设计的组划分。

ReplicateSnapshot 公开：

- source：原 LedgerSnapshot；revision、rows、row_ids、control_names 是它的只读投影。
- weights：按 rows 顺序排列的总体基础权重。
- replicate_ids：复制的稳定名称及顺序。
- replicate_weights[r][i]：第 r 个复制中第 i 行的权重，含被删除的零权重行。
- coefficients[r]：与该复制同位置的非负系数。

权重序列中的行位置始终由该对象的 row_ids 解释，不能与另一个修订的位置混用。复制及系数被完整冻结，后续更正不影响旧对象。

也可直接调用 `ReplicateSnapshot(source, replicate_ids, replicate_weights, coefficients)` 接纳外部已给定的复制设计。支持 1 至 8 个唯一复制名称，矩阵形状必须完整；权重及系数为有限非负数，均被复制为不可变序列。直接构造器只验证形状和数值，不判断外部设计的统计适用性。生成器可能产生超过 100 的复制权重，这是基础权重放大的正常结果。

## 分域报告

`Reporter.raw_ratios(replicate_snapshot, domains)` 接受域名到行 ID 可迭代对象的映射。域名为非空字符串；重复成员按集合解释，域可以重叠。未知行 ID 使整个调用抛 ValueError。成员集合按 ID 解释，独立于行顺序。

返回不可变 RawReport，含 revision 和按域映射顺序排列的 estimates。每个 DomainRatio 含 name、按 ID 排序的 row_ids、count、numerator、denominator、ratio。numerator 为域内 Σ(weight × value)，denominator 为 Σweight，采用 math.fsum；ratio 为两者之比。空域的分子、分母、行数为 0，比例为 None。

当前报告提供总体基础权重的点估计。复制权重与系数作为公开调查设计数据保存。

## 示例

```python
from survey import Ledger, Row, DeleteGroups, Reporter

ledger = Ledger(("population",), [
    Row("a", 1, (1,), 0, "g1"),
    Row("b", 1, (1,), 1, "g2"),
])
frozen = DeleteGroups(("g1", "g2")).build(ledger.snapshot())
report = Reporter.raw_ratios(frozen, {"all": {"a", "b"}, "empty": set()})
assert report.estimates[0].ratio == 0.5
assert report.estimates[1].ratio is None
assert frozen.replicate_weights == ((0.0, 2.0), (2.0, 0.0))
```

所有公开序列输出均不可变。领域约束错误抛 ValueError；非可迭代参数等 Python 类型误用可能抛 TypeError。组件不创建工作文件，不启动后台进程。报告对每个域扫描当前固定快照，行账更新和复制生成也采用完整构建，行为确定且适用于重复核账。
