# `3d_cf_from_exp.cpp` 最小修复清单

面向实现 agent 的执行说明。本文档只覆盖以下 3 项最小修复，其余 review 中提到的问题本轮一律忽略，不扩 scope。

目标文件：
- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp`

参考文档：
- `/Users/allenzhou/ALICE/scripts/macros/femtoep/doc/3d_cf_from_exp.md`

## 本轮范围

只做以下 3 项：

1. 修复 PML 拟合对空 bin 的处理
2. 给 full 3D Levy fit 的 `R^2` 非对角元加最小物理约束
3. 把“每个 slice / histogram 是否重新打开一次输出 ROOT 文件”做成开关

明确不做：

- 不修改 CF 的整体归一化定义
- 不修改 Coulomb / Gamow 近似
- 不修改 `phi` 谐波形式
- 不重构 legacy 1D helper
- 不顺手修其他 review 问题

## 实施顺序

建议顺序：

1. 先修 PML 空 bin
2. 再加 `R^2` 约束
3. 最后做输出文件 reopen 开关

这样便于先保证拟合统计定义正确，再处理 full model 的物理参数空间，最后再动 I/O 路径。

## 任务 1: 修复 PML 空 bin 问题

### 背景

当前 `Levy3DPMLFCN()` 在以下逻辑中直接跳过了任一侧为 0 的 bin：

- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:1553`
- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:1557`

`CountPMLUsableBins()` 也只统计 `SE>0 && ME>0` 的 bin：

- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:1588`
- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:1608`

这会让单边零计数 bin 完全不进入 likelihood，也会让 `ndf` 偏小。

### 最小修复要求

- 只在 `sameCounts == 0 && mixedCounts == 0` 时跳过该 bin。
- 若仅一边为 0，该 bin 仍然必须贡献到 `-2 ln L`。
- 保留当前两边都大于 0 时的公式。
- 对单边零计数情况，使用稳定的极限形式，并优先使用 `std::log1p(...)`。
- `CountPMLUsableBins()` 改为统计 `sameCounts + mixedCounts > 0` 的 bin。
- `ndf` 的计算逻辑保持原样，但其输入 bin 计数要与新的 FCN 逻辑一致。

### 建议实现方式

在 `Levy3DPMLFCN()` 中分 4 支：

1. `sameCounts < 0 || mixedCounts < 0`：直接报错/罚掉该参数点
2. `sameCounts == 0 && mixedCounts == 0`：`continue`
3. `sameCounts == 0 && mixedCounts > 0`：使用单边零计数极限式
4. `mixedCounts == 0 && sameCounts > 0`：使用单边零计数极限式
5. 两边都大于 0：走现有公式

### 验收标准

- PML 路径下，单边零计数 bin 不再被直接忽略。
- `CountPMLUsableBins()` 与 FCN 的有效 bin 定义一致。
- 在低统计输入上，PML 不会因为大量 `SE==0` bin 被跳过而明显高估相关强度。
- 不改变 `fitUsePML=false` 时的行为。

## 任务 2: 给 `R^2` 非对角元加最小约束

### 背景

当前 full model 中：

- 非对角元参数在 `BuildFullLevyFitFunction()` 中是自由的
- `EvaluateFullLevyCF()` 对负的二次型只做了 `max(argument, 0.0)` 裁剪

相关位置：

- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:892`
- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:897`
- `/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:1026`

这会让 fitter 进入非物理参数区，然后靠裁剪继续跑。

### 本轮目标

本轮不做大改，不改成 Cholesky 参数化；只做“最小可用”的物理约束，确保 full fit 不再静默接受明显非物理的 `R^2` 组合。

### 最小修复要求

- 新增一个 helper，用于检查 full model 的 `R^2` 矩阵是否为半正定。
- 检查对象为：
  - `Rout2`
  - `Rside2`
  - `Rlong2`
  - `Routside2`
  - `Routlong2`
  - `Rsidelong2`
- 至少检查以下量：
  - 三个对角元 `>= 0`
  - 三个 `2x2` 主子式 `>= 0`
  - 整个 `3x3` 行列式 `>= 0`
- `EvaluateFullLevyCF()` 不允许再单纯依赖 `max(argument, 0.0)` 掩盖参数非法。
- PML 路径下，若 `R^2` 矩阵非法，FCN 应直接返回大罚值。
- `chi2` 路径下，也要让非法参数点变成明显差的拟合点，而不是被静默裁剪成可继续拟合的形状。

### 建议实现方式

- 增加一个类似 `IsFullR2MatrixPositiveSemiDefinite(...)` 的 helper。
- 在 `EvaluateFullLevyCF()` 开头先检查矩阵合法性。
- 对于仅由浮点舍入导致的极小负值，可以给一个很小的容差，例如 `-1e-10` 级别；不要用宽松容差掩盖真实非法点。
- 仅在矩阵已通过半正定检查后，`argument` 才允许做接近 0 的数值保护。
- 保持 diagonal model 不受影响。

### 验收标准

- full model 不再依赖 `max(argument, 0.0)` 隐藏明显非物理参数点。
- PML 和 `chi2` 两条路径都不会接受明显非半正定的 `R^2` 组合。
- diagonal fit 的数值结果不受影响。
- 不引入 full model 之外的拟合接口变化。

## 任务 3: 把输出文件 reopen 行为做成开关

### 背景

当前 3D CF 构造和拟合写出时，都会在每个 slice / histogram 内部重新打开一次输出 ROOT 文件：

- CF 构造阶段：`/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:605`
- 拟合写出阶段：`/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp:1862`

这会造成明显 I/O 开销。目标不是强行改默认行为，而是把它做成可切换。

### 最小修复要求

- 增加一个布尔开关，控制“是否在每个 slice / histogram 写出时重新打开 ROOT 文件”。
- 该开关应同时覆盖：
  - `CFCalc3D()` 写出路径
  - 拟合结果写出路径
- 默认值保持当前行为兼容，建议默认 `true`。
- 当开关为 `false` 时：
  - CF 构造阶段的输出文件应在单次 `CFCalc3D()` 调用内尽量复用
  - 拟合写出阶段的输出文件应在单次 fit driver 调用内尽量复用

### 建议实现方式

- 在顶层 `_3d_cf_from_exp()` 增加一个显式配置项，例如：
  - `bool reopenOutputFilePerSlice = true;`
- 向下透传到：
  - `CFCalc3D(...)`
  - `FitCF3DWithSelectedBins(...)`
  - 如有必要，再透传到 `FitAndWriteSingleCFHistogram(...)`
- 为避免改动过大，可以保留现有“内部自己开文件”的路径作为默认分支。
- 新分支中只需要做到“单次 stage 复用一个已打开的 `TFile*`”；不要求进一步做目录缓存或更大规模重构。

### 验收标准

- 开关为 `true` 时，输出行为与当前版本兼容。
- 开关为 `false` 时，不再对每个 slice / histogram 单独 reopen 输出文件。
- 两种模式生成的对象命名和目录布局保持一致。
- 不改变物理结果，只改变写出方式与 I/O 性能。

## 建议交付内容

本轮交付应至少包含：

- 代码改动
- 一次最小编译/加载验证
- 一段简短实现说明，说明：
  - PML 空 bin 现在如何处理
  - full model 的 `R^2` 约束如何实现
  - reopen 开关名称、默认值、作用范围

## 最低验证要求

至少完成以下验证：

1. ROOT/Cling 能正常加载修改后的宏
2. `fitUsePML=true` 时能跑通，不因单边零计数 bin 崩掉
3. full fit 在明显非法的非对角参数组合下会被拒绝/罚掉，而不是被静默裁剪
4. `reopenOutputFilePerSlice=true/false` 两种模式都能完成一次写出流程

## 备注

如果实现过程中发现第 2 项必须改成更大的参数化重构才能稳定工作，需要先停下来汇报，不要在本轮自行扩大到“整体重写 full Levy 参数化”。
