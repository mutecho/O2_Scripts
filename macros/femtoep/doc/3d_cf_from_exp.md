## 说明程序
3D_cf_from_exp.cpp，类型为ROOT宏

## 功能

1、从femtoscopy实验数据中提取3Dfemto相关的关键观测量（根据cent，pairmt，phi-psi区分的3D correlationfunction，相应的1d投影等），存入中间结果文件

2、而后对其进行分区间levy-stable拟合，得到相关6个维度的$R^2$拟合值和对应的误差；使用$R^2$和phi-psi的关系绘制成曲线，与拟合结果一同保存在拟合结果文件中


## 输入输出格式和要求

### 输入
ROOT文件，sameevent和mixedevent的实验数据保存在各自的ndhisto中

### 输出
#### 中间结果文件：
ROOT文件：根据cent，pairmt，phi-psi区分的3D correlationfunction，相应的1d投影等和SE/ME的raw 3Dhisto。其中1d投影应满足CF的基本要求（如高q区间值靠近1）

#### 最终文件：
ROOT文件：根据cent，pairmt，phi-psi区分的对3D correlationfunction的拟合结果（以及对应的参数值，误差等），相应拟合结果在1d上的投影（和1dCF投影一起展现）；以及$R^2$和phi-psi的关系份区间的曲线。拟合应该与数据贴合良好。

TXT：各区间拟合得到的参数，误差，chi^2,ndf,拟合质量等等。