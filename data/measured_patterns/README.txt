README — 上行业务1970TXT 数据格式说明
=====================================

一、目录内容
------------
本目录是「上行业务1970」目录下全部 78 个 .xlsx 波束方向图文件的 .txt 版本，
文件名一一对应（仅扩展名不同），例如：

    上行业务1970/pattern_Pitch_0_Azimuth_0.xlsx
 -> 上行业务1970TXT/pattern_Pitch_0_Azimuth_0.txt

转换程序：工作区根目录 convert_xlsx_to_txt.py（python + openpyxl，
转换日期 2026-10-20）。转换结果已逐文件校验：矩阵结构、轴标签、
采样增益值与源 xlsx 完全一致（78/78 通过）。

二、文件名含义
--------------
    pattern_Pitch_<p>_Azimuth_<a>.txt

  <p> = 俯仰角 Pitch，取值 0, 10, 20, 30, 40, 50（度）
  <a> = 方位角 Azimuth，取值 0, 30, 60, ..., 360（度）

共 6 × 13 = 78 个文件，为完整扫描网格，无缺失。
业务类型：上行业务；频率：1970 MHz（取自源目录名「上行业务1970」）。

三、txt 文件内部格式
--------------------
  编码：UTF-8；分隔符：Tab；行尾：\n；共 97 行。

  第 1~5 行   ：以 # 开头的说明头（来源、坐标轴、数值含义、版式），
                纯数值解析时请跳过。
  第 6 行     ：列头行，格式为  theta\phi + 0  1  2  ...  359
                即第 1 个字段是占位标签 theta\phi，其后 360 个字段是
                phi 方位角列标签（度）。
  第 7~97 行  ：数据行，每行 361 个字段：
                第 1 个字段 = theta 离轴角行标签（度，0~90，整数，1° 步进）；
                其后 360 个字段 = 该 (theta, phi) 下的波束增益值。

  与源 xlsx 的对应关系（逐格一致）：
      txt 第 6 行    = xlsx 第 1 行（列头 phi）
      txt 第 7~97 行 = xlsx 第 2~92 行
      txt 每行第 1 个字段 = xlsx 该行第 1 列（行头 theta）

四、数值含义
------------
  单位：dB（归一化，0 dB = 该方向图峰值；依据项目根目录 ReadMe.txt
  「所提供波束数据均按照归一化数值提供」）。
  取值范围：-65 ~ 0 dB。
  小数位：最多 2 位（与源 xlsx 一致，转换无损）。
  注意：theta = 0 行（波束轴向上）在所有文件中均为 0 dB，
  即峰值增益位于离轴角 0 处，属正常现象。

五、读取示例（Python / numpy）
------------------------------
  import numpy as np
  f = 'pattern_Pitch_0_Azimuth_0.txt'
  # 前 5 行 # 注释；第 6 行列头第 1 字段是 theta\phi 文本，用 usecols 跳过
  phi = np.loadtxt(f, skiprows=5, max_rows=1, delimiter='\t',
                   usecols=range(1, 361))            # 方位角 0..359，共 360 个
  grid = np.loadtxt(f, skiprows=6, delimiter='\t')   # 91 x 361（仅数据行）
  theta = grid[:, 0]                                 # 离轴角 0..90，共 91 个
  gain = grid[:, 1:]                                 # 91 x 360，归一化增益 dB

  # pandas 等价：pd.read_csv(f, sep='\t', skiprows=6, header=None)

六、快速自检
------------
  wc -l pattern_*.txt              # 每个文件应输出 97
  awk -F'\t' 'FNR>5 && NF!=361' pattern_*.txt   # 第 6 行起字段数不等于 361 的行（应无输出）

七、备注
--------
* 源 xlsx 每个文件仅含一个工作表 Sheet1，无公式、无合并单元格，
  数据区无空值、无非数值字符。
* 如需重新生成或转换其他目录（信令展宽 1970 / 下行 2160 等），运行：
      python convert_xlsx_to_txt.py <源目录> <输出目录>
