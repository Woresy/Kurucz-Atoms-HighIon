"""Generate a Chinese HTML explanation of the Ti II .trans workflow."""

from pathlib import Path


OUT = Path("reports/meeting-2026-07-17/TiII_lines_agafgf_states_to_trans_workflow.html")


HTML = r"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>Ti II：从 .lines、.agafgf 和 .states 生成 .trans 的完整流程</title>
<style>
@page { size:A4; margin:17mm 16mm 18mm; @bottom-center { content:counter(page); font-size:9pt; color:#64748b; } }
body { max-width:1120px; margin:0 auto; padding:28px; font-family:"Noto Sans CJK SC","Microsoft YaHei",sans-serif; color:#172033; font-size:15px; line-height:1.72; }
h1 { color:#123b5d; text-align:center; font-size:30px; margin:30px 0 8px; }
.subtitle { text-align:center; color:#607184; font-size:17px; margin-bottom:34px; }
h2 { color:#145374; border-bottom:2px solid #4aa3a2; padding-bottom:6px; margin-top:34px; }
h3 { color:#17627d; margin-top:26px; }
.key { background:#e8f4ff; border:2px solid #3b82a0; border-radius:7px; padding:16px 20px; font-size:17px; }
.note { background:#f2f7fa; border-left:5px solid #4aa3a2; padding:13px 17px; }
.warning { background:#fff6e8; border-left:5px solid #d88920; padding:13px 17px; }
.ok { background:#edf8f0; border-left:5px solid #3b8c58; padding:13px 17px; }
.flow { text-align:center; font-weight:600; line-height:2.15; background:#f8fafc; border:1px solid #a9b8c6; border-radius:6px; padding:16px; }
.arrow { color:#237a83; font-size:23px; }
table { width:100%; border-collapse:collapse; margin:15px 0 22px; font-size:14px; }
th { background:#174f73; color:white; }
th,td { border:1px solid #aebdca; padding:9px 10px; text-align:left; vertical-align:top; }
tr:nth-child(even) td { background:#f7fafc; }
pre { background:#18222d; color:#edf4f7; padding:15px 18px; border-radius:6px; overflow:auto; font:13px/1.65 "DejaVu Sans Mono",monospace; }
code { font-family:"DejaVu Sans Mono",monospace; background:#eef2f6; padding:1px 4px; border-radius:3px; }
.equation { text-align:center; font:17px "DejaVu Serif",serif; padding:9px; }
.small { color:#607184; font-size:13px; }
ol li, ul li { margin:5px 0; }
</style></head><body>
<h1>Ti II：从 .lines、.agafgf 和 .states 生成 .trans</h1>
<div class="subtitle">以第一条有效跃迁为例的完整数据流说明</div>

<div class="key"><b>先记住三个文件的分工：</b><br>
<code>.lines</code> 回答“哪两个能级发生跃迁”；<code>.states</code> 回答“这两个能级在 ExoMol 中编号多少”；<code>.agafgf</code> 回答“这次跃迁的 A 系数和波数是多少”。最后把四个结果写入 <code>.trans</code>。</div>

<h2>1. 总体数据流</h2>
<div class="flow">
.lines：E₁、J₁、label₁；E₂、J₂、label₂<br>
<span class="arrow">↓ 用 (|E|, J, 原始 label) 查找能级</span><br>
中间 states 映射表：得到 state ID 1 和 4413<br><br>
.agafgf：读取同一条跃迁的波数 ν̃ 和 log(A)<br>
<span class="arrow">↓ 计算 A = 10<sup>log(A)</sup></span><br>
将 state ID、A 和波数合并<br>
<span class="arrow">↓ 按 ExoMol 格式输出</span><br>
.trans：4413&nbsp;&nbsp;1&nbsp;&nbsp;1.0162e+03&nbsp;&nbsp;2.174206e+05
</div>

<h2>2. 第一步：从 .lines 读取跃迁的两个端点</h2>
<p><code>gf2201.lines</code> 是固定宽度文件，不是普通的空格分隔表。项目使用前 80 个字符中的 9 个字段：</p>
<table><tr><th>字符范围（Python）</th><th>字段</th><th>含义</th><th>是否用于 state ID 匹配</th></tr>
<tr><td>0:11</td><td>wl</td><td>波长（nm）</td><td>否</td></tr>
<tr><td>11:18</td><td>log_gf</td><td>加权振子强度的常用对数</td><td>否</td></tr>
<tr><td>18:24</td><td>ele</td><td>元素与电离阶段编码；2201 表示 Ti II</td><td>否</td></tr>
<tr><td>24:36</td><td>E1</td><td>第一个能级的能量</td><td>是</td></tr>
<tr><td>36:41</td><td>J1</td><td>第一个能级的总角动量</td><td>是</td></tr>
<tr><td>41:52</td><td>label1</td><td>第一个能级的 Kurucz 原始标签</td><td>是</td></tr>
<tr><td>52:64</td><td>E2</td><td>第二个能级的能量</td><td>是</td></tr>
<tr><td>64:69</td><td>J2</td><td>第二个能级的总角动量</td><td>是</td></tr>
<tr><td>69:80</td><td>label2</td><td>第二个能级的 Kurucz 原始标签</td><td>是</td></tr></table>

<p>第一条有效记录中，与端点映射有关的内容为：</p>
<pre>端点 1：E1 =       0.000 cm⁻¹   J1 = 1.5   label1 = (3F)4s a4F
端点 2：E2 = -217420.587 cm⁻¹   J2 = 2.5   label2 = m(1S)2F344</pre>
<div class="note"><b>为什么要同时使用 E、J 和 label？</b> 只用能量或 J 可能遇到重复、简并或近似相同的能级。三元组 <code>(|E|, J, label)</code> 才是项目用于定位原始能级的键。</div>

<h2>3. 第二步：从原始能级映射到 .states 的 state ID</h2>
<p>最终 <code>.states</code> 是能级表。每个能级有一个唯一的整数编号，称为 state ID。编号是在能级清理、按能量排序后从 1 开始分配的。因此，<b>1 和 4413 是能级编号，不是源文件行号，也不是能量。</b></p>

<h3>3.1 能量符号的处理</h3>
<p>Kurucz 数据可能以负号标记计算能级。生成 states 和执行匹配时使用能量绝对值：</p>
<div class="equation">E1: |0.000| = 0.000；　E2: |−217420.587| = 217420.587 cm⁻¹</div>
<p>负号所表达的来源类别不会完全丢失：最终 states 的 <code>Abbr</code> 列用 <code>CA</code> 标识计算能级，用 <code>NI</code> 标识观测能级。</p>

<h3>3.2 使用保留原始 label 的中间映射表</h3>
<p>state ID 的查找不是直接拿 `.lines` 的紧凑标签与最终 `.states` 的展示字符串比较。程序使用一份保留原始 label 的中间 states 表：</p>
<table><tr><th>原始匹配键</th><th>匹配结果</th></tr>
<tr><td><code>(0.000, 1.5, "(3F)4s a4F")</code></td><td>state ID = <b>1</b></td></tr>
<tr><td><code>(217420.587, 2.5, "m(1S)2F344")</code></td><td>state ID = <b>4413</b></td></tr></table>

<h3>3.3 为什么 4413 的最终 label 看起来不同？</h3>
<p>原始紧凑标签在输出 `.states` 时被解析为 Configuration 和 Term：</p>
<pre>原始 label：m(1S)2F344
             │ │   └── 344：Kurucz 标签中的序号信息，最终不输出
             │ └────── 2F：谱项，写入 Term 列
             └──────── m(1S)：组态部分；m 通过 Ti II 映射表展开为 s211f

最终显示：Configuration = s211f(1S)
          Term          = 2F</pre>
<p>因此下面两种写法指向同一个能级：</p>
<table><tr><th>阶段</th><th>能量</th><th>J</th><th>标签形式</th><th>ID</th></tr>
<tr><td>.lines / 中间映射</td><td>217420.587</td><td>2.5</td><td>m(1S)2F344</td><td>4413</td></tr>
<tr><td>最终 .states</td><td>217420.587</td><td>2.5</td><td>s211f(1S) + 2F</td><td>4413</td></tr></table>

<h2>4. 第三步：从 .agafgf 取得波数和 Einstein A</h2>
<p><code>gf2201.agafgf</code> 的物理第一行是文字表头，表头后的第一行才是第一条有效数值记录。对于本例，该记录给出：</p>
<pre>wn     = 217420.587 cm⁻¹
log(A) = 3.007...</pre>
<p>文件保存的是 <code>log10(A)</code>，所以要转换为普通的 A 系数：</p>
<div class="equation">A = 10<sup>log(A)</sup> ≈ 1.0162 × 10³ s⁻¹</div>
<table><tr><th>.agafgf 字段</th><th>处理</th><th>写入 .trans</th></tr>
<tr><td>wn = 217420.587</td><td>取波数大小</td><td>2.174206e+05 cm⁻¹</td></tr>
<tr><td>log(A) = 3.007...</td><td>A = 10<sup>log(A)</sup></td><td>1.0162e+03 s⁻¹</td></tr></table>

<h2>5. 第四步：把两个来源的信息合并</h2>
<p>正确的同序关系是：</p>
<table><tr><th>来源</th><th>第一条有效记录提供的内容</th></tr>
<tr><td>.lines[1] → states 映射</td><td>端点 state ID：1 和 4413</td></tr>
<tr><td>.agafgf[1]</td><td>A = 1.0162×10³ s⁻¹；ν̃ = 217420.587 cm⁻¹</td></tr></table>
<p>合并前还要判断哪个是上态、哪个是下态。能量更高的是上态：</p>
<pre>state 4413：E = 217420.587 cm⁻¹ → 上态
state    1：E =      0.000 cm⁻¹ → 下态</pre>
<p>ExoMol transition 文件的列顺序为：</p>
<table><tr><th>第1列</th><th>第2列</th><th>第3列</th><th>第4列（本项目扩展）</th></tr>
<tr><td>上态 ID</td><td>下态 ID</td><td>Einstein A (s⁻¹)</td><td>波数 (cm⁻¹)</td></tr></table>
<p>所以最终输出为：</p>
<pre>        4413            1  1.0162e+03    2.174206e+05
          ↑            ↑       ↑              ↑
       上态ID        下态ID   A系数           波数</pre>

<h2>6. 用能量差做独立物理校验</h2>
<p>在写出 transition 之前，可以用 states 的能量差检查 agafgf 波数：</p>
<div class="equation">|E(4413) − E(1)| = |217420.587 − 0| = 217420.587 cm⁻¹</div>
<p>它与第一条 <code>.agafgf</code> 的波数完全一致：</p>
<table><tr><th>量</th><th>结果</th></tr>
<tr><td>states 能量差</td><td>217420.587 cm⁻¹</td></tr>
<tr><td>agafgf 记录波数</td><td>217420.587 cm⁻¹</td></tr>
<tr><td>残差</td><td>0.000 cm⁻¹</td></tr></table>
<div class="ok"><b>校验通过：</b>这说明 `.lines` 第一条端点与 `.agafgf` 第一条物理参数属于同一次跃迁。</div>

<h2>7. 一行错位问题发生在哪里？</h2>
<p>如果错误地跳过 <code>.agafgf</code> 的第一条有效数值，就会变成：</p>
<pre>.lines[1]   → 仍然得到 state 4413 和 state 1
.agafgf[2] → A = 5.3951×10²，wn = 217326.5

错误组合：4413  1  5.3951×10²  217326.5</pre>
<p>此时 state ID 没有错，但 A 和波数属于下一条跃迁：</p>
<div class="equation">|217326.5 − 217420.587| = 94.087 cm⁻¹ ≠ 0</div>
<div class="warning"><b>所以问题本质是：</b>“能级端点”和“跃迁物理参数”来自不同的记录。不是 4413 或 1 编错了，也不是最终展示 label 的转换造成了错误。</div>

<h2>8. 旧式按行拼接与当前稳健做法</h2>
<table><tr><th>做法</th><th>机制</th><th>风险/优点</th></tr>
<tr><td>旧式按位置拼接</td><td><code>.lines[n]</code> 的端点与 <code>.agafgf[n]</code> 的 A、波数直接横向拼接</td><td>任一文件多跳或少跳一行，后续记录就可能整体错位</td></tr>
<tr><td>当前稳健解析</td><td>优先从扩展 `.agafgf` 同一条记录中同时读取端点、A 和波数，再映射 states</td><td>端点和物理参数来自同一物理行，不依赖两个文件的位置同步</td></tr>
<tr><td>必要校验</td><td>对每条记录比较 agafgf 波数与 states 能量差</td><td>即使格式或来源改变，也能发现物理不一致</td></tr></table>
<div class="note">如果处理的是只含波数与 A、不含端点的简化版 <code>.agafgf</code>，仍需与 <code>.lines</code> 同序配对，但必须先分别识别有效数值记录，并对配对后的每条跃迁执行能量差校验。</div>

<h2>9. 对应的伪代码</h2>
<pre># 1. 建立原始能级到 state ID 的索引
state_index[(abs(E), J, raw_label)] = state_id

# 2. 逐条读取 .lines 的跃迁端点
for n, line_record in enumerate(valid_lines_records):
    id1 = state_index[(abs(E1), J1, label1)]
    id2 = state_index[(abs(E2), J2, label2)]

    # 3. 读取同序 .agafgf 有效记录
    wn, log_A = valid_agafgf_records[n]
    A = 10 ** log_A

    # 4. 根据能量确定上态和下态
    upper, lower = order_by_energy(id1, id2)

    # 5. 物理一致性检查
    residual = abs(wn - abs(E[upper] - E[lower]))
    assert residual < tolerance

    # 6. 写出 .trans
    write(upper, lower, A, wn)</pre>

<h2>10. 本例的完整追踪</h2>
<pre>.lines[1]
├─ (0.000, 1.5, "(3F)4s a4F")
│      └─→ state ID 1，E=0.000，下态
└─ (-217420.587, 2.5, "m(1S)2F344")
       └─ abs(E)，按原始label查中间states
              └─→ state ID 4413，E=217420.587，上态
                   └─ 最终展示为 s211f(1S) / 2F

.agafgf[1]
├─ wn = 217420.587 cm⁻¹
└─ log(A) = 3.007... → A = 1.0162×10³ s⁻¹

合并并校验
└─→ 4413  1  1.0162e+03  2.174206e+05</pre>

<p class="small">本文说明的是 Ti II 数据转换中的字段关系和可复现映射逻辑。原始 label 用于内部匹配；最终 `.states` 中的 Configuration/Term 是对该 label 的解析展示。</p>
</body></html>"""


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(HTML, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
