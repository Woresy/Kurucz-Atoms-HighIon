"""Generate the Ti II alignment issue meeting report as HTML and DOCX."""

from __future__ import annotations

import html
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape


OUT = Path("reports/meeting-2026-07-17")
TITLE = "Ti II Kurucz 跃迁数据行对齐问题核查报告"
SUBTITLE = "用于 2026 年 7 月 17 日线下讨论"


def h(text: str) -> str:
    return html.escape(text)


def make_html() -> str:
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>{TITLE}</title>
<style>
@page {{ size: A4; margin: 18mm 17mm 18mm 20mm; @bottom-center {{ content: counter(page); font-size: 9pt; color: #666; }} }}
body {{ font-family: 'Noto Sans CJK SC', sans-serif; color:#172033; font-size:10.5pt; line-height:1.65; }}
h1 {{ font-size:22pt; color:#133b5c; text-align:center; margin-top:45mm; }}
.subtitle {{ text-align:center; font-size:13pt; color:#506273; margin:10mm 0 55mm; }}
.meta {{ border-top:1px solid #9bacbb; padding-top:6mm; color:#506273; }}
h2 {{ font-size:15pt; color:#145374; border-bottom:2px solid #4aa3a2; padding-bottom:2mm; margin-top:9mm; }}
h3 {{ font-size:12pt; color:#145374; margin-top:6mm; }}
.summary {{ background:#eef7f7; border-left:5px solid #2b7a78; padding:4mm 5mm; }}
.fact {{ background:#f3f6fa; padding:3mm 5mm; border-radius:2mm; }}
.warning {{ background:#fff5e6; border-left:5px solid #d88920; padding:4mm 5mm; }}
.key {{ font-size:12pt; background:#e8f4ff; border:2px solid #3b82a0; padding:5mm; border-radius:2mm; }}
.flow {{ background:#f8fafc; border:1px solid #9bacbb; padding:4mm 5mm; line-height:1.9; }}
table {{ width:100%; border-collapse:collapse; margin:4mm 0 6mm; font-size:9.5pt; }}
th {{ background:#145374; color:white; }} th,td {{ border:1px solid #9bacbb; padding:2.2mm; vertical-align:top; }}
tr:nth-child(even) td {{ background:#f7fafc; }}
code,pre {{ font-family:'DejaVu Sans Mono',monospace; }}
pre {{ background:#18222d; color:#eef5f7; padding:4mm; border-radius:2mm; white-space:pre-wrap; font-size:8.8pt; }}
.equation {{ text-align:center; font-family:'DejaVu Serif',serif; font-size:12pt; margin:4mm; }}
a {{ color:#075985; word-break:break-all; }}
.pagebreak {{ page-break-before:always; }}
ul,ol {{ padding-left:7mm; }}
.small {{ font-size:8.8pt; color:#5b6875; }}
</style></head><body>
<h1>{TITLE}</h1><div class="subtitle">{SUBTITLE}</div>
<div class="meta"><b>项目：</b>Kurucz 原子数据到 ExoMol 格式的复现与核验<br>
<b>讨论对象：</b>Ti II（Kurucz 代码 2201）<br>
<b>报告性质：</b>可复现的技术核查；区分已验证事实与原因推测</div>

<div class="pagebreak"></div><h2>1. 先说结论：到底发现了什么？</h2>
<div class="key"><b>一句话结论：</b>Ti II 的两份 Kurucz 源文件本应逐行配对，但 ExoAtom 发布版看起来把其中一份文件向上错开了一行：第一条有效记录丢失，后面的 A 系数和波数可能配到了前一条能级跃迁上。</div>
<p><b>最容易理解的类比：</b>把两份文件想成 Excel 中并排的两列。左列写“哪两个能级发生跃迁”，右列写“这次跃迁的波数和 A 系数”。正确做法是第 1 行配第 1 行；当前现象却像是左列第 1 行配了右列第 2 行。</p>
<div class="flow"><b>正确配对：</b>lines 第1条 ←→ agafgf 第1条　✓ 波数吻合<br>
<b>疑似发布版：</b>lines 第1条 ←→ agafgf 第2条　✗ 波数不吻合<br>
<b>直接结果：</b>897,985 条 → 少 1 条 → 发布版剩 897,984 条</div>
<h3>只需记住三个数字</h3>
<ul><li><b>897,985</b>：两份 Kurucz 源文件各自的有效记录数。</li>
<li><b>897,984</b>：ExoAtom 发布版的记录数，少 1 条。</li>
<li><b>0 与 88.586 cm⁻¹</b>：正确配对的首条残差为 0；错开一行后，前 20 条的平均绝对残差为 88.586 cm⁻¹。</li></ul>
<p>因此，最稳妥的表述是“<b>发布版呈现 off-by-one（一行错位）特征</b>”。具体是哪一步代码造成的，目前仍是推测，需要 ExoAtom 作者结合中间文件确认。</p>
<table><tr><th>数据集</th><th>有效跃迁数</th><th>第一条有效数据处理</th><th>波数一致性</th></tr>
<tr><td>Kurucz 原始文件</td><td>897,985</td><td>应保留</td><td>同一行严格一致</td></tr>
<tr><td>本地修正版</td><td>897,985</td><td>保留</td><td>通过能级差检查</td></tr>
<tr><td>ExoAtom 发布文件</td><td>897,984</td><td>未出现</td><td>首个能级对使用下一行数值</td></tr></table>

<h2>2. 为什么能判断发生了错位？</h2>
<p>判断依据很简单：同一次跃迁的波数，应当等于两个能级能量之差。也就是说，<code>.lines</code> 自己就能给出一个“预期波数”，再拿它与 <code>.agafgf</code> 记录的波数比较。两者一致，说明配对正确；差得很大，说明配错了行。</p>
<div class="equation">检查量：残差 = |agafgf 记录波数 − lines 上下能级差|</div>
<table><tr><th>配对方法</th><th>前 20 条平均绝对残差</th><th>如何解释</th></tr>
<tr><td>第 n 条 lines + 第 n 条 agafgf</td><td>约 2.9×10⁻¹² cm⁻¹</td><td>等于 0（仅有浮点误差），配对正确</td></tr>
<tr><td>第 n 条 lines + 第 n+1 条 agafgf</td><td>88.586 cm⁻¹</td><td>明显不一致，说明错开一行</td></tr></table>

<h2>3. 用第一条数据走一遍</h2>
<p><code>gf2201.agafgf</code> 的第 1 个物理行只是列名；下一行才是第 1 条真正的数据。处理时应当“跳过列名，但保留第 1 条数值”。</p>
<pre>lines 第1条给出的能级：       0.000 → −217420.587 cm⁻¹
由能级差算出的预期波数：     217420.587 cm⁻¹
agafgf 第1条记录的波数：      217420.587 cm⁻¹
比较结果：                    完全一致，残差 = 0</pre>
<p>这证明 agafgf 的第 1 条数值不是应被删除的“多余行”，而是一条有效跃迁记录。</p>

<h2>4. ExoAtom 发布版具体哪里不一样？</h2>
<p><b>先解释 1 和 4413：</b>它们是 ExoMol <code>.states</code> 文件为每个能级分配的唯一编号（state ID），不是原始文件的行号，也不是能量值。生成 states 文件时，能级按能量排序并从 1 开始编号。</p>
<table><tr><th>state ID</th><th>对应能量</th><th>在本次跃迁中的角色</th><th>与上文的关系</th></tr>
<tr><td>1</td><td>0.000 cm⁻¹</td><td>下能级</td><td>对应上文的 E₁ = 0.000</td></tr>
<tr><td>4413</td><td>217420.587 cm⁻¹</td><td>上能级</td><td>对应上文取绝对值后的 E₂ = −217420.587</td></tr></table>
<p>因此，<code>4413  1</code> 的意思是“从编号 4413 的上能级跃迁到编号 1 的下能级”。两者的能量差正好是：</p>
<div class="equation">E(4413) − E(1) = 217420.587 − 0 = 217420.587 cm⁻¹</div>
<p>本地结果与 ExoAtom 比较的是<b>同一个能级对（4413 → 1）</b>，差别不在 state ID，而在分配给这个能级对的 A 系数和波数：</p>
<table><tr><th>版本</th><th>上态 ID</th><th>下态 ID</th><th>A 系数 (s⁻¹)</th><th>波数 (cm⁻¹)</th></tr>
<tr><td>正确同序配对</td><td>4413</td><td>1</td><td>1.0162×10³</td><td>217420.6</td></tr>
<tr><td>ExoAtom 发布版</td><td>4413</td><td>1</td><td>5.3951×10²</td><td>217326.5</td></tr></table>
<p>ExoAtom 保留了第 1 条 <code>.lines</code> 所指向的 state ID（4413 和 1），却为它们配上了下一条 <code>.agafgf</code> 的 <b>A=5.3951×10²、波数=217326.5</b>。这才是“一行错位”的指纹。</p>

<h2>5. 数据来源</h2>
<h3>5.1 Kurucz 官方原始文件</h3>
<ul><li>目录：<a href="http://kurucz.harvard.edu/atoms/2201/">http://kurucz.harvard.edu/atoms/2201/</a></li>
<li>能级端点：<a href="http://kurucz.harvard.edu/atoms/2201/gf2201.lines">gf2201.lines</a></li>
<li>波数与 A：<a href="http://kurucz.harvard.edu/atoms/2201/gf2201.agafgf">gf2201.agafgf</a></li></ul>
<p>依据 ExoAtom 论文的 Kurucz 数据说明，<code>.lines</code> 提供跃迁上下能级，<code>.agafgf</code> 提供波数和 log(A)。两者需要正确保持记录对应关系。</p>
<h3>5.2 ExoAtom 发布文件</h3>
<ul><li>States：<a href="https://www.exomol.com/exoatom/db/Ti_p/Kurucz/Ti_p__Kurucz.states">Ti_p__Kurucz.states</a></li>
<li>Transitions：<a href="https://www.exomol.com/exoatom/db/Ti_p/Kurucz/Ti_p__Kurucz.trans">Ti_p__Kurucz.trans</a></li>
<li>论文：<a href="https://academic.oup.com/rasti/article/doi/10.1093/rasti/rzaf065/8404150">ExoAtom: a database of atomic spectra in ExoMol format</a></li></ul>

<h2>6. 完整数值核验</h2>
<p><code>gf2201.agafgf</code> 的物理第一行是文字表头，第二行才是第一条有效数值记录。正确逻辑是“跳过文字表头，但保留第一条有效数据”。</p>
<pre>表头：wl(nm)  wn(cm-1)  log gf  log f  log fe  log A  log gA
第一条有效 agafgf：-45.9938  217420.587 ... log A = 3.007 ...
第一条有效 lines： E1 = 0.000, J1 = 1.5; E2 = -217420.587, J2 = 2.5</pre>
<div class="equation">ν̃<sub>expected</sub> = ||E₂| − |E₁|| = |217420.587 − 0| = 217420.587 cm⁻¹</div>
<table><tr><th>量</th><th>数值</th><th>来源</th></tr>
<tr><td>预期波数</td><td>217420.587 cm⁻¹</td><td><code>.lines</code> 能级差</td></tr>
<tr><td>记录波数</td><td>217420.587 cm⁻¹</td><td>第一条有效 <code>.agafgf</code></td></tr>
<tr><td>Einstein A</td><td>约 1.0162×10³ s⁻¹</td><td>第一条有效 <code>.agafgf</code></td></tr>
<tr><td>残差</td><td>0.000 cm⁻¹</td><td>二者之差</td></tr></table>

<div class="pagebreak"></div><h3>6.1 前 20 条偏移量对照实验</h3>
<p>对前20条有效记录分别测试两种配对：</p>
<table><tr><th>方案</th><th>配对方式</th><th>平均绝对残差</th><th>最大绝对残差</th></tr>
<tr><td>A：不偏移</td><td>lines[n] + agafgf[n]</td><td>约 2.9×10⁻¹² cm⁻¹</td><td>约 2.9×10⁻¹¹ cm⁻¹</td></tr>
<tr><td>B：跳过首条有效 A/波数</td><td>lines[n] + agafgf[n+1]</td><td>88.586 cm⁻¹</td><td>394.361 cm⁻¹</td></tr></table>
<p>方案 A 的误差仅来自浮点表示；方案 B 出现几十至数百 cm⁻¹ 的残差。因此，第一条有效 <code>.agafgf</code> 不应被跳过。</p>

<h3>6.2 与 ExoAtom 文件的直接对应</h3>
<p>本地按物理一致性配对后的首条输出为：</p>
<pre>4413  1  1.0162e+03  2.174206e+05</pre>
<p>ExoAtom 中没有这条组合。ExoAtom 对相同能级对给出：</p>
<pre>4413  1  5.395100E+02  2.173265E+05</pre>
<p>其中 A=5.3951×10²、ν̃=217326.5 恰好对应下一条有效 <code>.agafgf</code>。这一现象符合“端点仍从第 n 条 <code>.lines</code> 获取，而 A/波数从第 n+1 条获取”的模式。</p>

<h2>7. 哪些是事实，哪些只是原因推测？</h2>
<div class="fact"><b>已经验证：</b><ul><li>第一条原始 agafgf 物理行是文字表头。</li><li>表头后的第一条数值记录有效，波数与能级差完全一致。</li><li>本地保留 897,985 条；ExoAtom 为 897,984 条。</li><li>ExoAtom 首个相关能级对使用下一条 agafgf 的数值。</li></ul></div>
<div class="warning"><b>合理推测，尚需作者确认：</b>旧流程可能先通过 <code>read_fwf/to_csv</code> 或人工步骤移除了原始表头，随后在读取已经无表头的中间 CSV 时再次执行 <code>skiprows=1</code>，从而删掉第一条有效数据。生成网站文件所用的中间 CSV 未公开，因此不能把具体发生步骤写成确定事实。</div>

<h2>8. 这会产生什么影响？</h2>
<ul><li>不仅是总数少一条；若偏移系统性存在，A 系数和波数可能被分配给错误的能级对。</li>
<li>Einstein A 影响谱线强度、寿命与辐射转移计算。</li><li>波数与能级差不一致会影响谱线归属与可追溯性。</li>
<li>网站文件已按波数排序，肉眼查看最终文件不容易发现；必须回溯能级对才能识别。</li></ul>

<div class="pagebreak"></div><h2>9. 建议会议当场决定的事项</h2>
<ol><li>是否请组内另一位成员独立复现前20条残差测试？</li>
<li>本地 Ti II 是否保留物理一致的 897,985 条版本？建议：保留。</li>
<li>是否联系 ExoAtom 作者，请其检查 Ti II 生成时的中间 <code>AGAFGF.csv</code> 和 <code>skiprows</code> 设置？</li>
<li>是否对“Ti 及之后”的其他物种批量执行同样的波数一致性检查？</li>
<li>论文/毕业报告中如何记录本地修订与网站发布版之间的差异？</li></ol>

<h2>10. 推荐对外表述</h2>
<blockquote>在复现 ExoAtom 的 Ti II Kurucz 数据时，我们观察到发布版 <code>.trans</code> 相比 Kurucz 两个源文件的有效记录少一条。基于上下能级差与 <code>.agafgf</code> 波数的独立校验，同序配对在浮点精度内一致，而偏移一行会产生显著残差。发布版首个相关能级对使用了下一条 <code>.agafgf</code> 的 A 和波数。该现象提示中间预处理可能存在 off-by-one 对齐问题，建议结合原始中间文件进一步确认。</blockquote>

<h2>附录 A：最小复现方法</h2>
<pre>1. 下载 gf2201.lines 与 gf2201.agafgf。
2. 从 .lines 固定列读取 E1、J1、label1、E2、J2、label2。
3. 从 .agafgf 固定列读取 wavenumber 与 log(A)，自动忽略文字表头。
4. 计算 residual = abs(wavenumber - abs(abs(E2)-abs(E1)))。
5. 分别测试 agafgf offset=0 和 offset=1。
6. 将结果与 Ti_p__Kurucz.trans 中相同能级对比较。</pre>
<p class="small">报告生成日期：2026-07-16。数据核验针对报告生成时上述官方链接提供的文件。</p>
</body></html>"""


def make_html_en() -> str:
    """Return a clear English version for meetings and external review."""
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Investigation of a Row-Alignment Issue in Ti II Kurucz Transition Data</title>
<style>
@page { size: A4; margin: 18mm 17mm 18mm 20mm; @bottom-center { content: counter(page); font-size: 9pt; color: #666; } }
body { font-family: Arial, 'Noto Sans', sans-serif; color:#172033; font-size:10.5pt; line-height:1.65; }
h1 { font-size:22pt; color:#133b5c; text-align:center; margin-top:45mm; }
.subtitle { text-align:center; font-size:13pt; color:#506273; margin:10mm 0 55mm; }
.meta { border-top:1px solid #9bacbb; padding-top:6mm; color:#506273; }
h2 { font-size:15pt; color:#145374; border-bottom:2px solid #4aa3a2; padding-bottom:2mm; margin-top:9mm; }
h3 { font-size:12pt; color:#145374; margin-top:6mm; }
.key { font-size:12pt; background:#e8f4ff; border:2px solid #3b82a0; padding:5mm; border-radius:2mm; }
.flow { background:#f8fafc; border:1px solid #9bacbb; padding:4mm 5mm; line-height:1.9; }
.fact { background:#f3f6fa; padding:3mm 5mm; border-radius:2mm; }
.warning { background:#fff5e6; border-left:5px solid #d88920; padding:4mm 5mm; }
table { width:100%; border-collapse:collapse; margin:4mm 0 6mm; font-size:9.5pt; }
th { background:#145374; color:white; } th,td { border:1px solid #9bacbb; padding:2.2mm; vertical-align:top; }
tr:nth-child(even) td { background:#f7fafc; }
code,pre { font-family:'DejaVu Sans Mono',monospace; }
pre { background:#18222d; color:#eef5f7; padding:4mm; border-radius:2mm; white-space:pre-wrap; font-size:8.8pt; }
.equation { text-align:center; font-family:'DejaVu Serif',serif; font-size:12pt; margin:4mm; }
a { color:#075985; word-break:break-all; }
.pagebreak { page-break-before:always; }
ul,ol { padding-left:7mm; }
.small { font-size:8.8pt; color:#5b6875; }
</style></head><body>
<h1>Investigation of a Row-Alignment Issue in Ti II Kurucz Transition Data</h1>
<div class="subtitle">Prepared for the in-person discussion on 17 July 2026</div>
<div class="meta"><b>Project:</b> Reproduction and validation of Kurucz atomic data in ExoMol format<br>
<b>Species:</b> Ti II (Kurucz code 2201)<br>
<b>Scope:</b> A reproducible technical investigation that separates verified evidence from hypotheses</div>

<div class="pagebreak"></div><h2>1. The finding in one sentence</h2>
<div class="key"><b>Bottom line:</b> The two Ti II Kurucz source files should be paired row by row. The ExoAtom release, however, appears to shift one file upward by one row: the first valid record is absent, and subsequent Einstein A coefficients and wavenumbers may have been assigned to the preceding energy-level pair.</div>
<p><b>A simple analogy:</b> Imagine two adjacent columns in a spreadsheet. The left column identifies the two energy levels involved in a transition; the right column gives the transition's wavenumber and A coefficient. The correct operation matches row 1 with row 1. The observed release behaves as though row 1 on the left were matched with row 2 on the right.</p>
<div class="flow"><b>Correct:</b> lines record 1 ↔ agafgf record 1 &nbsp; ✓ wavenumbers agree<br>
<b>Observed pattern:</b> lines record 1 ↔ agafgf record 2 &nbsp; ✗ wavenumbers disagree<br>
<b>Record count:</b> 897,985 → one record lost → 897,984 in the release</div>
<h3>Three numbers to remember</h3>
<ul><li><b>897,985:</b> the number of valid records in each Kurucz source file.</li>
<li><b>897,984:</b> the number of records in the ExoAtom transition file—one fewer.</li>
<li><b>0 versus 88.586 cm⁻¹:</b> the first correctly matched record has zero residual; shifting by one row gives a mean absolute residual of 88.586 cm⁻¹ over the first 20 records.</li></ul>
<p>The safest conclusion is therefore that the release <b>shows the signature of an off-by-one row-alignment issue</b>. The exact processing step responsible has not been verified and requires confirmation from the ExoAtom authors.</p>

<h2>2. How can the shift be detected?</h2>
<p>For a given transition, the wavenumber must equal the energy difference between its upper and lower levels. The <code>.lines</code> file therefore provides an independently calculated expected wavenumber, which can be compared with the wavenumber stored in <code>.agafgf</code>. Agreement indicates a correct match; a large difference indicates that the rows have been paired incorrectly.</p>
<div class="equation">residual = |agafgf wavenumber − energy-level difference from lines|</div>
<table><tr><th>Pairing</th><th>Mean absolute residual, first 20 records</th><th>Interpretation</th></tr>
<tr><td>lines[n] + agafgf[n]</td><td>≈ 2.9×10⁻¹² cm⁻¹</td><td>Effectively zero; correct pairing</td></tr>
<tr><td>lines[n] + agafgf[n+1]</td><td>88.586 cm⁻¹</td><td>Clear physical inconsistency; one-row shift</td></tr></table>

<h2>3. Walking through the first record</h2>
<p>The physical first line of <code>gf2201.agafgf</code> is a text header. The following line is the first valid numerical record. The correct rule is therefore: skip the text header, but retain the first numerical record.</p>
<pre>Energy levels in lines record 1:       0.000 → −217420.587 cm⁻¹
Expected wavenumber from the levels:    217420.587 cm⁻¹
Wavenumber in agafgf record 1:          217420.587 cm⁻¹
Result:                                 exact agreement; residual = 0</pre>
<p>This confirms that the first numerical <code>.agafgf</code> record is valid transition data, not an extra line that should be removed.</p>

<h2>4. What is different in the ExoAtom release?</h2>
<p><b>What do 1 and 4413 mean?</b> They are unique state IDs assigned to energy levels in the ExoMol <code>.states</code> file. They are neither source-file row numbers nor energy values. The states are sorted by energy and numbered from 1.</p>
<table><tr><th>State ID</th><th>Energy</th><th>Role</th><th>Connection to Section 3</th></tr>
<tr><td>1</td><td>0.000 cm⁻¹</td><td>Lower state</td><td>The level E₁ = 0.000</td></tr>
<tr><td>4413</td><td>217420.587 cm⁻¹</td><td>Upper state</td><td>The absolute-energy form of E₂ = −217420.587</td></tr></table>
<p>Thus, <code>4413  1</code> means a transition from upper state 4413 to lower state 1, whose energy difference is:</p>
<div class="equation">E(4413) − E(1) = 217420.587 − 0 = 217420.587 cm⁻¹</div>
<p>The local result and the ExoAtom release refer to the <b>same state pair (4413 → 1)</b>. What differs is the A coefficient and wavenumber assigned to that pair:</p>
<table><tr><th>Version</th><th>Upper-state ID</th><th>Lower-state ID</th><th>A coefficient (s⁻¹)</th><th>Wavenumber (cm⁻¹)</th></tr>
<tr><td>Correct row-wise pairing</td><td>4413</td><td>1</td><td>1.0162×10³</td><td>217420.6</td></tr>
<tr><td>ExoAtom release</td><td>4413</td><td>1</td><td>5.3951×10²</td><td>217326.5</td></tr></table>
<p>ExoAtom retains the state IDs from the first <code>.lines</code> record, but assigns them <b>A = 5.3951×10² s⁻¹ and ν̃ = 217326.5 cm⁻¹</b> from the next valid <code>.agafgf</code> record. This—not the ordering of the state IDs—is the direct signature of a one-row shift.</p>

<h2>5. Data sources</h2>
<h3>5.1 Official Kurucz files</h3>
<ul><li>Directory: <a href="http://kurucz.harvard.edu/atoms/2201/">http://kurucz.harvard.edu/atoms/2201/</a></li>
<li>Transition endpoints: <a href="http://kurucz.harvard.edu/atoms/2201/gf2201.lines">gf2201.lines</a></li>
<li>Wavenumbers and A values: <a href="http://kurucz.harvard.edu/atoms/2201/gf2201.agafgf">gf2201.agafgf</a></li></ul>
<h3>5.2 ExoAtom release</h3>
<ul><li>States: <a href="https://www.exomol.com/exoatom/db/Ti_p/Kurucz/Ti_p__Kurucz.states">Ti_p__Kurucz.states</a></li>
<li>Transitions: <a href="https://www.exomol.com/exoatom/db/Ti_p/Kurucz/Ti_p__Kurucz.trans">Ti_p__Kurucz.trans</a></li>
<li>Paper: <a href="https://academic.oup.com/rasti/article/doi/10.1093/rasti/rzaf065/8404150">ExoAtom: a database of atomic spectra in ExoMol format</a></li></ul>

<div class="pagebreak"></div><h2>6. Full numerical checks</h2>
<h3>6.1 First-record quantities</h3>
<table><tr><th>Quantity</th><th>Value</th><th>Source</th></tr>
<tr><td>Expected wavenumber</td><td>217420.587 cm⁻¹</td><td>Energy-level difference in <code>.lines</code></td></tr>
<tr><td>Recorded wavenumber</td><td>217420.587 cm⁻¹</td><td>First valid <code>.agafgf</code> record</td></tr>
<tr><td>Einstein A</td><td>≈ 1.0162×10³ s⁻¹</td><td>First valid <code>.agafgf</code> record</td></tr>
<tr><td>Residual</td><td>0.000 cm⁻¹</td><td>Difference between the two wavenumbers</td></tr></table>
<h3>6.2 Offset test over the first 20 records</h3>
<table><tr><th>Case</th><th>Pairing</th><th>Mean absolute residual</th><th>Maximum absolute residual</th></tr>
<tr><td>No offset</td><td>lines[n] + agafgf[n]</td><td>≈ 2.9×10⁻¹² cm⁻¹</td><td>≈ 2.9×10⁻¹¹ cm⁻¹</td></tr>
<tr><td>One-row offset</td><td>lines[n] + agafgf[n+1]</td><td>88.586 cm⁻¹</td><td>394.361 cm⁻¹</td></tr></table>

<h2>7. Verified facts versus inferred cause</h2>
<div class="fact"><b>Verified:</b><ul><li>The first physical line of the original <code>.agafgf</code> file is a text header.</li><li>The first numerical record after that header is valid and its wavenumber exactly matches the energy-level difference.</li><li>The local physically consistent result has 897,985 records; the ExoAtom release has 897,984.</li><li>For the first relevant level pair, the ExoAtom release uses values from the next valid <code>.agafgf</code> record.</li></ul></div>
<div class="warning"><b>Plausible but unverified explanation:</b> An earlier processing step may have removed the original header when creating an intermediate CSV, after which <code>skiprows=1</code> may have been applied again to that already headerless CSV. This would delete the first valid record. The intermediate file used to generate the website release is not public, so this mechanism must not be presented as established fact.</div>

<h2>8. Why does this matter scientifically?</h2>
<ul><li>The issue is not limited to one missing line. If the shift is systematic, A coefficients and wavenumbers may be assigned to the wrong energy-level pairs.</li>
<li>Einstein A coefficients affect line intensities, radiative lifetimes, and radiative-transfer calculations.</li>
<li>A mismatch between wavenumber and energy difference affects line identification and traceability.</li>
<li>Because the published transition file is sorted by wavenumber, visual inspection alone is unlikely to expose the problem; the level pairs must be reconstructed.</li></ul>

<div class="pagebreak"></div><h2>9. Decisions proposed for the meeting</h2>
<ol><li>Should another group member independently reproduce the residual test for the first 20 records?</li>
<li>Should the local Ti II dataset retain the physically consistent 897,985-record version? <b>Recommendation: yes.</b></li>
<li>Should the ExoAtom authors be contacted and asked to inspect the intermediate <code>AGAFGF.csv</code> file and any <code>skiprows</code> setting used for Ti II?</li>
<li>Should the same wavenumber-consistency test be run in bulk for Ti and later species?</li>
<li>How should the difference between the local correction and the website release be documented in the thesis or final report?</li></ol>

<h2>10. Suggested wording for external communication</h2>
<blockquote>While reproducing the Ti II Kurucz data in ExoAtom format, we found that the released transition file contains one fewer valid record than the two corresponding Kurucz source files. An independent comparison of energy-level differences with the wavenumbers in <code>.agafgf</code> shows that row-wise pairing agrees to floating-point precision, whereas a one-row offset produces substantial residuals. For the first relevant level pair, the release uses the A coefficient and wavenumber from the next valid <code>.agafgf</code> record. This pattern suggests a possible off-by-one alignment issue during intermediate preprocessing. Confirmation against the original intermediate files is recommended.</blockquote>

<h2>Appendix A: Minimal reproduction procedure</h2>
<pre>1. Download gf2201.lines and gf2201.agafgf.
2. Read E1, J1, label1, E2, J2, and label2 from the fixed-width lines file.
3. Read wavenumber and log(A) from agafgf, ignoring only the text header.
4. Calculate residual = abs(wavenumber - abs(abs(E2) - abs(E1))).
5. Compare agafgf offsets 0 and 1.
6. Compare the result with the same level pair in Ti_p__Kurucz.trans.</pre>
<p class="small">Report generated on 16 July 2026. The checks refer to the files available from the official links above at the time of analysis.</p>
</body></html>"""


def run(text: str, bold: bool = False, size: int | None = None) -> str:
    props = []
    if bold:
        props.append("<w:b/>")
    if size:
        props.append(f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>')
    rpr = f"<w:rPr>{''.join(props)}</w:rPr>" if props else ""
    return f'<w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r>'


def para(text: str = "", style: str | None = None, bold: bool = False, align: str | None = None) -> str:
    ppr = []
    if style:
        ppr.append(f'<w:pStyle w:val="{style}"/>')
    if align:
        ppr.append(f'<w:jc w:val="{align}"/>')
    return f'<w:p><w:pPr>{"".join(ppr)}</w:pPr>{run(text, bold)}</w:p>'


def table(rows: list[list[str]]) -> str:
    parts = ['<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/></w:tblPr>']
    for ri, row in enumerate(rows):
        parts.append("<w:tr>")
        for cell in row:
            shade = '<w:shd w:fill="DCE6F1"/>' if ri == 0 else ""
            parts.append(f'<w:tc><w:tcPr>{shade}</w:tcPr>{para(cell, bold=(ri == 0))}</w:tc>')
        parts.append("</w:tr>")
    parts.append("</w:tbl>")
    return "".join(parts)


def make_docx(path: Path) -> None:
    body = [
        para(TITLE, "Title", align="center"), para(SUBTITLE, "Subtitle", align="center"),
        para("项目：Kurucz 原子数据到 ExoMol 格式的复现与核验"),
        para("讨论对象：Ti II（Kurucz 代码 2201）"),
        para("报告性质：可复现的技术核查；区分已验证事实与原因推测"),
        para("1. 先说结论：到底发现了什么？", "Heading1"),
        para("一句话结论：Ti II 的两份 Kurucz 源文件本应逐行配对，但 ExoAtom 发布版看起来把其中一份文件向上错开了一行：第一条有效记录丢失，后面的 A 系数和波数可能配到了前一条能级跃迁上。", bold=True),
        para("最容易理解的类比：把两份文件想成 Excel 中并排的两列。左列写“哪两个能级发生跃迁”，右列写“这次跃迁的波数和 A 系数”。正确做法是第1行配第1行；当前现象却像是左列第1行配了右列第2行。"),
        para("正确配对：lines 第1条 ↔ agafgf 第1条，波数吻合。"),
        para("疑似发布版：lines 第1条 ↔ agafgf 第2条，波数不吻合。"),
        para("三个关键数字：Kurucz 为897,985条；ExoAtom为897,984条；错开一行后前20条的平均绝对残差为88.586 cm⁻¹。"),
        table([["数据集","有效跃迁数","第一条有效数据","波数一致性"],["Kurucz 原始文件","897,985","应保留","同一行严格一致"],["本地修正版","897,985","保留","通过能级差检查"],["ExoAtom 发布文件","897,984","未出现","首个能级对使用下一行数值"]]),
        para("2. 为什么能判断发生了错位？", "Heading1"),
        para("同一次跃迁的波数应当等于两个能级的能量差。因此可以用 lines 算出预期波数，再与 agafgf 记录的波数比较。两者一致说明配对正确，差得很大说明配错了行。"),
        para("检查量：残差 = |agafgf记录波数 − lines上下能级差|"),
        table([["配对方法","前20条平均绝对残差","解释"],["lines[n]+agafgf[n]","约2.9×10⁻¹² cm⁻¹","等于0，仅有浮点误差"],["lines[n]+agafgf[n+1]","88.586 cm⁻¹","明显不一致，错开一行"]]),
        para("3. 用第一条数据走一遍", "Heading1"),
        para("agafgf 的第1个物理行只是列名；下一行才是第1条真正的数据。处理时应跳过列名，但保留第1条数值。"),
        para("lines 给出的能级为0.000→−217420.587 cm⁻¹，因此预期波数为217420.587 cm⁻¹；agafgf 第1条记录也正好是217420.587 cm⁻¹，残差为0。"),
        para("4. ExoAtom 发布版具体哪里不一样？", "Heading1"),
        para("1 和 4413 是 ExoMol .states 文件为能级分配的唯一编号（state ID），不是原始文件行号或能量值。能级按能量排序并从1开始编号。"),
        table([["state ID","对应能量","角色","与上文的关系"],["1","0.000 cm⁻¹","下能级","E₁=0.000"],["4413","217420.587 cm⁻¹","上能级","取绝对值后的E₂=−217420.587"]]),
        para("因此 4413  1 表示从编号4413的上能级跃迁到编号1的下能级；能量差为217420.587−0=217420.587 cm⁻¹。"),
        table([["版本","上态ID","下态ID","A系数(s⁻¹)","波数(cm⁻¹)"],["正确同序配对","4413","1","1.0162×10³","217420.6"],["ExoAtom发布版","4413","1","5.3951×10²","217326.5"]]),
        para("两版比较的是同一个能级对。ExoAtom 保留了 lines 第1条的 state ID，却配上了 agafgf 下一条的 A=5.3951×10² 和波数217326.5，这才是“一行错位”的直接特征。"),
        para("5. 数据来源", "Heading1"),
        para("Kurucz Ti II 目录：http://kurucz.harvard.edu/atoms/2201/"),
        para("gf2201.lines：http://kurucz.harvard.edu/atoms/2201/gf2201.lines"),
        para("gf2201.agafgf：http://kurucz.harvard.edu/atoms/2201/gf2201.agafgf"),
        para("ExoAtom states：https://www.exomol.com/exoatom/db/Ti_p/Kurucz/Ti_p__Kurucz.states"),
        para("ExoAtom trans：https://www.exomol.com/exoatom/db/Ti_p/Kurucz/Ti_p__Kurucz.trans"),
        para("论文：https://academic.oup.com/rasti/article/doi/10.1093/rasti/rzaf065/8404150"),
        para("6. 完整数值核验", "Heading1"),
        para("gf2201.agafgf 的物理第一行是文字表头，第二行才是第一条有效数值记录。正确逻辑是跳过文字表头，但保留第一条有效数据。"),
        para("第一条 lines：E1=0.000 cm⁻¹，E2=-217420.587 cm⁻¹。"),
        para("第一条有效 agafgf：波数=217420.587 cm⁻¹，A≈1.0162×10³ s⁻¹。"),
        para("因此 ||E₂|−|E₁||=217420.587 cm⁻¹，与 agafgf 完全一致，残差为 0。"),
        table([["量","数值","来源"],["预期波数","217420.587 cm⁻¹","lines 能级差"],["记录波数","217420.587 cm⁻¹","第一条有效 agafgf"],["Einstein A","约 1.0162×10³ s⁻¹","第一条有效 agafgf"],["残差","0.000 cm⁻¹","二者之差"]]),
        para("6.1 前20条偏移量对照实验", "Heading1"),
        table([["方案","配对方式","平均绝对残差","最大绝对残差"],["不偏移","lines[n]+agafgf[n]","约 2.9×10⁻¹² cm⁻¹","约 2.9×10⁻¹¹ cm⁻¹"],["偏移一行","lines[n]+agafgf[n+1]","88.586 cm⁻¹","394.361 cm⁻¹"]]),
        para("结论：不偏移时误差仅来自浮点表示；偏移一行产生显著物理残差。"),
        para("6.2 与 ExoAtom 的直接对应", "Heading1"),
        para("本地物理一致的首条输出：4413  1  1.0162e+03  2.174206e+05"),
        para("ExoAtom 对相同能级对给出：4413  1  5.395100E+02  2.173265E+05"),
        para("后者的 A 和波数对应下一条有效 agafgf，符合端点取第 n 行而 A/波数取第 n+1 行的模式。"),
        para("7. 哪些是事实，哪些只是原因推测？", "Heading1"),
        para("已验证事实：第一条数值记录有效；本地897,985条，ExoAtom 897,984条；ExoAtom首个相关能级对使用下一条agafgf数值。"),
        para("合理推测：中间 CSV 已经去除表头后，又执行 skiprows=1，导致第一条有效数据被删除。由于网站生成时的中间 CSV 未公开，该具体步骤需要作者确认。"),
        para("8. 这会产生什么影响？", "Heading1"),
        para("若偏移系统性存在，Einstein A 与波数会被分配给错误能级对，可能影响谱线强度、寿命、辐射转移、谱线位置及归属。"),
        para("9. 建议会议当场决定的事项", "Heading1"),
        para("1）是否安排独立复现前20条残差测试？"), para("2）本地是否保留物理一致的897,985条版本？建议保留。"),
        para("3）是否联系ExoAtom作者检查Ti II中间AGAFGF.csv与skiprows设置？"),
        para("4）是否对Ti及之后物种批量执行波数一致性检查？"), para("5）毕业报告中如何记录本地修订与网站版差异？"),
        para("10. 推荐对外表述", "Heading1"),
        para("在复现 ExoAtom 的 Ti II Kurucz 数据时，我们观察到发布版 trans 相比 Kurucz 两个源文件的有效记录少一条。基于上下能级差与 agafgf 波数的独立校验，同序配对在浮点精度内一致，而偏移一行会产生显著残差。该现象提示中间预处理可能存在 off-by-one 对齐问题，建议结合原始中间文件进一步确认。"),
        para("附录：最小复现", "Heading1"),
        para("下载两个Kurucz文件；从lines读取端点，从agafgf读取波数和log(A)；自动忽略文字表头；计算 residual=abs(wavenumber-abs(abs(E2)-abs(E1)))；分别测试offset=0和offset=1；与ExoAtom相同能级对比较。"),
        para("报告生成日期：2026-07-16。", align="right"),
    ]
    sect = '<w:sectPr><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1020" w:right="964" w:bottom="1020" w:left="1134"/></w:sectPr>'
    document = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>' + "".join(body) + sect + '</w:body></w:document>'
    styles = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:rPr><w:rFonts w:eastAsia="Microsoft YaHei"/><w:sz w:val="21"/></w:rPr><w:pPr><w:spacing w:after="120" w:line="330" w:lineRule="auto"/></w:pPr></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:basedOn w:val="Normal"/><w:rPr><w:b/><w:color w:val="133B5C"/><w:sz w:val="40"/></w:rPr><w:pPr><w:spacing w:before="1200" w:after="360"/></w:pPr></w:style>
<w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:basedOn w:val="Normal"/><w:rPr><w:color w:val="506273"/><w:sz w:val="26"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:rPr><w:b/><w:color w:val="145374"/><w:sz w:val="30"/></w:rPr><w:pPr><w:keepNext/><w:spacing w:before="360" w:after="180"/></w:pPr></w:style>
<w:style w:type="table" w:styleId="TableGrid"><w:name w:val="Table Grid"/><w:tblPr><w:tblBorders><w:top w:val="single" w:sz="4"/><w:left w:val="single" w:sz="4"/><w:bottom w:val="single" w:sz="4"/><w:right w:val="single" w:sz="4"/><w:insideH w:val="single" w:sz="4"/><w:insideV w:val="single" w:sz="4"/></w:tblBorders></w:tblPr></w:style></w:styles>'''
    content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/></Types>'''
    rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>'''
    doc_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'''
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("word/document.xml", document)
        z.writestr("word/styles.xml", styles)
        z.writestr("word/_rels/document.xml.rels", doc_rels)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "TiII_alignment_issue_report.html").write_text(make_html(), encoding="utf-8")
    (OUT / "TiII_alignment_issue_report_EN.html").write_text(make_html_en(), encoding="utf-8")
    make_docx(OUT / "TiII_alignment_issue_report.docx")
    print(OUT)


if __name__ == "__main__":
    main()
