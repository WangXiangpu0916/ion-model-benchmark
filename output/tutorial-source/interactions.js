
'use strict';
const D=JSON.parse(document.getElementById('benchmark-data').textContent);
const SNIPPETS=JSON.parse(document.getElementById('source-snippets').textContent);
const FIGURES=JSON.parse(document.getElementById('figure-data').textContent);
const $=id=>document.getElementById(id), cm='#187996',iod='#b45623',muted='#63747c';
const mean=a=>a.reduce((s,x)=>s+x,0)/a.length;
const sd=a=>Math.sqrt(a.reduce((s,x)=>s+(x-mean(a))**2,0)/(a.length-1));
const fmt=(x,n=3)=>Number(x).toFixed(n),pct=x=>fmt(x*100,2)+'%';
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const table=(heads,rows)=>'<table><thead><tr>'+heads.map(h=>'<th>'+h+'</th>').join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+r.map(c=>'<td>'+c+'</td>').join('')+'</tr>').join('')+'</tbody></table>';
const cutoff=ion=>ion==='Mg'?.30:.35;
const label=c=>c.id.replaceAll('_',' · ');
const addOptions=(el,cases)=>el.innerHTML=cases.map(c=>`<option value="${c.id}">${label(c)}</option>`).join('');
function svgText(x,y,t,extra=''){return `<text x="${x}" y="${y}" ${extra}>${esc(t)}</text>`}
function chart(id,series,opt={}){
 const svg=$(id),W=Math.max(290,svg.clientWidth||700),H=svg.clientHeight||330,L=62,R=16,T=22,B=54,w=W-L-R,h=H-T-B;
 const points=series.flatMap(s=>s.values).filter(p=>p.every(Number.isFinite));
 let xmin=opt.xmin??Math.min(...points.map(p=>p[0])),xmax=opt.xmax??Math.max(...points.map(p=>p[0]));
 let ymin=opt.ymin??Math.min(...points.map(p=>p[1])),ymax=opt.ymax??Math.max(...points.map(p=>p[1]));
 if(ymax===ymin){ymin-=.5;ymax+=.5}else if(opt.ymin===undefined||opt.ymax===undefined){let pad=(ymax-ymin)*.12;if(opt.ymin===undefined)ymin-=pad;if(opt.ymax===undefined)ymax+=pad}
 const x=v=>L+(v-xmin)/(xmax-xmin)*w,y=v=>T+h-(v-ymin)/(ymax-ymin)*h;
 let out=`<defs><clipPath id="clip-${id}"><rect x="${L}" y="${T}" width="${w}" height="${h}"/></clipPath></defs>`;
 if(!opt.hideY)for(let i=0;i<=4;i++){let v=ymin+(ymax-ymin)*i/4,yy=y(v);out+=`<line class="gridline" x1="${L}" x2="${W-R}" y1="${yy}" y2="${yy}"/>`+svgText(L-9,yy+4,fmt(v,opt.ydigits??1),'text-anchor="end"')}
 let ticks=W<430?3:4;
 for(let i=0;i<=ticks;i++){let v=xmin+(xmax-xmin)*i/ticks,xx=x(v);out+=svgText(xx,T+h+22,fmt(v,opt.xdigits??1),`text-anchor="${i===0?'start':i===ticks?'end':'middle'}"`)}
 out+=svgText(L+w/2,H-6,opt.xlabel||'','class="axis-title" text-anchor="middle"');
 out+=`<text class="axis-title" transform="translate(16 ${T+h/2}) rotate(-90)" text-anchor="middle">${esc(opt.ylabel||'')}</text>`;
 out+=`<g clip-path="url(#clip-${id})">`;
 if(opt.bands)for(const b of opt.bands)out+=`<rect x="${x(b[0])}" y="${T}" width="${x(b[1])-x(b[0])}" height="${h}" fill="#c9d7cb" opacity=".35"/>`;
 if(opt.vline!==undefined)out+=`<line x1="${x(opt.vline)}" x2="${x(opt.vline)}" y1="${T}" y2="${T+h}" stroke="${muted}" stroke-dasharray="5 4"/>`;
 if(opt.hline!==undefined)out+=`<line x1="${L}" x2="${W-R}" y1="${y(opt.hline)}" y2="${y(opt.hline)}" stroke="${muted}" stroke-dasharray="5 4"/>`;
 series.forEach(s=>{let path=s.values.map((p,i)=>`${i?'L':'M'}${x(p[0]).toFixed(2)},${y(p[1]).toFixed(2)}`).join(' ');if(!s.dotsOnly)out+=`<path class="line" d="${path}" stroke="${s.color||cm}" ${s.dash?`stroke-dasharray="${s.dash}"`:''} opacity="${s.opacity??1}"/>`;if(s.dots||s.dotsOnly)out+=s.values.map(p=>`<circle cx="${x(p[0])}" cy="${y(p[1])}" r="4" fill="${s.color||cm}"><title>${esc(s.name||'')} ${fmt(p[0],3)}, ${fmt(p[1],3)}</title></circle>`).join('')});
 if(opt.cursor!==undefined)out+=`<line x1="${x(opt.cursor)}" x2="${x(opt.cursor)}" y1="${T}" y2="${T+h}" stroke="#1b3040" stroke-width="1.3"/>`;
 out+='</g>';svg.setAttribute('viewBox',`0 0 ${W} ${H}`);svg.innerHTML=out;
}
function renderPBC(){let p=+$('pbc').value,closest=Math.min(p,3-p);$('pbc-value').textContent=fmt(p,2)+' nm';let svg=$('pbc-plot'),w=Math.max(290,svg.clientWidth),xx=40+(w-80)*p/3;svg.setAttribute('viewBox',`0 0 ${w} 120`);svg.innerHTML=`<line x1="40" x2="${w-40}" y1="52" y2="52" stroke="#b5c3c3" stroke-width="3"/><circle cx="40" cy="52" r="12" fill="${cm}"/><circle cx="${xx}" cy="52" r="7" fill="${iod}"/>`+svgText(40,90,'0 nm · 离子','text-anchor="start"')+svgText(w-40,90,'3 nm · 周期边界','text-anchor="end"');$('pbc-text').innerHTML=`直接相减：<strong>${fmt(p,2)} nm</strong>；最近周期距离：<strong>${fmt(closest,2)} nm</strong>。${p>1.5?'沿边界绕过去更近。':'盒内直接路径已是最近路径。'}`}


function getHyd(ion,model){return D.hydration.filter(c=>c.summary.ion===ion&&(!model||c.summary.model===model))}

function renderShell(){const c=D.hydration.find(c=>c.id===$('shell-case').value),snap=c.snapshots[+$('shell-frame').value],scan=c.cutoff_scan[+$('shell-cutoff').value],cut=scan.cutoff_nm,inside=snap.waters.filter(p=>p[1]<cut);$('shell-cutoff-label').textContent=fmt(cut,3)+' nm';const series=[{values:inside.map((p,i)=>[p[1],.3+(i%5)*.1]),dotsOnly:true,color:cm},{values:snap.waters.filter(p=>p[1]>=cut).map((p,i)=>[p[1],.3+(i%5)*.1]),dotsOnly:true,color:'#a3b0ae'}];chart('shell-plot',series,{xmin:0,xmax:.6,ymin:0,ymax:1,xlabel:'真实离子–水氧距离 / nm',ylabel:'仅横轴有物理意义',xdigits:2,vline:cut,hideY:true});$('shell-instant').textContent=inside.length;$('shell-mean').textContent=fmt(scan.mean_CN,4);$('shell-default').textContent=fmt(cutoff(c.summary.ion),2);$('shell-reading').innerHTML=`<strong>${c.id}，${snap.time_ps} ps：</strong>蓝色 ${inside.length} 个氧在阈值内，灰色在外。整条轨迹同阈值均值 ${fmt(scan.mean_CN,4)}，与单帧整数计数不同。`}
function amideCases(){return D.amide.filter(c=>c.summary.ion===$('amide-ion').value&&c.summary.model===$('amide-model').value)}


const codeNotes=[
['把文献参数交给 OpenMM','<ul class="explain-list"><li>输入是 protocol.json 中当前模型 / 离子的两个数字。</li><li>0.2 同时完成“Rmin/2 → Rmin”和“Å → nm”。</li><li>除以 2**(1/6) 得到 σ，乘 4.184 转换 ε 单位。</li><li>返回二元组 (sigma, epsilon)。</li></ul>','CM / Mg 输入 1.360 Å、0.01020237 kcal/mol → σ≈0.242324 nm、ε≈0.042687 kJ/mol。之后才模拟，测得 RDF 峰 0.203 nm。'],
['测最近周期距离','<ul class="explain-list"><li>positions[others] 取一组氧坐标，减离子坐标得到位移。</li><li>减去整倍数盒长，折回最近周期方向。</li><li>norm(..., axis=-1) 把三维向量转成距离。</li></ul>','RDF、CN、酰胺距离都依赖这个函数。不处理跨边界问题，就会把很近的邻居误算很远。'],
['加水并替换离子 LJ 参数','<ul class="explain-list"><li>addSolvent 加 TIP3P 水并中和电荷，不额外加盐。</li><li>createSystem 建立 PME、截断和约束。</li><li>选择离子粒子，按 CM / IOD 改 σ、ε，保留 +2e。</li><li>代码还核对内置 CM 基线。</li></ul>','相同水与相互作用约定，仅比较离子参数。这里 1 nm 力计算 cutoff 与 CN 分析 cutoff 不同。'],
['从设置推导步数和帧数','<ul class="explain-list"><li>最小化缓解不合理接触，再分配初速度。</li><li>预平衡步数 = 100/0.002 = 50,000。</li><li>sample_steps = 1/0.002 = 500；n_frames = 2000/1 = 2000。</li><li>DCD / 热力学日志每 10 个采样间隔保存。</li></ul>','曲线从 1 到 2000 ps，共 2000 个保存样本。每个样本之间还有 500 个积分步。'],
['将每帧位置变为数组','<ul class="explain-list"><li>先前进 500 步，再取坐标和能量。</li><li>distances_nm[frame] 是一帧内全部水氧距离，每列始终对应同一个水氧。</li><li>另记时间、体积、温度、势能和最近 Cl⁻ 距离。</li><li>有酰胺则保存羰基氧距离；水合 case 的该指标为 NaN，表示不适用。</li></ul>','本页散点、CN 滑条、酰胺曲线都源于这些数组。NPZ 是原料，CSV / JSON 是统计摘要。'],
['True / False 变成配位数','<ul class="explain-list"><li>distances 行是帧，列是水氧。</li><li>distances &lt; cutoff 逐格判断，等于阈值不算壳内。</li><li>sum(axis=1) 数每行的 True，得到每帧 CN。</li><li>再把这些 CN 求均值，得到一条轨迹的平均配位数。</li></ul>','Mg 平均水 CN 为 6，Ca 约 8。酰胺氧不在 water 距离矩阵内，不会混入水配位数。'],
['RDF 为何不是普通直方图','<ul class="explain-list"><li>越远的球壳体积越大，均匀水也会有更多氧。</li><li>searchsorted 把距离分到 bin。</li><li>volumes/nwater 是每帧水数密度的倒数，压控时它会变化。</li><li>4π/3 × diff(edges**3) 算球壳体积。</li><li>加权计数除帧数与壳体积得 g(r)，累积 CN 则用未加密度权重的计数。</li></ul>','g(r) 表示相对平均密度的局部水氧密度。峰位置、峰高度、壳内距离均值是不同指标。'],
['用一行代码算接触帧比例','<ul class="explain-list"><li>amide 是 2000 个羰基氧距离。</li><li>amide &lt; cutoff 得到 2000 个是/否。</li><li>对布尔数组求均值等于接触帧数 / 总帧数。</li><li>之后也比较前后半段均距，没有自由能或速率拟合。</li></ul>','0% / 100% 描述有限时保存帧。不同起点的分歧说明采样依赖性，不是已求得平衡结合概率。']];
function renderCode(){let i=+$('code-step').value,s=SNIPPETS[i];$('code-title').textContent=codeNotes[i][0];$('code-source').textContent=`${s.file} · 第 ${s.start}–${s.start+s.lines.length-1} 行`;$('code-source').href='../'+s.file;$('code-content').innerHTML=s.lines.map((line,j)=>`<span class="code-line"><span class="line-no">${s.start+j}</span>${esc(line)}</span>`).join('');$('code-explain').innerHTML=codeNotes[i][1];$('code-result').textContent=codeNotes[i][2]}
const quizzes=[
['CM 换成 IOD，主要改了什么？',['把镁换成钙','离子 LJ 参数','水模型'],1,'两种离子都分别比较 CM / IOD，水仍为 TIP3P，电荷仍为 +2e。'],
['σ≈0.2423 nm 能当成 Mg 的 RDF 峰吗？',['可以，都是距离','不可以，前者是输入，后者是模拟结构指标'],1,'混合 LJ 参数、静电和周围溶液共同影响结构，实际 CM/Mg 峰约 0.203 nm。'],
['2,000 帧是否等于 2,000 个独立重复？',['是，时间不同','否，连续帧相关，每条件有 3 条重复轨迹'],1,'帧是同一部电影的连续画面。重复间差异用三条 seeded trajectories 描述。'],
['Mg 的 contact 起点 100%，separated 起点 0%，最直接说明？',['已证明平衡亲和力高','有限时结果仍依赖初态','永远不能自行结合'],1,'2 ns 内起点记忆没有消失，不能由此得到平衡亲和力或无限时间行为。'],
['CN cutoff 调大，配位数变化，意味着？',['原子被滑条移了位置','把更远的水氧纳入计数','运行了新的模型'],1,'这里只重算保存距离的定义，未改变轨迹或运行 OpenMM。'],
['项目校验通过证明什么？',['模型准确且充分平衡','被检查的数据、元数据和派生结果一致'],1,'一致性保证可追踪、可复算；准确性和充分采样需要各自的证据。']];
let answers={};

function runSelfCheck(){const check=(ok,msg)=>{if(!ok)throw Error(msg)};check(D.hydration.length===12&&D.amide.length===24,'轨迹数量错误');check([.3,.299].filter(x=>x<.3).length===1,'严格阈值错误');check(sd([1,2,3])===1,'样本SD错误');let zero=0;for(const c of D.amide){check(c.distance.length===2000,'酰胺帧数错误');check(c.distance.every(Number.isFinite),'距离无效');let f=c.distance.filter(x=>x<cutoff(c.summary.ion)).length/2000;check(Math.abs(f-c.summary.amide_fraction_within_diagnostic_water_cutoff)<1e-12,'接触比例不匹配 '+c.id);if(c.summary.initial_state==='separated'&&f===0)zero++}check(zero===12,'分离起点计数错误');for(const c of D.hydration){check(c.cn.length===2000,'CN帧数错误');check(Math.abs(mean(c.cn)-c.summary.mean_CN)<1e-12,'CN均值不匹配');check(c.rdf.length===D.radii.length,'RDF长度错误');check(c.cutoff_scan.length===9,'扫描长度错误');check(Math.abs(c.cutoff_scan[4].mean_CN-c.summary.mean_CN)<1e-12,'默认半径不匹配');for(const snap of c.snapshots){let n=snap.waters.filter(p=>p[1]<cutoff(c.summary.ion)).length;check(n===c.cn[snap.time_ps-1],'快照CN错误')}}return '通过：36 条轨迹、48,000 个酰胺距离、24,000 个水合 CN、84 个快照及接触比例一致。'}

// Browser teaching calculations; Python remains the source of scientific summaries.
function pairParameters(ion,model,scaleS=1,scaleE=1){
 const p=D.protocol.models[model][ion],w=D.water_oxygen;
 const sigma=.2*p.rmin_half_angstrom/2**(1/6)*scaleS,epsilon=4.184*p.epsilon_kcal_mol*scaleE;
 return {sigma,epsilon,pairSigma:(sigma+w.sigma_nm)/2,pairEpsilon:Math.sqrt(epsilon*w.epsilon_kj_mol)};
}
function lj(r,s,e){return 4*e*((s/r)**12-(s/r)**6)}
function renderParams(){
 const ion=$('param-ion').value,model=$('param-model').value,p=D.protocol.models[model][ion];
 const base=pairParameters(ion,model),a=+$('sigma-scale').value,b=+$('epsilon-scale').value,edit=pairParameters(ion,model,a,b),w=D.water_oxygen;
 $('sigma-label').textContent=fmt(a,2);$('epsilon-label').textContent=fmt(b,2);
 $('param-table').innerHTML=table(['参数对象','σ / nm 或 Rmin/2 / Å','ε / kJ mol⁻¹ 或 kcal mol⁻¹'],[
 ['Amber 离子原始输入',`Rmin/2=${p.rmin_half_angstrom} Å`,`${p.epsilon_kcal_mol} kcal/mol`],
 ['OpenMM 离子原始值',fmt(base.sigma,9)+' nm',fmt(base.epsilon,9)+' kJ/mol'],
 ['实际 TIP3P 水氧',fmt(w.sigma_nm,9)+' nm',fmt(w.epsilon_kj_mol,9)+' kJ/mol'],
 ['原始离子–水氧成对参数',fmt(base.pairSigma,9)+' nm',fmt(base.pairEpsilon,9)+' kJ/mol']]);
 $('param-explain').textContent=`${ion} / ${model}：下图先把离子参数与水氧按 Lorentz–Berthelot 混合。水氧 q=${w.charge_e}e；离子 q=+2e。`;
 const well=$('potential-view').value==='well';
 chart('lj-plot',[{values:Array.from({length:401},(_,i)=>{const r=.18+i*.001;return[r,lj(r,base.pairSigma,base.pairEpsilon)]}),color:cm},
 {values:Array.from({length:401},(_,i)=>{const r=.18+i*.001;return[r,lj(r,edit.pairSigma,edit.pairEpsilon)]}),color:iod,dash:'7 4'}],
 {xmin:well?.25:.18,xmax:well?.5:.58,ymin:well?-.8:-1,ymax:well?.8:20,xlabel:'离子–TIP3P 水氧距离 r / nm',ylabel:'仅 LJ 势能 / kJ mol⁻¹',xdigits:2,hline:0,ydigits:well?2:1});
 $('potential-reading').textContent=`原始成对 LJ 最低点 r=${fmt(2**(1/6)*base.pairSigma,6)} nm，U=${fmt(-base.pairEpsilon,6)} kJ/mol。教学探索最低点 r=${fmt(2**(1/6)*edit.pairSigma,6)} nm，U=${fmt(-edit.pairEpsilon,6)} kJ/mol。${a===1&&b===1?'当前两线重合；调整倍率后比较。':'倍率已改变：这是教学探索，没有新的 MD 结果。'}`;
}
function renderRDF(){
 const ion=$('rdf-ion').value,mode=$('rdf-view').value,model=$('rdf-model').value;
 const all=getHyd(ion).filter(c=>model==='both'||c.summary.model===model),series=[];
 for(const m of ['CM','IOD']){
  const rows=all.filter(c=>c.summary.model===m);if(!rows.length)continue;
  if(mode==='mean')series.push({values:D.radii.map((r,i)=>[r,mean(rows.map(c=>c.rdf[i]))]),color:m==='CM'?cm:iod,dash:m==='IOD'?'7 3':''});
  else rows.filter(c=>mode==='replicas'||c.summary.replica===+mode).forEach(c=>series.push({values:D.radii.map((r,j)=>[r,c.rdf[j]]),color:m==='CM'?cm:iod,dash:['','7 3','2 3'][c.summary.replica-1]}));
 }
 chart('rdf-plot',series,{xmin:0,xmax:.6,ymin:0,ymax:Math.max(...series.flatMap(s=>s.values.filter(p=>p[0]<=.6).map(p=>p[1])))*1.08,xlabel:'离子–水氧距离 r / nm',ylabel:'g(r)（无量纲）',xdigits:2,vline:ion==='Mg'?.209:.246,bands:ion==='Mg'?[[.205,.213]]:[]});
 const rows=all.filter(c=>['mean','replicas'].includes(mode)||c.summary.replica===+mode);
 $('rdf-table').innerHTML=table(['真实案例','RDF 峰 / nm','条件均距 / nm','平均 CN'],rows.map(c=>[c.id,fmt(c.summary.rdf_peak_nm,3),fmt(c.summary.mean_shell_distance_nm,6),fmt(c.summary.mean_CN,4)]));
 $('rdf-group-table').innerHTML=table(['组（n=3）','三条峰均值 ± SD / nm','壳内均距均值 / nm','CN 均值 ± SD'],D.groups.filter(g=>g.ion===ion&&(model==='both'||g.model===model)).map(g=>[g.model,`${fmt(g.rdf_peak_nm_mean,4)} ± ${fmt(g.rdf_peak_nm_SD,4)}`,fmt(g.mean_shell_distance_nm_mean,6),`${fmt(g.mean_CN_mean,4)} ± ${fmt(g.mean_CN_SD,4)}`]));
 $('rdf-reading').textContent=ion==='Mg'?'Mg：CM→IOD 的峰均值向外移动 0.0040 nm，平均 CN 均为 6。三条峰落在相同 bin，并不消除分箱、模型和慢过程的不确定性。':'Ca：CM→IOD 的峰均值向内移动 0.0040 nm。CM/IOD 平均 CN 为 8.0420/8.0025；壳内条件均距均比各自 RDF 峰更大，因为它平均的是全部壳内距离。';
}
const shellBase=renderShell;
renderShell=function(){shellBase();const c=D.hydration.find(c=>c.id===$('shell-case').value),cut=+$('shell-cutoff').value;
 $('cn-distribution').innerHTML=table(['默认 cutoff 下 CN','全部 2000 帧占比','帧数'],Object.entries(c.summary.CN_distribution).map(([n,f])=>[n,pct(f),Math.round(f*2000)]));
 $('shell-definition-note').textContent=cut===4?'当前使用默认 cutoff，下面的正式分布与滑块定义相同。':'滑块已改变定义；下面仍列默认 cutoff 的正式 CN 分布，不能用它代替新阈值下的分布。';
 const before=c.snapshots.find(s=>s.time_ps===167),after=c.snapshots.find(s=>s.time_ps===168),rc=cutoff(c.summary.ion),a=before.waters.filter(w=>w[1]<rc).map(w=>w[0]),b=after.waters.filter(w=>w[1]<rc).map(w=>w[0]);
 $('identity-example').textContent=`${c.id}：167 ps 壳内水列 ID [${a.join(', ')}]；168 ps [${b.join(', ')}]。离开 [${a.filter(i=>!b.includes(i)).join(', ')||'无'}]，进入 [${b.filter(i=>!a.includes(i)).join(', ')||'无'}]。阈值 ${fmt(rc,2)} nm；这里只比较两个保存帧。`;
};
function timeBounds(value){return value==='first'?[0,1000]:value==='second'?[1000,2000]:[0,2000]}
function blocks(values,lo,hi,size){const result=[];for(let i=lo;i<hi;i+=size){const end=Math.min(i+size,hi);result.push({time:(i+1+end)/2,mean:mean(values.slice(i,end)),n:end-i})}return result}
function renderAmide(){
 const cases=amideCases(),rep=+$('amide-rep').value,pair=['contact','separated'].map(st=>cases.find(c=>c.summary.replica===rep&&c.summary.initial_state===st)),ion=$('amide-ion').value,cut=cutoff(ion),ti=+$('amide-time').value;
 const [lo,hi]=timeBounds($('amide-window').value);
 chart('amide-plot',pair.map((c,i)=>({values:c.distance.slice(lo,hi).map((d,j)=>[(lo+j+1)/1000,d]),color:i===0?cm:iod,dash:i===0?'':'6 3'})),{xmin:lo/1000,xmax:hi/1000,ymin:0,ymax:Math.max(...pair.flatMap(c=>c.distance))*1.05,xlabel:'生产时间 / ns',ylabel:'离子–羰基氧距离 / nm',hline:cut,cursor:(ti+1)/1000});
 $('amide-time-label').textContent=`${ti+1} ps = ${fmt((ti+1)/1000,3)} ns${ti<lo||ti>=hi?'（游标位于所选窗口外）':''}`;
 $('amide-instant').innerHTML=table(['起点','当前距离 / nm',`是否 < ${fmt(cut,2)} nm`,'选定窗口接触比例','全程接触帧'],pair.map(c=>[c.summary.initial_state,fmt(c.distance[ti],5),c.distance[ti]<cut?'是 · 接触':'否',pct(c.distance.slice(lo,hi).filter(d=>d<cut).length/(hi-lo)),`${Math.round(c.summary.amide_fraction_within_diagnostic_water_cutoff*2000)} / 2000 (${pct(c.summary.amide_fraction_within_diagnostic_water_cutoff)})`]));
 $('contact-bars').innerHTML=table(['初态 / replica','全程接触比例','全程平均水 CN','全程羰基氧均距 / nm'],['contact','separated'].flatMap(st=>[1,2,3].map(r=>{const c=cases.find(c=>c.summary.initial_state===st&&c.summary.replica===r);return[`${st} · r${r}`,pct(c.summary.amide_fraction_within_diagnostic_water_cutoff),fmt(c.summary.mean_CN,4),fmt(c.summary.mean_amide_distance_nm,5)]})));
 const f=cases.filter(c=>c.summary.initial_state==='contact').map(c=>c.summary.amide_fraction_within_diagnostic_water_cutoff);
 $('amide-reading').textContent=`${ion} / ${$('amide-model').value}：contact 三条全程比例 ${f.map(pct).join('、')}；均值 ${pct(mean(f))}，SD ${fmt(sd(f)*100,2)} 个百分点。separated 三条均为 0%。先解释初态与重复差异，再考虑均值；这些比例不构成平衡亲和力排序。`;
 const prediction=$('amide-predict').selectedIndex;
 $('prediction-feedback').textContent=prediction===0?'先预测，再比较 Ca / CM / r3 的 297→298 ps 和两个半段；操作后用一句话说明初态影响。':prediction===1?'这个预测与观察相符：两种起点的结果在2 ns中仍不同；同一条件的三条重复也有差异。请用表中的具体比例支持解释。':'温度稳定只约束热运动诊断，不说明配体已反复进出、不同初态已探索相同状态。比较 contact 与 separated 的真实曲线，修正你的解释。';
}
function renderDiag(){
 const c=D.hydration.find(c=>c.id===$('diag-case').value),s=c.summary,key=$('diag-metric').value,series=key==='CN'?c.cn:c[key],[lo,hi]=timeBounds($('diag-window').value),size=+$('diag-block').value,bs=blocks(series,lo,hi,size);
 const names={CN:'水配位数 CN',temperature_K:'温度 / K',water_number_density_nm3:'水氧数密度 / nm⁻³'};
 chart('diag-plot',[{values:bs.map(b=>[b.time/1000,b.mean]),dots:true,color:s.model==='CM'?cm:iod}],{xmin:lo/1000,xmax:hi/1000,xlabel:'生产时间 / ns',ylabel:names[key],ydigits:key==='temperature_K'?1:3,hline:key==='temperature_K'?298.15:undefined});
 const h=s.halves[key];$('diag-table').innerHTML=table(['前 1 ns 均值','后 1 ns 均值','后 − 前','当前窗口均值'],[[fmt(h.first,4),fmt(h.second,4),fmt(h.second-h.first,4),fmt(mean(series.slice(lo,hi)),4)]]);
 const rows=getHyd(s.ion,s.model),values=rows.map(r=>mean((key==='CN'?r.cn:r[key]).slice(lo,hi)));
 $('diag-replicas').innerHTML=table(['相同条件的真实轨迹','相同窗口指标均值'],rows.map((r,i)=>[r.id,fmt(values[i],5)]).concat([['三条均值 ± 样本 SD',`${fmt(mean(values),5)} ± ${fmt(sd(values),5)}`]]));
 $('diag-block-table').innerHTML=table(['块内保存时间 / ps','帧数','块均值'],bs.map((b,i)=>[`${lo+i*size+1}–${Math.min(lo+(i+1)*size,hi)}`,b.n,fmt(b.mean,5)]));
 $('diag-reading').textContent=`窗口 ${lo+1}–${hi} ps，按 ${size} 帧分块；统计来自每个完整保存值，不从图上抽点。全程共见 ${s.distinct_waters_seen_in_shell} 个不同壳内水氧，身份改变间隔 ${s.shell_identity_changed_intervals}/${s.identity_observation_intervals}。该轨迹最近 Cl⁻ <0.40 nm 的保存帧占比为 ${pct(s.chloride_contact_fraction_lt_0p4nm)}。纵轴自动缩放；平坦或小幅变化不构成慢过程收敛证明。`;
}
function histToy(pooled=false){
 const counts=[[1,1],[2,0]],density=[2,1],shells=[4*Math.PI/3*.2**3,4*Math.PI/3*(.4**3-.2**3)];
 const weighted=counts[0].map((_,j)=>mean(counts.map((r,i)=>r[j]/density[i])));
 return {counts,density,shells,weighted,g:shells.map((v,j)=>(pooled?mean(counts.map(r=>r[j]))/mean(density):weighted[j])/v),cn:[1.5,2]};
}
function renderHist(){const step=+$('hist-step').value,pooled=$('hist-mode').value==='pooled',t=histToy(pooled);const lines=[
 '<p><strong>输入：</strong>距离 shape=(2,2)，单位 nm；体积 shape=(2,)，单位 nm³。</p><pre>distances = [[0.10, 0.30],\n             [0.10, 0.15]]\nvolumes = [1.0, 2.0]\nedges = [0.0, 0.2, 0.4]</pre>',
 '<p><strong>分箱：</strong>[0,0.2) 和 [0.2,0.4)。0.30 落外壳；其他三个距离落内壳。</p><pre>counts (帧 × bin) = [[1, 1],\n                    [2, 0]]</pre>',
 `<p><strong>球壳体积：</strong>[${t.shells.map(x=>fmt(x,7)).join(', ')}] nm³；密度 [2,1] nm⁻³。</p><pre>inverse_density[:, None] = [[0.5], [1.0]]\nbroadcast → [[0.5, 0.5], [1.0, 1.0]]</pre>`,
 `<p><strong>${pooled?'反例':'项目定义'}：</strong>${pooled?'平均计数 [1.5,0.5] / 平均密度 1.5 → [1,0.333333] nm³。它不等于逐帧密度加权的平均。':'每帧 counts / density → [[0.5,0.5],[2,0]] nm³；对帧平均 → [1.25,0.25] nm³。'}</p>`,
 table(['bin / nm','g（无量纲）','累计平均 CN'],[[ '[0,0.2)',fmt(t.g[0],6),'1.5'],['[0.2,0.4)',fmt(t.g[1],6),'2.0']])+`<p>${pooled?'这个 g 使用错误的汇总顺序，与项目定义不同。累计 CN 仍直接来自原始计数，未随错误的密度归一化变化。':'累计 CN 来自计数平均 [1.5,0.5] 的累积；它没有用密度权重。'}</p>`];
 $('hist-output').innerHTML=lines.slice(0,step+1).join('');}
function arrayToy(){const xyz=[[[2.9,0,0],[0,.2,0],[.4,0,0]],[[2.9,0,0],[0,.4,0],[.2,0,0]]];const delta=xyz.map(frame=>frame.map(v=>v.map(x=>x-3*Math.round(x/3)))),d=delta.map(frame=>frame.map(v=>Math.hypot(...v))),mask=d.map(row=>row.map(x=>x<.3));return{xyz,delta,d,mask,cn:mask.map(row=>row.filter(Boolean).length),columns:mask[0].map((_,j)=>mask.reduce((s,row)=>s+Number(row[j]),0)),flips:mask[0].filter((v,j)=>v!==mask[1][j]).length,conditional:mean(d.flat().filter(v=>v<.3))}}
function renderArray(){const t=arrayToy(),i=+$('array-step').value;const coordinateText=a=>a.map((f,k)=>`frame ${k}:\n`+f.map((v,j)=>`  water ${j}: [${v.map(x=>fmt(x,3)).join(', ')}]`).join('\n')).join('\n\n');const rows=[
 ['输入：两个保存帧，每帧三个水氧','xyz.shape = (2,3,3)；离子坐标固定 (0,0,0)，box=[3,3,3] nm。',coordinateText(t.xyz),'positions[others] - positions[index]'],
 ['折回最近镜像','2.9/3≈0.9667，round→1；2.9−3=−0.1。每帧位移 shape 仍是 (3,3)。',coordinateText(t.delta),'delta -= box * np.round(delta / box)'],
 ['取向量长度','最后一维的 xyz 分量被 norm 消去，输出 d.shape=(2,3)，行是帧、列是固定水身份。',t.d.map((r,k)=>`frame ${k}: [${r.map(x=>fmt(x,3)).join(', ')}]`).join('\n'),'d = np.linalg.norm(delta, axis=-1)'],
 ['严格阈值变成布尔矩阵','0.30 nm 内为 True；比较不改变 shape=(2,3)。',t.mask.map((r,k)=>`frame ${k}: [${r.map(x=>x?'True ':'False').join(', ')}]`).join('\n'),'mask = d < 0.30'],
 ['沿哪一轴求和？','axis=1 消去水列，得到每帧 CN；axis=0 消去时间行，得到每个水在壳内的帧数。',`CN shape=(2,): ${JSON.stringify(t.cn)}\n每水帧数 shape=(3,): ${JSON.stringify(t.columns)}`,'cn = mask.sum(axis=1)\nper_water = mask.sum(axis=0)'],
 ['数值稳定不等于身份不变',`两帧 CN 均为 2，但第2列离开、第3列进入，共 ${t.flips} 次成员翻转。壳内四个距离的均值 ${fmt(t.conditional,6)} nm。`,'d[mask] shape=(4,) → [0.1, 0.2, 0.1, 0.2]','d[mask].mean()\nnp.count_nonzero(mask[1:] != mask[:-1], axis=1)']];
 const row=rows[i];$('array-output').innerHTML=`<h4>${i+1}/6 · ${row[0]}</h4><p>${row[1]}</p><pre>${esc(row[2])}</pre><p>对应的简化 Python 操作：</p><pre>${esc(row[3])}</pre>`;$('array-prev').disabled=i===0;$('array-next').disabled=i===5;
}
function renderQuiz(){
 $('quiz-list').innerHTML=quizzes.map((q,i)=>`<div class="question"><p>${i+1}. ${q[0]}</p><div class="answers" role="group" aria-label="第 ${i+1} 题">${q[1].map((a,j)=>`<button type="button" data-q="${i}" data-a="${j}" aria-pressed="false">${String.fromCharCode(65+j)}. ${a}</button>`).join('')}</div><div id="feedback-${i}" class="feedback" aria-live="polite"></div></div>`).join('');
 $('quiz-score').textContent='已尝试 0/6 题。';
}
function numericFeedback(){const value=$('numeric-answer').value,n=Number(value);$('numeric-feedback').textContent=value===''?'先输入你算出的百分数。':Math.abs(n-14.85)<.005?'正确：297÷2000×100=14.85%。它是这段 2 ns 保存帧的接触占比，不能把 297 帧称作 297 次结合事件。':Math.abs(n-.1485)<.001?'你算出了比例 0.1485；题目要求百分数，还需要乘 100。':'先让接触帧数除以总帧数，再乘 100。297/2000，而不是2000/297；分母也不是3条重复。'}
function extendedCheck(){
 const check=(b,m)=>{if(!b)throw Error(m)},near=(a,b,tol=1e-10)=>Math.abs(a-b)<tol;
 runSelfCheck();const p=pairParameters('Mg','CM');check(near(p.sigma,.2423244513341723),'单位转换');check(near(p.pairSigma,.27869984599584235)&&near(p.pairEpsilon,.16476463653334547),'与保存体系水氧参数的混合对照');check(near(lj(2**(1/6)*.3,.3,.5),-.5),'LJ 势阱');
 const t=histToy();check(near(t.g[0]*t.shells[0],1.25)&&near(t.g[1]*t.shells[1],.25),'RDF 逐帧归一化');check(!near(t.g[0],histToy(true).g[0]),'反例未区分');
 const a=arrayToy();check(near(a.d[0][0],.1)&&a.cn.join()==='2,2'&&a.columns.join()==='2,1,1'&&a.flips===2&&near(a.conditional,.15),'距离/轴/身份教学样例');
 for(const c of D.hydration){check(c.cutoff_scan.every((r,i,s)=>!i||r.mean_CN>=s[i-1].mean_CN),'CN 必须单调');for(const key of ['CN','temperature_K','water_number_density_nm3']){const v=key==='CN'?c.cn:c[key],bs=blocks(v,0,2000,250);check(bs.every((b,i)=>near(b.mean,c.summary.block_means[i][key],1e-9)),'Python/JS 分块不匹配 '+c.id+' '+key)}}
 check(near(mean([.999,.842,.1485])*100,66.31666666666666),'重复均值');check(near(sd([.999,.842,.1485])*100,45.257439535763986),'重复SD');
 return '通过：36 条完整轨迹、84 个快照、单位/混合/LJ、RDF 手算、数组与身份、12×3 指标分块和阈值单调性。';
}
window.tutorialCheck=extendedCheck;
const panels=[...document.querySelectorAll('main>section')],initialized=new Set();let currentPanel='question',keyboardTarget=null;
const defaults={potential:{'param-ion':'Mg','param-model':'CM','potential-view':'well','sigma-scale':'1','epsilon-scale':'1'},rdf:{'rdf-ion':'Mg','rdf-model':'both','rdf-view':'mean'},shell:{'shell-case':'CM_Ca_r1_hydrated','shell-frame':'0','shell-cutoff':'4'},amide:{'amide-ion':'Ca','amide-model':'CM','amide-rep':'3','amide-time':'296','amide-window':'all','amide-predict':'unsure'},diag:{'diag-case':'CM_Ca_r1_hydrated','diag-metric':'CN','diag-window':'all','diag-block':'250'},hist:{'hist-step':'0','hist-mode':'frame'},array:{'array-step':'0'},pbc:{pbc:'2.9'},code:{'code-step':'0'}};
const renders={potential:renderParams,rdf:renderRDF,shell:renderShell,amide:renderAmide,diag:renderDiag,hist:renderHist,array:renderArray,pbc:renderPBC,code:renderCode};
function reset(which){if(defaults[which]){Object.entries(defaults[which]).forEach(([id,v])=>$(id).value=v);renders[which]()}else if(which==='quiz'){answers={};renderQuiz()}else if(which==='numeric'){$('numeric-answer').value='';$('numeric-feedback').textContent=''}else if(which==='draft')$('explain-draft').value='';}
function drawActive(){if(currentPanel==='concepts')renderPBC();if(currentPanel==='models')renderParams();if(currentPanel==='metrics')renderHist();if(currentPanel==='results'){renderRDF();renderShell();renderAmide();renderDiag()}if(currentPanel==='python'){renderArray();renderCode()}}
function initialize(id){if(initialized.has(id))return;initialized.add(id);
 if(id==='question')$('study-inventory').textContent=`已核对 ${D.hydration.length} 条水合、${D.amide.length} 条乙酰胺轨迹，全部完成；每条 ${D.protocol.production_ps/1000} ns，正式生产共 ${D.timing.production_ns} ns。每条保留 ${D.time_ps.length} 个距离样本，首帧 ${D.time_ps[0]} ps，末帧 ${D.time_ps.at(-1)} ps。`;
 if(id==='concepts')$('concentration-note').textContent=`从36条轨迹的瞬时体积计算，各轨迹平均阳离子浓度为 ${fmt(Math.min(...D.hydration.concat(D.amide).map(c=>c.mean_cation_molar)),4)}–${fmt(Math.max(...D.hydration.concat(D.amide).map(c=>c.mean_cation_molar)),4)} M；Cl⁻ 浓度是其两倍。这是有限浓度体系。`;
 if(id==='models'){const v=D.conversion,errors=v.geometries.flatMap(g=>Object.values(g.by_force));$('conversion-result').textContent=`现存验证记录：${v.geometries.length} 个结构全部通过；最大逐力项能量误差 ${Math.max(...errors.map(e=>e.absolute_energy_error_kj_mol)).toExponential(3)} kJ/mol，最大原子力分量误差 ${Math.max(...errors.map(e=>e.max_absolute_force_error_kj_mol_nm)).toExponential(3)} kJ/(mol·nm)。`;}
 if(id==='design')$('run-inventory').innerHTML=table(['campaign / case','seed','水数','平台','最短盒边 / nm','生产 / ns'],D.hydration.concat(D.amide).map(c=>[c.summary.campaign+'/'+c.id,c.run.seed,c.run.n_water,c.run.platform,fmt(c.run.min_box_length_nm,5),c.summary.production_ps/1000]));
 if(id==='results'){addOptions($('shell-case'),D.hydration);addOptions($('diag-case'),D.hydration);$('shell-case').value='CM_Ca_r1_hydrated';$('diag-case').value='CM_Ca_r1_hydrated';$('amide-rep').value='3';$('amide-time').value='296'}
 if(id==='practice'){renderQuiz();const tm=D.timing;$('timing-note').textContent=`正式生产速度 ${fmt(tm.min_ns_per_day,1)}–${fmt(tm.max_ns_per_day,1)} ns/day；总生产量除以生产耗时之和 ${fmt(tm.production_ns_per_day,1)} ns/day；各 case 总耗时相加 ${fmt(tm.case_wall_seconds/3600,2)} 小时。`;$('provenance-info').textContent=`${D.provenance.sources.length} 个路径与 SHA-256 已内嵌。网页无需读取旁边的 JSON；导出时与原始数据、正式汇总核对。`;$('original-figures').innerHTML=FIGURES.map(f=>`<figure><img class="figure" loading="lazy" src="${f.data}" alt="${esc(f.caption)}"><figcaption>${f.caption}</figcaption></figure>`).join('')}
}
function route(){
 let raw;try{raw=decodeURIComponent(location.hash.slice(1)||'question')}catch{raw='question'}const target=$(raw),panel=target?.closest('main>section')||$('question');currentPanel=panel.id;
 panels.forEach(p=>p.hidden=p!==panel);document.querySelectorAll('[role=tab]').forEach(t=>{const active=t.getAttribute('aria-controls')===panel.id;t.setAttribute('aria-selected',active);t.tabIndex=active?0:-1});
 initialize(panel.id);drawActive();const headings=[...panel.querySelectorAll('h3[id]')].filter(h=>!h.closest('details'));
 $('chapter-toc').innerHTML=headings.map(h=>`<a href="#${h.id}" ${h.id===raw?'aria-current="location"':''}>${esc(h.textContent)}</a>`).join('');$('chapter-name').textContent=panel.dataset.title;
 const i=panels.indexOf(panel);$('prev-chapter').hidden=i===0;$('next-chapter').hidden=i===panels.length-1;if(i>0){$('prev-chapter').href='#'+panels[i-1].id;$('prev-chapter').textContent='← '+panels[i-1].dataset.title}if(i<panels.length-1){$('next-chapter').href='#'+panels[i+1].id;$('next-chapter').textContent=panels[i+1].dataset.title+' →'}
 if(target&&panel.contains(target)){for(let p=target.parentElement;p&&p!==panel;p=p.parentElement)if(p.tagName==='DETAILS')p.open=true;target.scrollIntoView({block:'start'});}else window.scrollTo(0,0);if(keyboardTarget){keyboardTarget.focus({preventScroll:true});keyboardTarget=null;}
}
function fitHeader(){const height=document.querySelector('.book-header').offsetHeight+12;document.documentElement.style.setProperty('--header-height',height+'px')}
new ResizeObserver(fitHeader).observe(document.querySelector('.book-header'));
document.addEventListener('click',e=>{const r=e.target.closest('[data-reset]');if(r)reset(r.dataset.reset);const b=e.target.closest('button[data-q]');if(b){const i=+b.dataset.q,j=+b.dataset.a,q=quizzes[i];answers[i]=j;b.closest('.answers').querySelectorAll('button').forEach(el=>el.setAttribute('aria-pressed',String(el===b)));$('feedback-'+i).className='feedback '+(j===q[2]?'correct':'wrong');$('feedback-'+i).textContent=(j===q[2]?'正确。':'再想一步。')+q[3];$('quiz-score').textContent=`已尝试 ${Object.keys(answers).length}/6 题，当前答对 ${Object.entries(answers).filter(([i,a])=>a===quizzes[i][2]).length} 题。`;}});
[['pbc','input',renderPBC],['param-ion','change',renderParams],['param-model','change',renderParams],['sigma-scale','input',renderParams],['epsilon-scale','input',renderParams],['rdf-ion','change',renderRDF],['rdf-model','change',renderRDF],['rdf-view','change',renderRDF],['shell-case','change',renderShell],['shell-frame','change',renderShell],['shell-cutoff','input',renderShell],['amide-ion','change',renderAmide],['amide-model','change',renderAmide],['amide-rep','change',renderAmide],['amide-window','change',renderAmide],['amide-time','input',renderAmide],['diag-case','change',renderDiag],['diag-metric','change',renderDiag],['diag-window','change',renderDiag],['diag-block','change',renderDiag],['code-step','change',renderCode],['array-step','change',renderArray],['hist-step','change',renderHist],['hist-mode','change',renderHist]].forEach(([id,event,fn])=>$(id).addEventListener(event,fn));
$('array-prev').onclick=()=>{$('array-step').value=String(Math.max(0,+$('array-step').value-1));renderArray()};$('array-next').onclick=()=>{$('array-step').value=String(Math.min(5,+$('array-step').value+1));renderArray()};
$('potential-view').onchange=renderParams;$('amide-predict').onchange=renderAmide;
$('numeric-check').onclick=numericFeedback;$('numeric-answer').addEventListener('keydown',e=>{if(e.key==='Enter')numericFeedback()});
$('self-check').onclick=()=>{$('self-check-status').textContent=extendedCheck()};
document.querySelector('[role=tablist]').addEventListener('keydown',e=>{const tabs=[...document.querySelectorAll('[role=tab]')],i=tabs.indexOf(document.activeElement);if(i<0)return;let j=i;if(e.key==='ArrowRight')j=(i+1)%tabs.length;else if(e.key==='ArrowLeft')j=(i+tabs.length-1)%tabs.length;else if(e.key==='Home')j=0;else if(e.key==='End')j=tabs.length-1;else return;e.preventDefault();keyboardTarget=tabs[j];tabs[j].focus();const next=tabs[j].getAttribute('aria-controls');if(location.hash==='#'+next)keyboardTarget=null;else location.hash=next});
window.addEventListener('hashchange',route);let resizeTimer;window.addEventListener('resize',()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(drawActive,120)});
history.scrollRestoration='manual';fitHeader();route();extendedCheck();
