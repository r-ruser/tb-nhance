import fs from 'node:fs/promises';
import {chromium} from 'playwright';
const root='./EHP_R2_20261118';
const src=await fs.readFile(root+'/Figures/Figure1_flow.mmd','utf8');
const panels=await Promise.all(['Figure1_B_states','Figure1_C_nutrition'].map(n=>fs.readFile(root+'/Figures/'+n+'.svg','utf8')));
const browser=await chromium.launch({headless:true,args:['--no-sandbox']});
try {
 const page=await browser.newPage({viewport:{width:1440,height:1340},deviceScaleFactor:3});
 await page.setContent('<html><body style="margin:0;background:white"><div id="out"></div></body></html>');
 await page.evaluate(async()=>{window.mm=(await import('https://cdn.jsdelivr.net/npm/mermaid@11.12.0/dist/mermaid.esm.min.mjs')).default;window.mm.initialize({startOnLoad:false,securityLevel:'strict'});});
 const svg=await page.evaluate(async({src,panels})=>{
  const ns='http://www.w3.org/2000/svg',root=document.createElementNS(ns,'svg');
  root.setAttribute('xmlns',ns);root.setAttribute('viewBox','0 0 1400 1300');root.setAttribute('width','1400');root.setAttribute('height','1300');
  const bg=document.createElementNS(ns,'rect');bg.setAttribute('width','1400');bg.setAttribute('height','1300');bg.setAttribute('fill','white');root.append(bg);
  const flow=(await window.mm.render('flowStudy',src)).svg;
  function insert(text,x,y,w,h,prefix){const wrap=document.createElement('div');wrap.innerHTML=text;const e=wrap.querySelector('svg');if(!e)throw new Error('Missing SVG');
   if(prefix!=='f')for(const a of e.querySelectorAll('[id]')){const id=a.id;for(const el of e.querySelectorAll('*'))for(const at of [...el.attributes])if(at.value.includes('#'+id))el.setAttribute(at.name,at.value.replaceAll('#'+id,'#'+prefix+id));a.id=prefix+id;}
   for(const[k,v]of Object.entries({x,y,width:w,height:h,preserveAspectRatio:'xMidYMin meet'}))e.setAttribute(k,v);e.removeAttribute('style');root.append(document.importNode(e,true));}
  insert(flow,10,70,545,1140,'f');insert(panels[0],585,40,800,280,'b');insert(panels[1],585,370,800,865,'c');
  function text(x,y,t,size=26){const e=document.createElementNS(ns,'text');e.setAttribute('x',x);e.setAttribute('y',y);e.setAttribute('font-family','Arial');e.setAttribute('font-size',size);e.setAttribute('fill','#263746');e.textContent=t;root.append(e);}
  text(20,35,'A');text(590,35,'B');text(590,355,'C');
  document.querySelector('#out').append(root);return new XMLSerializer().serializeToString(root);
 },{src,panels});
 await fs.writeFile(root+'/Figures/Figure1_composite.svg',svg);
 await page.locator('#out>svg').screenshot({path:root+'/Figures/Figure1_composite.png'});
 await page.pdf({path:root+'/Figures/Figure1_composite.pdf',width:'1400px',height:'1300px',printBackground:true,margin:{top:'0',bottom:'0',left:'0',right:'0'}});
 console.log('Mermaid schematic and R quantitative panels composed.');
} finally {await browser.close();}
