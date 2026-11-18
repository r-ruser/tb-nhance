import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {chromium} from 'playwright';
const dir=fileURLToPath(new URL('../Figures/',import.meta.url));const browser=await chromium.launch({headless:true});
for(const [num,files] of [[1,['Figure1_composite','Figure1D_thresholds']],[3,['Figure3_circular_forest','Figure3B_modules']]]){
 const page=await browser.newPage({deviceScaleFactor:3});
 const parts=await Promise.all(files.map(f=>fs.readFile(path.join(dir,f+'.svg'),'utf8')));
 const body=await page.evaluate(({parts,num})=>{
 const ns='http://www.w3.org/2000/svg';let top=0,root=document.createElementNS(ns,'svg');root.setAttribute('xmlns',ns);
 for(const [i,txt] of parts.entries()){
 const el=new DOMParser().parseFromString(txt,'image/svg+xml').documentElement;
 const [,,w,h]=el.getAttribute('viewBox').split(/\s+/).map(Number);
 el.setAttribute('x',0);el.setAttribute('y',top);el.setAttribute('width','800');el.setAttribute('height',800*h/w);
 if(i){const map=new Map([...el.querySelectorAll('[id]')].map(n=>[n.id,'sub_'+n.id]));for(const n of el.querySelectorAll('*'))for(const a of [...n.attributes]){if(a.name==='id')n.setAttribute('id',map.get(a.value));else if(a.value.includes('url(#'))n.setAttribute(a.name,a.value.replace(/url\(#([^)]+)\)/g,(m,id)=>'url(#'+(map.get(id)||id)+')'))}}
 root.append(el);top+=800*h/w+4;
 } root.setAttribute('viewBox','0 0 800 '+Math.ceil(top));root.setAttribute('width',800);root.setAttribute('height',Math.ceil(top));return [new XMLSerializer().serializeToString(root),Math.ceil(top)];
 },{parts,num});
 const out=path.join(dir,'Figure'+num+'_complete');await fs.writeFile(out+'.svg',body[0]);await page.setContent('<style>*{margin:0}body{background:white}</style>'+body[0]);await page.locator('body>svg').screenshot({path:out+'.png'});await page.pdf({path:out+'.pdf',width:'800px',height:(body[1]+2)+'px',printBackground:true});await page.close()
}await browser.close();
