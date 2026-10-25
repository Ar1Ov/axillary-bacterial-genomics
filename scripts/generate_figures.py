from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,Patch
import pandas as pd,numpy as np,fitz,zipfile,json
R=Path(__file__).resolve().parents[1]; ROOT=R; SRC=R/'data/figures'
import argparse,os,tempfile
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=R/'figures');parser.add_argument('--font-dir',type=Path);args=parser.parse_args()
OUT=args.output.resolve();OUT.mkdir(parents=True,exist_ok=True)
TMP=Path(tempfile.mkdtemp(prefix='axillary_figures_'))
# Arial can be supplied locally; no proprietary font is distributed.
if args.font_dir:
 for fp in args.font_dir.glob('*.ttf'):font_manager.fontManager.addfont(fp)
try:
 regular_font=font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False);family='Arial'
except ValueError:
 regular_font=font_manager.findfont(font_manager.FontProperties(family='Liberation Sans'));family=font_manager.FontProperties(fname=regular_font).get_name()
print('Figure font:',family)
plt.rcParams.update({'font.family':family,'font.size':7.5,'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
ORANGE=orange='#c76b22';GRAY=gray='#666666'
REFS=['DtpT','PepV','PatB','AgaA','TpdA','AecD'];SPEC=['Staphylococcus epidermidis','Staphylococcus haemolyticus','Staphylococcus hominis','Corynebacterium amycolatum','Corynebacterium axilliensis','Corynebacterium gottingense','Corynebacterium jamesii','Corynebacterium kefirresidentii','Corynebacterium marquesiae','Corynebacterium tuberculostearicum']
ff=fitz.Font(fontfile=regular_font)
def read(n):return pd.read_csv(SRC/n,sep='\t',dtype={'genome_id':str})
def abbr(s):return s[0]+'. '+s.split()[1]
def save(f,n,title):
 # Trim unused canvas, then append a compact, consistently inset title band.
 f.canvas.draw();f.savefig(TMP/f'raw{n}.pdf',bbox_inches='tight',pad_inches=.06)
 d=fitz.open(TMP/f'raw{n}.pdf');p=d[0]
 # Preserve true 6.5-inch full-figure width and keep the title near the content.
 w=468;inset=9;h=p.rect.height*(w-2*inset)/p.rect.width
 title=title if title.startswith('Figure ') else f'Figure {n}. {title}'
 title_x=45 if n==2 else inset
 lines=[];line=''
 for word in title.split():
  trial=(line+' '+word).strip()
  if ff.text_length(trial,fontsize=10)>w-title_x-inset and line:lines.append(line);line=word
  else:line=trial
 lines.append(line)
 out=fitz.open();q=out.new_page(width=w,height=h+13*len(lines)+14)
 q.show_pdf_page(fitz.Rect(inset,0,w-inset,h),d,0)
 q.insert_font(fontname='ArialTitle',fontfile=regular_font)
 for i,line in enumerate(lines):q.insert_text((title_x,h+13+i*13),line,fontname='ArialTitle',fontsize=10)
 out.save(OUT/f'Figure{n}.pdf',deflate=True)
 q.get_pixmap(dpi=400,alpha=False).pil_save(OUT/f'Figure{n}.jpg',quality=97,dpi=(400,400))
 q.get_pixmap(matrix=fitz.Matrix(1.6,1.6)).save(TMP/f'preview{n}.png')
 spans=[s for b in q.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['bbox'][1]>=h]
 assert spans and all(abs(s['size']-10)<.01 and abs(s['bbox'][0]-title_x)<.05 for s in spans)
 plt.close(f)




from Bio import Phylo
import textwrap
def read(name):return pd.read_csv(SRC/name,sep='\t')
def dumb(ax,a,b):
 for y,(x,z) in enumerate(zip(a,b)):
  ax.plot([x,z],[y-.055,y+.055],color='#ccc',lw=.8,zorder=1)
  ax.plot(z,y+.055,'o',mfc='white',mec=gray,ms=3.6,mew=.7,zorder=2)
  ax.plot(x,y-.055,'o',color=orange,ms=3.6,zorder=3)
def axis(ax,n):
 ax.set_xlim(-.04,1.04);ax.set_ylim(n-.5,-.5);ax.set_xticks([0,.25,.5,.75,1]);ax.tick_params(labelsize=6.5);ax.set_xlabel('Prevalence',fontsize=7);ax.grid(axis='x',color='#eee',lw=.6);ax.set_axisbelow(True)
def key():return [Line2D([],[],marker='o',ls='',mfc='white',mec=gray,label='Other Skin',ms=4),Line2D([],[],marker='o',ls='',color=orange,label='Axilla',ms=4)]
strict=read('staph_epi_strict_primary_cohort.tsv')['group'].value_counts(); same=read('staph_epi_prjna803478_final_cohort.tsv')['group'].value_counts()
d=read('staph_epi_prjna803478_same_study_fdr_axilla_candidates.tsv').sort_values(['OTHER_SKIN_prevalence','pgfam_id'])
f,aa=plt.subplots(1,2,figsize=(6.5,6.8),sharey=True);f.subplots_adjust(left=.345,right=.975,top=.90,bottom=.20,wspace=.22)
labels=[str(r.representative_product).replace(' of Caudovirales','').replace('FIG01109124: ','')+'\n'+r.pgfam_id for _,r in d.iterrows()]
labels=['\n'.join(textwrap.wrap(x.split('\n')[0],29))+'\n'+x.split('\n')[1] for x in labels]
for ax,a,b,title in [(aa[0],'AXILLARY_GROIN_prevalence','OTHER_SKIN_prevalence',f'A  Same-Study Comparison\n{same["AXILLARY_GROIN"]} Axillary/Groin / {same["OTHER_SKIN"]} Other Skin'),(aa[1],'strict_AXILLA_prevalence','strict_OTHER_SKIN_prevalence',f'B  Strict-Cohort Comparison\n{strict["AXILLA"]} Axilla / {strict["OTHER_SKIN"]} Other Skin')]:
 dumb(ax,d[a],d[b]);axis(ax,14);ax.set_title(title,fontsize=8,pad=12,loc='left',weight='bold');ax.set_xlabel('PGFam Prevalence')
aa[0].set_yticks(range(14),labels,fontsize=6.3)
h=key();h[1].set_label('Axillary/Groin (A); Axilla (B)');f.legend(handles=h,loc='lower center',bbox_to_anchor=(.61,.108),ncol=2,frameon=False,fontsize=7)
save(f,2,'Figure 2. The 14 Same-Study FDR-Significant PGFams Show Markedly Weaker Separation in the Strict Cohort')

d=read('cory_exact_site_stratified_results.tsv').set_index('pgfam_id');ids=['PGF_10347576','PGF_08044753','PGF_03285992','PGF_00267745','PGF_03133617'];labels=['ABC Antimicrobial-Peptide\nTransport ATPase','Putative Nisin-Resistance\nProtein','ECF Transporter\nATPase Component','ECF Transporter\nTransmembrane Component','AecD/MetC-Like β-Lyase\n(Known-Pathway Control)']
f=plt.figure(figsize=(6.5,7.6));aa=[f.add_axes([.30+i*.174,.55,.148,.33]) for i in range(3)]
f.text(.045,.94,'A  Candidate Prevalence Within Species',fontsize=8.5,fontweight='bold')
cohort=read('cory_exact_site_final_cohort.tsv'); ccounts=cohort.groupby(['species','group']).size()
for ax,sp in zip(aa,['amycolatum','kefirresidentii','tuberculostearicum']):
 counts=f'{ccounts[("Corynebacterium "+sp,"AXILLA")]} AX / {ccounts[("Corynebacterium "+sp,"OTHER_SKIN")]} OS'
 dumb(ax,d.loc[ids,sp+'_AXILLA_prevalence'],d.loc[ids,sp+'_OTHER_SKIN_prevalence']);axis(ax,5);ax.set_xticks([0,.5,1]);ax.set_title('C. '+sp+'\n'+counts,fontsize=6.6,pad=10);ax.set_yticks(range(5),labels if ax==aa[0] else ['']*5,fontsize=6.5)
# Preserve the original exact-P and BH-FDR columns.
t=f.add_axes([.837,.55,.15,.33]);t.set_ylim(4.5,-.5);t.set_xlim(0,1);t.axis('off');t.text(.16,1.07,'Exact P',transform=t.transAxes,ha='center',fontsize=6.4);t.text(.78,1.07,'BH FDR',transform=t.transAxes,ha='center',fontsize=6.4)
for i,(p,q) in enumerate(zip([f'{d.loc[k,"exact_two_sided_p"]:.4f}' if pd.notna(d.loc[k,"exact_two_sided_p"]) else 'Not Tested' for k in ids],[f'{d.loc[k,"fdr_bh"]:.2f}' if pd.notna(d.loc[k,"fdr_bh"]) else '—' for k in ids])):t.text(.16,i,p,ha='center',va='center',fontsize=6.2);t.text(.78,i,q,ha='center',va='center',fontsize=6.2)
f.text(.30,.492,'AX, Axilla; OS, Other Skin',fontsize=6.6,color=gray)
f.text(.045,.435,'B  Genomic Neighborhood Context',fontsize=8.5,fontweight='bold');ax=f.add_axes([.30,.18,.67,.23]);m=read('candidate_mobile_context_summary.tsv');m=m[m.analysis_arm=='cory'].set_index(['target_pgfam_id','group']);mid=[ids[0],ids[2],ids[3]]
dumb(ax,[m.loc[(i,'AXILLA'),'fraction_mobile_context'] for i in mid],[m.loc[(i,'OTHER_SKIN'),'fraction_mobile_context'] for i in mid]);axis(ax,3);ax.set_yticks(range(3),[labels[i] for i in [0,2,3]],fontsize=6.5);ax.set_xlabel('Fraction of Neighborhoods With a Mobile/Phage Annotation',fontsize=7)
f.legend(handles=key(),loc='lower center',bbox_to_anchor=(.59,.085),ncol=2,frameon=False,fontsize=7)
save(f,3,'Figure 3. Exploratory Corynebacterium Loci and Genomic Context')
# Figure 4 preserves the original 2-row/3-column arrangement and site-wise summaries.
d=read('malodor_best_hit_per_genome_6.tsv');d['detected']=d.detected_hit.astype(str).str.lower().isin(['true','1'])
f,axs=plt.subplots(2,3,figsize=(6.5,6.5),sharex=True,sharey=True);f.subplots_adjust(left=.255,right=.985,top=.95,bottom=.19,wspace=.13,hspace=.20)
for ax,ref,letter in zip(axs.flat,REFS,'ABCDEF'):
 for i,sp in enumerate(SPEC):
  ss=d[(d.reference==ref)&(d.species==sp)].sort_values('genome_id')
  for case in [True,False]:
   group=ss[ss.source_group.isin(['AXILLA','AXILLARY_GROIN'])==case];y0=i+(-.13 if case else .13)
   vals=group[group.detected&(group.query_coverage_pct>=70)].pident.to_numpy()
   if len(vals):
    med=np.median(vals);ax.plot([med,med],[y0-.12,y0+.12],lw=1.15,color='black',zorder=10)
   for j,(_,r) in enumerate(group.iterrows()):
    y=y0+((j%5)-2)*.025
    if not r.detected:ax.plot(-15,y,'s',mfc='white',mec='#999',ms=2.2,mew=.4)
    elif r.query_coverage_pct<70:ax.plot(r.pident,y,'x',color='#888',ms=2.5,mew=.5)
    else:ax.plot(r.pident,y,'o',mfc=ORANGE if case else 'white',mec=ORANGE if case else '#777',ms=2.5,mew=.4,alpha=.8)
 ax.set_title(letter+'  '+ref,loc='left',weight='bold',fontsize=8.5);ax.set_xlim(-23,105);ax.set_ylim(9.5,-.5);ax.set_xticks([-15,0,50,100],['ND','0','50','100']);ax.tick_params(axis='x',labelsize=6);ax.axhline(2.5,color='#ddd',lw=.6)
 for y in range(10):ax.axhline(y,color='#f0f0f0',lw=.4,zorder=0)
axs[0,0].set_yticks(range(10),[abbr(sp)+f' (n={len(d[(d.reference==REFS[0])&(d.species==sp)])})' for sp in SPEC],fontsize=6.3)
for ax in axs[1]:ax.set_xlabel('Amino-Acid Identity (%)',fontsize=6.6)
handles=[Line2D([],[],marker='o',ls='',color=ORANGE,label='Axilla / Axillary-Groin',ms=3),Line2D([],[],marker='o',ls='',mfc='white',mec='#777',label='Other Skin',ms=3),Line2D([],[],marker='x',ls='',color='#888',label='Query Coverage <70%',ms=3),Line2D([],[],marker='s',ls='',mfc='white',mec='#999',label='ND: No Reported Hit',ms=3),Line2D([],[],marker='|',ls='',color='black',label='Median',ms=6)]
f.legend(handles=handles,loc='lower center',bbox_to_anchor=(.52,.085),ncol=3,frameon=False,fontsize=6.3)
save(f,4,'Best-Hit Similarity to Six Malodor-Related Reference Proteins')



def key(case='Axilla',other='Other Skin'):
 return [Line2D([],[],marker='o',ls='',color=ORANGE,label=case,markersize=4),Line2D([],[],marker='o',ls='',color=GRAY,markerfacecolor='white',label=other,markersize=4)]
for n,genus,stem,h in [(5,'Cory','cory',6.6),(6,'Staph','staph',8.8)]:
 meta=pd.read_csv(SRC/f'Supplementary{genus}.tsv',sep='\t',dtype={'genome_id':str}).sort_values('tip_order').reset_index(drop=True)
 tree=Phylo.read(SRC/f'{stem}_malodor_codon_tree_tree.nwk','newick');N=len(meta);lookup=dict(zip(meta.genome_id,range(N)))
 tip_names={t.name for t in tree.get_terminals()};assert tip_names==set(lookup),(n,tip_names-set(lookup))
 depths=tree.depths();ys={}
 def y(cl):
  ys[cl]=lookup[cl.name] if cl.is_terminal() else sum(y(c) for c in cl.clades)/len(cl.clades)
  return ys[cl]
 y(tree.root)
 f=plt.figure(figsize=(6.5,h));ax=f.add_axes([.04,.22,.35,.70]);lab=f.add_axes([.40,.22,.24,.70]);site=f.add_axes([.65,.22,.025,.70]);heat=f.add_axes([.70,.22,.27,.70])
 for c in tree.find_clades():
  if c.clades:
   ax.plot([depths[c],depths[c]],[min(ys[x] for x in c.clades),max(ys[x] for x in c.clades)],color='#333',lw=.5)
   for ch in c.clades:ax.plot([depths[c],depths[ch]],[ys[ch],ys[ch]],color='#333',lw=.5)
 ax.set_yticks([]);ax.set_xlabel('Branch Length');ax.set_xlim(0,max(depths.values())*1.03);ax.set_title('A  Codon Tree',loc='left',weight='bold')
 lab.axis('off');site.axis('off');lab.set_xlim(0,1);site.set_xlim(0,1)
 colors={'AXILLA':ORANGE,'OTHER_SKIN':'#d6d6d6','AXILLARY_GROIN':'#555555'}
 for i,r in meta.iterrows():site.add_patch(Rectangle((0,i-.5),1,1,color=colors[r.source_group],lw=0))
 start=0
 for i in range(1,N+1):
  if i==N or meta.species.iloc[i]!=meta.species.iloc[start]:
   lab.text(.02,(start+i-1)/2,abbr(meta.species.iloc[start]),va='center',fontsize=6.7)
   for a in [ax,lab,heat]:a.axhline(i-.5,color='#ddd',lw=.5)
   start=i
 arr=np.column_stack([pd.to_numeric(meta[r+'_identity_pct'],errors='coerce') for r in REFS]);cmap=plt.get_cmap('Greys').copy();cmap.set_bad('white')
 im=heat.imshow(arr,aspect='auto',vmin=0,vmax=100,cmap=cmap,interpolation='nearest')
 for i in range(N):
  for j,ref in enumerate(REFS):
   if np.isnan(arr[i,j]):heat.text(j,i,'×',ha='center',va='center',fontsize=5 if N==100 else 6,color='#555')
   elif float(meta.loc[i,ref+'_query_coverage_pct'])<70:heat.add_patch(Rectangle((j-.5,i-.5),1,1,fill=False,hatch='////',edgecolor='#b66b30',lw=.2))
 heat.set_xticks(range(6),REFS,rotation=90);heat.set_yticks([]);heat.set_title('B  Similarity',loc='left',weight='bold')
 for a in [ax,lab,site,heat]:a.set_ylim(N-.5,-.5)
 cb=f.add_axes([.72,.13,.23,.016]);f.colorbar(im,cax=cb,orientation='horizontal',ticks=[0,50,100]);cb.set_xlabel('Amino-Acid Identity (%)',fontsize=7,labelpad=2)
 handles=[Patch(fc=col,label=gr.replace('_',' ').title()) for gr,col in colors.items() if gr in set(meta.source_group)]
 handles += [Line2D([],[],marker='x',ls='',color='#555',label='No Reported Hit',ms=4),Patch(fc='white',ec='#b66b30',hatch='////',label='Query Coverage <70%')]
 f.legend(handles=handles,loc='lower left',bbox_to_anchor=(.025,.06),ncol=2,frameon=False,fontsize=6.5)
 save(f,n,('Corynebacterium' if n==5 else 'Staphylococcus')+' Phylogeny and Similarity')

p=read('Figure7_pairwise_exact_CDS_overlap_v3.tsv').sort_values('pair_order')
s=read('Figure7B_PatB_AecD_shared_species_summary_v4.tsv').sort_values('order')
f=plt.figure(figsize=(6.5,8.1));a=f.add_axes([.30,.56,.58,.36]);b=f.add_axes([.30,.16,.58,.32])
a.barh(range(15),p.pct_same_exact_CDS,color=[ORANGE if x else '#aaa' for x in p.same_exact_CDS],height=.55);a.plot(p.pct_same_exact_CDS,range(15),'o',color='#555',ms=2)
a.set_yticks(range(15),p.pair,fontsize=7);a.set_ylim(14.5,-.5);a.set_xlim(0,105);a.set_xlabel('Co-Detected Genomes With the Same CDS (%)');a.set_title('A  Exact-CDS Overlap',loc='left',weight='bold')
for i,r in p.iterrows():a.text(107,r.pair_order,f'{r.same_exact_CDS}/{r.both_detected}',va='center',fontsize=6.5)
a.text(107,-1.2,'Same / Both',fontsize=6.3)
for i,r in s.reset_index(drop=True).iterrows():
 b.plot([r.PatB_median,r.AecD_median],[i,i],color='#bbb',lw=1);b.plot(r.PatB_median,i,'o',mfc='white',mec=GRAY,ms=4);b.plot(r.AecD_median,i,'s',color=ORANGE,ms=4)
b.set_yticks(range(10),[abbr(r.species)+f' (n={r.n})' for _,r in s.iterrows()],fontsize=6.7);b.set_ylim(9.5,-.5);b.set_xlim(0,105);b.set_xlabel('Median Amino-Acid Identity (%)');b.set_title('B  Similarity of the Shared CDS',loc='left',weight='bold')
f.legend(handles=[Line2D([],[],marker='o',mfc='white',mec=GRAY,ls='',label='PatB Median'),Line2D([],[],marker='s',color=ORANGE,ls='',label='AecD Median')],loc='lower center',bbox_to_anchor=(.58,.055),ncol=2,frameon=False)
save(f,7,'Exact-CDS Overlap and Reference Similarity')

import shutil
shutil.rmtree(TMP)
print('Saved Figures 1–7 as PDFs and 400-dpi JPEGs to',OUT)
