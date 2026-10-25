"""Run archived analyses in an isolated workspace populated from bundled data."""
from pathlib import Path
import argparse, gzip, shutil, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('analysis',choices=['blast','sensitivity','strict','same-study','cory','haemolyticus','ecoli']);p.add_argument('--workdir',type=Path,required=True);a=p.parse_args()
w=a.workdir.resolve()
if w.exists() and any(w.iterdir()):raise SystemExit('Choose an empty work directory to avoid overwriting results.')
w.mkdir(parents=True,exist_ok=True)
for src in (ROOT/'data/archive').iterdir():
 if src.suffix=='.gz':
  with gzip.open(src,'rb') as i,open(w/src.stem,'wb') as o:shutil.copyfileobj(i,o)
 else:shutil.copyfile(src,w/src.name)
scripts={'blast':'analyze_targeted_malodor_blast_6.py','strict':'staph_epi_strict_primary_pgfam_analysis.py','same-study':'staph_epi_prjna803478_same_study_analysis.py','cory':'cory_exact_site_stratified_analysis.py','haemolyticus':'staph_haem_exploratory_analysis.py','ecoli':'ecoli_exploratory_analysis.py'}
if a.analysis=='sensitivity':
 cmd=[sys.executable,str(ROOT/'scripts/run_control_sensitivity.py'),'--archive',str(w),'--added-features',str(w/'staph_epi_four_missing_features.tsv'),'--output',str(w/'sensitivity_results')]
else:cmd=[sys.executable,str(ROOT/'scripts/archived'/scripts[a.analysis])]
subprocess.run(cmd,cwd=w,check=True)
print('Analysis workspace:',w)
