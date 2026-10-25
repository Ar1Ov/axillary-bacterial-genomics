"""Recompute strict PGFam tests; optionally add four controls from a complete export.
Usage: python run_control_sensitivity.py --archive PATH --output PATH
       [--added-features complete_four_controls.tsv]
Dependencies: pandas, numpy, scipy. Genome identifiers are always read as strings.
An added export must contain genome_id, feature.patric_id, feature.pgfam_id.
Supply complete feature exports, not a list of selected candidate families.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

ADDED = {'1282.4413', '1282.4414', '1282.4415', '1282.4416'}

def read(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)

def calculate(cohort, features):
    pairs = features.loc[features['feature.pgfam_id'].ne(''), ['genome_id','feature.pgfam_id']].drop_duplicates()
    assert set(pairs.genome_id) <= set(cohort.genome_id)
    ax = set(cohort.loc[cohort['group'].eq('AXILLA'), 'genome_id'])
    ct = set(cohort.loc[cohort['group'].eq('OTHER_SKIN'), 'genome_id'])
    rows = []
    for family, g in pairs.groupby('feature.pgfam_id'):
        ids = set(g.genome_id); a = len(ids & ax); c = len(ids & ct)
        tested = a+c >= 3 and len(ax)+len(ct)-a-c >= 3
        p = fisher_exact([[a,len(ax)-a],[c,len(ct)-c]], alternative='two-sided').pvalue if tested else np.nan
        rows.append(dict(pgfam_id=family, AXILLA_present=a, OTHER_SKIN_present=c,
                         prevalence_difference=a/len(ax)-c/len(ct), tested=tested, p_value=p))
    result = pd.DataFrame(rows)
    result['fdr_bh'] = np.nan
    mask = result.tested
    result.loc[mask,'fdr_bh'] = false_discovery_control(result.loc[mask,'p_value'].to_numpy(),method='bh')
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--added-features', type=Path)
    args = parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    stem = args.archive/'staph_epi_strict_primary'
    cohort = read(str(stem)+'_cohort.tsv'); features = read(str(stem)+'_features.tsv')
    assert cohort.genome_id.is_unique
    assert cohort['group'].value_counts().to_dict() == {'OTHER_SKIN':117,'AXILLA':5}
    baseline = calculate(cohort,features)
    saved = read(str(stem)+'_pgfam_results.tsv').set_index('pgfam_id')
    observed = baseline.set_index('pgfam_id').loc[saved.index]
    for col in ['AXILLA_present','OTHER_SKIN_present','prevalence_difference','p_value','fdr_bh']:
        np.testing.assert_allclose(observed[col],pd.to_numeric(saved[col].replace('',np.nan)),rtol=1e-10,atol=1e-12,equal_nan=True)
    baseline.to_csv(args.output/'baseline_recomputed.tsv',sep='\t',index=False)
    def summary(r):
        return {'observed_families':len(r),'tested_families':int(r.tested.sum()),
                'significant_axilla_enriched':int(((r.fdr_bh<.05)&(r.prevalence_difference>0)).sum()),
                'significant_axilla_depleted':int(((r.fdr_bh<.05)&(r.prevalence_difference<0)).sum())}
    report = {'baseline_matches_saved':True,'baseline':summary(baseline),'expanded_run_completed':False}
    if args.added_features:
        added = read(args.added_features)
        assert {'genome_id','feature.patric_id','feature.pgfam_id'} <= set(added.columns)
        assert set(added.genome_id)==ADDED, 'Export must contain exactly the four omitted genomes'
        assert added['feature.patric_id'].ne('').all() and added['feature.patric_id'].is_unique
        assert all(added.loc[added.genome_id.eq(g),'feature.pgfam_id'].ne('').any() for g in ADDED)
        cohort2=pd.concat([cohort,pd.DataFrame({'genome_id':sorted(ADDED),'group':'OTHER_SKIN','bioproject':'PRJNA692334'})],ignore_index=True)
        expanded=calculate(cohort2,pd.concat([features,added],ignore_index=True))
        expanded.to_csv(args.output/'expanded_5_vs_121.tsv',sep='\t',index=False)
        cohort2.to_csv(args.output/'expanded_cohort.tsv',sep='\t',index=False)
        comparison=baseline.merge(expanded,on='pgfam_id',how='outer',suffixes=('_baseline','_expanded'))
        comparison.to_csv(args.output/'all_family_comparison.tsv',sep='\t',index=False)
        report.update(expanded_run_completed=True,expanded=summary(expanded))
    (args.output/'run_status.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
