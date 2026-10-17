"""Reproduce Figure 04a from frozen case inputs; only aggregate outputs."""
from pathlib import Path
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
RUN = 'AGE-PRS-FOUR-GROUP-C-20261017-01'
sys.path.insert(0, str(ROOT))
from scripts.functions.plotting import PALETTE

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def main():
    audit_dir = ROOT / 'docs/age_interaction_runs' / RUN
    audit_dir.mkdir(parents=True, exist_ok=True)
    tables = ROOT / 'tables'
    tables.mkdir(parents=True, exist_ok=True)
    audit = audit_dir / 'Figure_04a_provenance.json'
    state = dict(run_id=RUN, status='running', started_utc=now(),
                 command='python scripts/Figure_04a_plot.py',
                 code_sha256=sha(__file__),
                 software=dict(python=platform.python_version(), numpy=np.__version__,
                               pandas=pd.__version__, matplotlib=matplotlib.__version__),
                 definition='Shared PRS311 score percentiles P25/P50/P75 across the full frozen case-control cohort (141204 participants); score boundaries right-inclusive; age bins left-closed',
                 figure04c_definition='Four carrier/PRS combinations per panel; median low/high and P20/P80 extremes; full-cohort reference',
                 participant_rows_written=False, models_refitted=False, supersedes='AGE-PRS-POPULATION-QUARTILES-20261017-01 Figure04c layout only; Figure04a unchanged')
    audit.write_text(json.dumps(state, indent=2)+'\n')
    try:
        lock = json.loads((ROOT/'data/analysis_runs/AGE-INT-20261016-01/input_lock.json').read_text())['files']
        keys = ['frozen_cohort', 'age_ids', 'scores']
        for key in keys:
            assert sha(ROOT/lock[key]['path']) == lock[key]['sha256']
        state['input_sha256'] = {k:lock[k]['sha256'] for k in keys}
        f = pd.read_csv(ROOT/lock['frozen_cohort']['path'], usecols=['person_id','is_case','age_at_diagnosis','PRS311'])
        ids = pd.read_csv(ROOT/lock['age_ids']['path']).person_id
        d = f.loc[f.person_id.isin(ids)].copy()
        assert len(d)==7996 and d.person_id.is_unique and d.is_case.eq(1).all()
        assert np.isfinite(d[['age_at_diagnosis','PRS311']]).all().all()
        labels = ['<35','35–<45','45–<50','50–<55','55–<60','60–<65','65–<70','≥70']
        totals = np.array([155,737,747,1004,1092,1240,1358,1663])
        assert len(f) == 141204 and f.person_id.is_unique
        assert np.isfinite(f.PRS311).all()
        thresholds = f.PRS311.quantile([.25,.5,.75]).to_numpy()
        state['reference_population_n'] = len(f)
        state['reference_population'] = 'Full frozen case-control cohort; cases and controls combined'
        np.testing.assert_allclose(thresholds, np.percentile(f.PRS311,[25,50,75]))
        levels = ['0–25%', '25–50%', '50–75%', '75–100%']
        d['quartile'] = pd.cut(d.PRS311, [-np.inf,*thresholds,np.inf], labels=levels)
        d['age_bin'] = pd.cut(d.age_at_diagnosis, [-np.inf,35,45,50,55,60,65,70,np.inf], right=False, labels=labels)
        counts = pd.crosstab(d.age_bin,d.quartile).reindex(index=labels,columns=levels,fill_value=0)
        np.testing.assert_array_equal(counts.sum(axis=1),totals)
        independent = np.histogram2d(d.age_at_diagnosis, d.PRS311,
            bins=[[-np.inf,35,45,50,55,60,65,70,np.inf],[-np.inf,*np.nextafter(thresholds,np.inf),np.inf]])[0]
        np.testing.assert_array_equal(counts.to_numpy(),independent)
        release = counts.ge(20).all(axis=1)
        plot_counts = counts.where(release, other=0, axis=0)
        scores = pd.read_csv(ROOT/lock['scores']['path'], usecols=['sample_id','carried_gene_count_C']).set_index('sample_id')
        assert scores.index.is_unique
        d['carrier'] = scores.loc[d.person_id,'carried_gene_count_C'].to_numpy() > 0
        age_rows = []
        for carrier in [False, True]:
            for level in levels:
                ages = d.loc[d.carrier.eq(carrier) & d.quartile.eq(level), 'age_at_diagnosis']
                assert len(ages) >= 20, 'Small age distribution group: withhold output'
                q = ages.quantile([.1,.25,.5,.75,.9])
                age_rows.append(dict(carrier=carrier, PRS_score_percentile=level, n=len(ages),
                    age_P10=q.loc[.1], age_Q1=q.loc[.25], age_median=q.loc[.5],
                    age_Q3=q.loc[.75], age_P90=q.loc[.9]))
        age_table = pd.DataFrame(age_rows)
        assert age_table.n.sum() == 7996
        plt.rcParams.update({'font.size':11,'svg.fonttype':'none','axes.grid':False})
        fig,ax = plt.subplots(figsize=(11,5.8))
        x=np.arange(8); bottom=np.zeros(8)
        colors=[PALETTE[k] for k in ['red','blue1','orange','dark_grey']]
        for q,color,legend in zip(counts.columns,colors,levels):
            values=plot_counts[q].to_numpy()
            ax.bar(x,values,bottom=bottom,width=.72,color=color,edgecolor='white',linewidth=.6,label=legend)
            for pos,value,base in zip(x,values,bottom):
                if value == 0:
                    continue
                ax.text(pos,base+value/2,f'{value:,}',ha='center',va='center',fontsize=6 if pos == 0 else 9,color='white' if color == PALETTE['dark_grey'] else '#222222')
            bottom+=values
        for pos,n,safe in zip(x,totals,release):
            if not safe:
                ax.bar(pos,n,width=.72,color='#D9D9D9',edgecolor='white',linewidth=.6)
                ax.text(pos,n+15,'PRS split\nwithheld',ha='center',va='bottom',fontsize=8)
        ax.set_xticks(x,[f'{label}\n(n={n:,})' for label,n in zip(labels,totals)])
        ax.set_xlabel('Age at diagnosis (years)',labelpad=13)
        ax.set_ylabel('Count')
        ax.set_ylim(0,1800); ax.set_xlim(-.65,7.65)
        ax.spines[['top','right']].set_visible(False)
        ax.legend(title='PRS311 score percentile (low → high)',ncol=2,frameon=False,loc='upper left')
        fig.text(.5,.025,('PRS score percentiles: full case–control cohort (n=141,204); bars: 7,996 cases.' + ('\nPRS breakdown withheld for age groups containing a cell <20.' if not release.all() else '')),ha='center',fontsize=10)
        fig.tight_layout(rect=[0,.07,1,1])
        prefix=OUT/'Figure_04a_PRS_quartile_counts_by_age'
        fig.savefig(prefix.with_suffix('.svg'),facecolor='white')
        fig.savefig('/tmp/Figure_04a_preview.png',dpi=150,facecolor='white')
        plt.close(fig)
        counts.where(release,axis=0).assign(n_cases=totals,PRS_breakdown_released=release).rename_axis('age_group').to_csv(tables/'Figure_04a_plot_data.csv')
        # Figure 04c: four carrier/PRS combinations in each of two panels.
        c_thresholds = f.PRS311.quantile([.2,.5,.8])
        c_rows = []
        fig,axes = plt.subplots(1,2,figsize=(12,5.6),sharey=True)
        for ax,kind,title,group_labels in zip(
            axes,['median','extremes'],
            ['PRS split at the cohort median','PRS extremes: bottom 20% versus top 20%'],
            [['Low','High'],['Bottom 20%','Top 20%']]):
            tick_labels = []
            for i,(carrier,high) in enumerate([(False,False),(False,True),(True,False),(True,True)]):
                if kind == 'median':
                    mask = d.PRS311.ge(c_thresholds.loc[.5]) if high else d.PRS311.lt(c_thresholds.loc[.5])
                else:
                    mask = d.PRS311.ge(c_thresholds.loc[.8]) if high else d.PRS311.le(c_thresholds.loc[.2])
                ages = d.loc[d.carrier.eq(carrier) & mask,'age_at_diagnosis']
                assert len(ages) >= 20, 'Small four-group age distribution: withhold output'
                q = ages.quantile([.1,.25,.5,.75,.9])
                color = colors[i]
                ax.bxp([dict(med=q.loc[.5],q1=q.loc[.25],q3=q.loc[.75],whislo=q.loc[.1],whishi=q.loc[.9])],
                    positions=[i+1],widths=.48,patch_artist=True,showfliers=False,manage_ticks=False,
                    boxprops=dict(facecolor=color,edgecolor=color,linewidth=1.4),
                    whiskerprops=dict(color=color,linewidth=1.5),capprops=dict(color=color,linewidth=1.5),
                    medianprops=dict(color='white',linewidth=2.2))
                carrier_label = 'Carrier' if carrier else 'Non-carrier'
                prs_label = group_labels[int(high)]
                tick_labels.append(f'{carrier_label}\nPRS {prs_label}\n(n={len(ages):,})')
                c_rows.append(dict(panel=kind,carrier=carrier,PRS_group=prs_label,n=len(ages),
                    age_P10=q.loc[.1],age_Q1=q.loc[.25],age_median=q.loc[.5],age_Q3=q.loc[.75],age_P90=q.loc[.9]))
            ax.set_xticks(range(1,5),tick_labels,fontsize=9)
            ax.set_xlim(.4,4.6); ax.set_ylim(30,85)
            ax.set_title(title,loc='left',fontsize=11,pad=14)
            ax.set_xlabel('Carrier status and PRS group',labelpad=12)
            ax.spines[['top','right']].set_visible(False)
        axes[0].set_ylabel('Recorded diagnosis age (years)')
        fig.text(.5,.025,'PRS thresholds: full case–control cohort (n=141,204). Right panel excludes middle 60%.\nWhite line: median; box: Q1–Q3; whiskers: P10–P90.',ha='center',fontsize=10)
        fig.tight_layout(rect=[0,.11,1,1],w_pad=3)
        fig.savefig(OUT/'Figure_04c_age_median_and_extreme_PRS.svg',facecolor='white')
        fig.savefig('/tmp/Figure_04c_preview.png',dpi=150,facecolor='white')
        plt.close(fig)
        c_table = pd.DataFrame(c_rows)
        assert c_table.loc[c_table.panel.eq('median'),'n'].sum() == 7996
        assert c_table.loc[c_table.panel.eq('extremes'),'n'].sum() == (d.PRS311.le(c_thresholds.loc[.2]) | d.PRS311.ge(c_thresholds.loc[.8])).sum()
        c_table.to_csv(tables/'Figure_04c_four_group_age_summary.csv',index=False)
        pd.DataFrame({'percentile':[20,50,80],'PRS311_score_threshold':c_thresholds.to_numpy(),'reference_n':len(f)}).to_csv(tables/'Figure_04c_PRS_group_thresholds.csv',index=False)
        pd.DataFrame({'percentile':[25,50,75],'PRS311_score_threshold':thresholds,'reference_n':len(f)}).to_csv(tables/'Figure_04_PRS_quartile_thresholds.csv',index=False)
        for key in keys:
            assert sha(ROOT/lock[key]['path'])==lock[key]['sha256']
        state.update(status='completed', n_cases=7996, all_released_cells_ge20=True, suppressed_age_groups=release.index[~release].tolist(),
                     verification='Exact eight age totals; independent histogram agrees; unchanged frozen input hashes',
                     output_sha256={p.name:sha(p) for p in [prefix.with_suffix('.svg'),tables/'Figure_04a_plot_data.csv',OUT/'Figure_04c_age_median_and_extreme_PRS.svg',tables/'Figure_04c_four_group_age_summary.csv',tables/'Figure_04c_PRS_group_thresholds.csv',tables/'Figure_04_PRS_quartile_thresholds.csv']},
                     scientific_baseline_changed=False)
    except Exception as exc:
        state.update(status='failed',failure_reason=type(exc).__name__)
        raise
    finally:
        state['ended_utc']=now()
        audit.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({'status':state['status'],'n_cases':7996,'all_released_cells_ge20':True}))

if __name__=='__main__':
    main()
