"""Restyle original four-band review Figure 04 as manuscript Figure 07b."""
from pathlib import Path
import sys,json,hashlib,shutil,platform
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.functions.plotting import PALETTE
RUN='AGE-HR-FOUR-STYLE-20261017-01'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    audit=ROOT/'docs/age_interaction_runs'/RUN
    audit.mkdir(parents=True,exist_ok=True)
    source=ROOT/'tables/62_review_1008_age_specific_HR.csv'
    curve_source=ROOT/'tables/62_review_1008_frailty_refitted_curve.csv'
    inputs={str(p.relative_to(ROOT)):sha(p) for p in [source,curve_source]}
    state=dict(run_id=RUN,status='running',started_utc=datetime.now(timezone.utc).isoformat(),command='python scripts/plot_figure07b_age_HR.py',code_sha256=sha(__file__),input_sha256=inputs,models_refitted=False,participant_data_accessed=False,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
    (audit/'provenance.json').write_text(json.dumps(state,indent=2)+'\n')
    frame=pd.read_csv(source)
    old_labels=['<45','45–54','55–64','≥65']
    labels=['<45','45–<55','55–<65','≥65']
    terms=['BRCA1/2','Any carrier','Other carriers','PRS311']
    legend=['BRCA1 or BRCA2 carriers','Any panel-gene carrier','Other panel-gene carriers','PRS311, per SD']
    colors=[PALETTE[k] for k in ['red','blue1','blue2','orange']]
    curve=pd.read_csv(curve_source)
    curve=curve.loc[np.isclose(curve.theta,.73)].sort_values('age')
    assert len(curve)>4 and not curve.age.duplicated().any() and curve.HR.gt(0).all()
    ages=[35,50,60,72]
    expected=np.interp(ages,curve.age,curve.HR)
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','axes.grid':False})
    fig,ax=plt.subplots(figsize=(8,6.7))
    x=np.arange(4)
    ax.plot(x,expected,linestyle='--',color=PALETTE['dark_grey'],lw=1.8,
        label='BRCA1/2 expected from\ndepletion alone (frailty 0.73)')
    plotted=[]
    for term,label,color,offset in zip(terms,legend,colors,[-.18,-.06,.06,.18]):
        sub=frame.loc[frame.term.eq(term)].set_index('band').loc[old_labels].copy()
        assert len(sub)==4 and np.isfinite(sub[['estimate','CI_low','CI_high']]).all().all()
        assert ((sub.CI_low>0)&(sub.CI_low<=sub.estimate)&(sub.estimate<=sub.CI_high)).all()
        np.testing.assert_allclose(sub.estimate,np.exp(sub.beta))
        np.testing.assert_allclose(sub.CI_low,np.exp(sub.beta-1.959963984540054*sub.se))
        np.testing.assert_allclose(sub.CI_high,np.exp(sub.beta+1.959963984540054*sub.se))
        ax.errorbar(x+offset,sub.estimate,yerr=np.array([sub.estimate-sub.CI_low,sub.CI_high-sub.estimate]),
            fmt='s',linestyle='--',color=color,lw=1.3,elinewidth=1.3,capsize=4,markersize=6,label=label)
        sub=sub.reset_index();sub['age_group']=labels;sub['color']=color
        plotted.append(sub)
    ax.set_xticks(x,labels)
    ax.set_xlabel('Attained age (years)',labelpad=12)
    ax.set_ylabel('Hazard ratio (95% CI)')
    ax.set_yscale('log');ax.set_xlim(-.5,3.5);ax.set_ylim(1,27)
    ax.yaxis.set_major_locator(FixedLocator([1,2,3,5,10,20]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(0,1.025),fontsize=10,ncol=2)
    fig.text(.5,.025,'Squares: HR estimates; error bars: 95% CIs. Dashed lines connect age-specific estimates.\nGray dashed curve: existing illustrative depletion scenario with frailty 0.73.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.1,1,.86])
    target=ROOT/'figures/Figure_07b_age_HR.svg'
    archive=ROOT/'manuscript/archive/2026-10-17_figure07b_original_style'
    archive.mkdir(parents=True,exist_ok=True)
    if target.exists() and not (archive/target.name).exists():shutil.copy2(target,archive/target.name)
    fig.savefig(target,facecolor='white')
    fig.savefig('/tmp/Figure_07b_preview.png',dpi=150,facecolor='white');plt.close(fig)
    table=ROOT/'tables/Figure_07b_plot_data.csv'
    pd.concat(plotted,ignore_index=True).to_csv(table,index=False)
    curveout=ROOT/'tables/Figure_07b_frailty_display.csv'
    pd.DataFrame(dict(age_group=labels,representative_age=ages,HR=expected,theta=.73)).to_csv(curveout,index=False)
    assert all(sha(ROOT/p)==digest for p,digest in inputs.items())
    state.update(status='completed',ended_utc=datetime.now(timezone.utc).isoformat(),age_bands=labels,colors=dict(zip(legend,colors)),gray_curve='Existing baseline-refitted theta=.73 illustrative scenario; original representative ages35/50/60/72 retained',verification='16 HR/95%CI points checked against exp(beta±1.95996398454SE); exact source age bands and hashes preserved',output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [target,table,curveout]},scientific_baseline_changed=False,limitations=['No boxplot distribution available: square point estimates and95%CIs','Original model-based Breslow Cox uncertainty retained','Frailty curve is an existing independently refitted scenario, not a newly identified theta estimate'])
    (audit/'provenance.json').write_text(json.dumps(state,indent=2)+'\n')
    print('Completed Figure07b: original four age bands; sixteen verified HR/CI estimates; five legend groups.')
if __name__=='__main__':main()
