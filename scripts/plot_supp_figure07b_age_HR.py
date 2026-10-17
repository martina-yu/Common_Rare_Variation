"""Eight attained-age bands for review Figure 04 / manuscript Figure 07b."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,hashlib,platform,shutil,warnings
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import scipy,statsmodels
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.review.review_1008 import prepare_survival,split_intervals
from scripts.functions.plotting import PALETTE
RUN='AGE-HR-SEVEN-SUPP-20261017-01'
LABELS=['<40','40–<45','45–<50','50–<55','55–<60','60–<65','≥65']
EDGES=[18,40,45,50,55,60,65,np.inf]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()

def draw(table):
    paths={'curve':ROOT/'tables/62_review_1008_frailty_refitted_curve.csv'}
    curve=pd.read_csv(paths['curve']);curve=curve.loc[np.isclose(curve.theta,.73)].sort_values('age')
    assert len(curve)>8 and curve.HR.gt(0).all()
    # Reuse existing gamma .73 scenario, evaluated at fixed display ages.
    representative=np.array([35,42.5,47.5,52.5,57.5,62.5,72])
    expected=np.interp(representative,curve.age,curve.HR)
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','axes.grid':False})
    fig,ax=plt.subplots(figsize=(12,6.8))
    x=np.arange(7)
    gray=PALETTE['dark_grey']
    ax.plot(x,expected,'--',color=gray,lw=1.8,label='BRCA1/2 expected from\ndepletion alone (frailty 0.73)',zorder=2)
    offsets=[-.18,-.06,.06,.18]
    for term,label,color,offset in zip(['BRCA1/2','Any carrier','Other carriers','PRS311'],['BRCA1 or BRCA2 carriers','Any panel-gene carrier','Other panel-gene carriers','PRS311, per SD'],[PALETTE[k] for k in ['red','blue1','blue2','orange']],offsets):
        sub=table.loc[table.term.eq(term)].set_index('age_group').loc[LABELS]
        values=sub.estimate.to_numpy(float);lo=sub.CI_low.to_numpy(float);hi=sub.CI_high.to_numpy(float)
        ax.errorbar(x+offset,values,yerr=np.array([values-lo,hi-values]),fmt='s',linestyle='--',color=color,elinewidth=1.3,lw=1.2,capsize=3,markersize=5,label=label,zorder=3)
    ax.set_xticks(x,LABELS)
    ax.set_xlabel('Attained age (years)',labelpad=12)
    ax.set_ylabel('Hazard ratio (95% CI)')
    ax.set_yscale('log');ax.set_xlim(-.55,6.55)
    ax.set_ylim(1,40)
    ax.yaxis.set_major_locator(FixedLocator([1,2,3,5,10,20,30]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(0,1.025),fontsize=10,ncol=2)
    fig.text(.5,.025,'Squares: HR estimates; error bars: 95% CIs. Dashed lines connect age-specific estimates.\nGray curve: existing illustrative frailty scenario.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.1,1,1])
    target=ROOT/'figures/supp_Figure_07b_seven_age_HR.svg'
    archive=ROOT/'manuscript/archive/2026-10-17_figure07b_four_age_bands'
    archive.mkdir(parents=True,exist_ok=True)
    if target.exists():shutil.copy2(target,archive/target.name)
    fig.savefig(target,facecolor='white')
    fig.savefig('/tmp/supp_Figure_07b_preview.png',dpi=150,facecolor='white');plt.close(fig)
    outtable=ROOT/'tables/supp_Figure_07b_seven_age_HR.csv'
    table.to_csv(outtable,index=False)
    curveout=ROOT/'tables/supp_Figure_07b_frailty_display.csv'
    pd.DataFrame(dict(age_group=LABELS,representative_age=representative,expected_HR=expected,theta=.73)).to_csv(curveout,index=False)
    return target,outtable,curveout

def main():
    private=ROOT/'data/analysis_runs'/RUN
    private.mkdir(parents=True,exist_ok=False)
    public=ROOT/'docs/age_interaction_runs'/RUN
    public.mkdir(parents=True,exist_ok=False)
    shutil.copy2(__file__,private/'source_snapshot.py')
    base=ROOT/'data/29_rowlands_extension'
    paths={'frozen':base/'section9/frozen_model_inputs.csv.gz','scores':base/'section8/comparable_scores_case_control.csv.gz','genes':base/'section7/unweighted_gene_carriers.csv.gz','dates':ROOT/'data/62_review_1008/dates.csv.gz','curve':ROOT/'tables/62_review_1008_frailty_refitted_curve.csv','source_model':ROOT/'scripts/review/review_1008.py'}
    hashes={k:sha(p) for k,p in paths.items()}
    state=dict(run_id=RUN,status='running',started_utc=now(),command='python scripts/plot_supp_figure07b_age_HR.py',code_sha256=sha(__file__),optimizer='trust-exact using score and Hessian of same Breslow Cox likelihood; no penalty',input_sha256=hashes,age_bands=LABELS,
        definitions='Original first_ehr_plus1 design; attained age interval-split Cox; PRS SD across full frozen cohort; PC1–10 and technical effects shared; age-varying PRS/carrier effects; Breslow ties/model-based CIs retained',
        software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,matplotlib=matplotlib.__version__),
        output_rule='Withhold HR/CI for class-specific event support below20; no sparse-event counts released',supersedes=None,failure_reason=None,models=[])
    (private/'status.json').write_text(json.dumps(state,indent=2)+'\n')
    (public/'provenance.json').write_text(json.dumps(state,indent=2)+'\n')
    try:
        d=pd.read_csv(paths['frozen']).set_index('person_id')
        s=pd.read_csv(paths['scores'],usecols=['sample_id','B_raw_C']).set_index('sample_id').loc[d.index]
        g=pd.read_csv(paths['genes'],usecols=['sample_id','BRCA1_main','BRCA2_main']).set_index('sample_id').loc[d.index]
        assert len(d)==141204 and d.index.is_unique
        p=((d.PRS311-d.PRS311.mean())/d.PRS311.std(ddof=1)).to_numpy()
        carrier=s.B_raw_C.to_numpy()>0
        brca=(g.BRCA1_main|g.BRCA2_main).to_numpy(bool)
        assert np.all(~brca|carrier)
        technical=pd.concat([pd.get_dummies(d.seq_site_id,drop_first=True),pd.get_dummies(d.seq_sequencer,drop_first=True)],axis=1).to_numpy(float)
        data=(d,None,None,None,None,d.is_case.to_numpy(),None,None,None,None,None,None,None)
        entry,exit,event,valid,_=prepare_survival(data)
        assert not np.isin(exit[valid&event],EDGES[1:-1]).any(), 'An event lies exactly on an age boundary: review interval closure'
        intervals=split_intervals(entry[valid],exit[valid],event[valid],EDGES)
        ix=intervals[:,0].astype(int);band=intervals[:,4].astype(int)
        assert int(intervals[:,3].sum())==int(event[valid].sum())
        np.testing.assert_allclose(np.bincount(ix,weights=intervals[:,2]-intervals[:,1],minlength=valid.sum()),exit[valid]-entry[valid])
        indicator=np.eye(7)[band]
        common=np.column_stack([d[[f'PC{i}' for i in range(1,11)]].to_numpy(float),technical])[valid][ix]
        rows=[]
        for model,variables,names,masks in [
            ('any_carrier',np.column_stack([p,carrier]),['PRS311','Any carrier'],[np.ones(len(d),bool),carrier]),
            ('gene_classes',np.column_stack([p,brca,carrier&~brca]),['PRS311_class_adjusted','BRCA1/2','Other carriers'],[np.ones(len(d),bool),brca,carrier&~brca])]:
            vx=variables[valid][ix]
            x=np.column_stack([indicator*vx[:,j,None] for j in range(vx.shape[1])]+[common])
            keep=np.ptp(x,axis=0)>1e-12
            assert keep[:7*len(names)].all()
            scale=x[:,keep].std(axis=0)
            assert event[valid].sum()/keep.sum()>=10
            ph=PHReg(intervals[:,2],x[:,keep]/scale,status=intervals[:,3].astype(int),entry=np.nextafter(intervals[:,1],np.inf),ties='breslow')
            print('Fitting '+model,flush=True)
            with warnings.catch_warnings(record=True) as caught:
                from scipy.optimize import minimize
                from types import SimpleNamespace
                n_event=float(intervals[:,3].sum())
                optimum=minimize(lambda b:-ph.loglike(b)/n_event,np.zeros(int(keep.sum())),
                    jac=lambda b:-ph.score(b)/n_event,hess=lambda b:-ph.hessian(b)/n_event,
                    method='trust-exact',options={'gtol':1e-9,'maxiter':200})
                information=-ph.hessian(optimum.x)
                covariance_scaled=np.linalg.inv(information)
                assert np.linalg.eigvalsh(information).min()>0
                fit=SimpleNamespace(params=optimum.x,bse=np.sqrt(np.diag(covariance_scaled)),cov_params=lambda:covariance_scaled)
            score=float(np.max(np.abs(ph.score(fit.params))))
            (private/(model+'_fit_diagnostics.json')).write_text(json.dumps(dict(score_max=score,finite_params=bool(np.isfinite(fit.params).all()),finite_SE=bool(np.isfinite(fit.bse).all()),warnings=[str(w.message) for w in caught])))
            assert np.isfinite(fit.params).all() and np.isfinite(fit.bse).all() and score<1e-3
            beta=fit.params/scale;se=fit.bse/scale
            covariance=fit.cov_params()/np.outer(scale,scale)
            np.testing.assert_allclose(se,np.sqrt(np.diag(covariance)))
            state['models'].append(dict(model=model,completed_utc=now(),score_max=score,n_parameters=int(keep.sum()),warnings=[type(w.message).__name__ for w in caught]))
            (private/(model+'_coefficients.json')).write_text(json.dumps(dict(beta=beta.tolist(),covariance=covariance.tolist())))
            for j,(term,mask) in enumerate(zip(names,masks)):
                for k,label in enumerate(LABELS):
                    loc=int(np.flatnonzero(np.flatnonzero(keep)==j*7+k)[0])
                    supported=int(np.sum(event[valid][ix].astype(bool)&intervals[:,3].astype(bool)&(band==k)&mask[valid][ix]))>=20
                    value=float(np.exp(beta[loc]));lo=float(np.exp(beta[loc]-norm.ppf(.975)*se[loc]));hi=float(np.exp(beta[loc]+norm.ppf(.975)*se[loc]))
                    assert 0<lo<=value<=hi and np.isfinite(hi)
                    rows.append(dict(model=model,age_group=label,term=term,estimate=value if supported else None,CI_low=lo if supported else None,CI_high=hi if supported else None,status='released' if supported else 'withheld_sparse_event_support'))
            (private/'status.json').write_text(json.dumps(state,indent=2)+'\n')
            print('Completed '+model,flush=True)
        table=pd.DataFrame(rows)
        target,outtable,curveout=draw(table)
        assert all(sha(paths[k])==v for k,v in hashes.items())
        state.update(status='completed',ended_utc=now(),verification='Events and person-time conserved by splitting; score convergence <1e-3; covariance/SE identity; finite ordered CIs; source hashes unchanged',output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [target,outtable,curveout]},limitations=['Post hoc seven-band extension of existing review model','Breslow ties and model-based SE retained; no new robust-SE/PH diagnostic claim','Gray .73 curve reused from existing independently refitted illustrative scenario; theta not newly estimated','Sparse-event HR/CI withheld; dashed connections stop at missing estimates','No raw score reconstruction, external weight adjudication, or full manuscript consistency audit'])
    except Exception as exc:
        state.update(status='failed',ended_utc=now(),failure_reason=type(exc).__name__+'; details only in controlled log')
        raise
    finally:
        (private/'status.json').write_text(json.dumps(state,indent=2)+'\n')
        (public/'provenance.json').write_text(json.dumps(state,indent=2)+'\n')
    print('Completed '+RUN,flush=True)
if __name__=='__main__':
    if '--plot-only' in sys.argv:
        draw(pd.read_csv(ROOT/'tables/supp_Figure_07b_seven_age_HR.csv'))
    else:
        with threadpool_limits(limits=1):main()
