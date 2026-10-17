"""Render carrier/non-carrier counts from the existing aggregate table."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
sys.path.insert(0,str(ROOT))
from scripts.functions.plotting import PALETTE

def main():
    source=ROOT/'tables/55_13gene_carrier_frequency_by_diagnosis_age.csv'
    frame=pd.read_csv(source)
    totals=np.array([155,737,747,1004,1092,1240,1358,1663])
    np.testing.assert_array_equal(frame.n_cases,totals)
    values=[frame.n_cases-frame.n_carriers,frame.n_carriers]
    assert all((v>=20).all() for v in values)
    target=OUT/'Figure_04b_carrier_rate_by_diagnosis_age.svg'
    archive=ROOT/'manuscript/archive/2026-10-17_figure04b_rate'
    archive.mkdir(parents=True,exist_ok=True)
    backup=archive/'Figure_04b_previous_rate.svg'
    tables=ROOT/'tables'
    tables.mkdir(parents=True,exist_ok=True)
    audit=ROOT/'docs/age_interaction_runs/AGE-PRS-BARS-20261017-01'
    audit.mkdir(parents=True,exist_ok=True)
    if target.exists() and not backup.exists():
        shutil.copy2(target,backup)
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','axes.grid':False})
    fig,ax=plt.subplots(figsize=(11,5.8))
    x=np.arange(8);bottom=np.zeros(8)
    for counts,label,color in zip(values,['Non-carrier','Carrier'],['#86B6D8','#D8893D']):
        ax.bar(x,counts,bottom=bottom,width=.72,color=color,edgecolor='white',linewidth=.6,label=label)
        for pos,n,base in zip(x,counts,bottom):
            ax.text(pos,base+n/2,f'{n:,}',ha='center',va='center',fontsize=9,color='#222222')
        bottom+=counts
    labels=['<35','35–<45','45–<50','50–<55','55–<60','60–<65','65–<70','≥70']
    ax.set_xticks(x,[f'{label}\n(n={n:,})' for label,n in zip(labels,totals)])
    ax.set_xlabel('Age at diagnosis (years)',labelpad=13)
    ax.set_ylabel('Count');ax.set_ylim(0,1800);ax.set_xlim(-.65,7.65)
    ax.spines[['top','right']].set_visible(False)
    ax.legend(ncol=2,frameon=False,loc='upper left')
    fig.text(.5,.025,'13-gene carrier status among 7,996 cases; segment labels show counts.',ha='center',fontsize=10)
    fig.tight_layout(rect=[0,.07,1,1])
    fig.savefig(target,facecolor='white')
    fig.savefig('/tmp/Figure_04b_preview.png',dpi=150,facecolor='white')
    plt.close(fig)
    pd.DataFrame({'age_group':labels,'non_carrier':values[0],'carrier':values[1],'n_cases':totals}).to_csv(tables/'Figure_04b_plot_data.csv',index=False)
    (audit/'Figure_04b_provenance.json').write_text(json.dumps(dict(
        run_id='AGE-PRS-BARS-20261017-01',status='completed',source=str(source.relative_to(ROOT)),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        verification='Exact requested age totals; carrier and non-carrier counts sum to age total; all cells >=20',
        previous_rate_figure_backup=str(backup.relative_to(ROOT)),models_refitted=False),indent=2)+'\n')
    print('Figure 04b completed; exact age totals verified.')

if __name__=='__main__':main()
