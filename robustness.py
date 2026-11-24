"""Exploratory uncertainty sweep; parameters are NOT literature estimates.

Keep calibrated kill parameter fixed to isolate transport sensitivity.
"""
from pathlib import Path
import sys, csv
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from tumour_chip.model import Scenario,SphericalTransport
from tumour_chip.control import load_published_schedule,calibrate_kill_rate,ControlProblem

def main():
    s=load_published_schedule(ROOT/'data/petreus2021_microfluidic_24h.csv')
    base=Scenario()
    k=calibrate_kill_rate(SphericalTransport(base),s,.60)
    vals=[]
    for d in [1e-13,1e-12,1e-11]:
        for p in [1e-8,1e-7,1e-6]:
            sc=Scenario(diffusion_m2_s=d,transfer_m_s=p)
            prob=ControlProblem(SphericalTransport(sc),s,k)
            constant=np.full(7,s.auc_nM_hours/24)
            vi_pub=prob.viability(s.sn38_nM)
            vi_const=prob.viability(constant)
            vals.append([f'{d:.1e}',f'{p:.1e}',f'{vi_pub:.6f}',f'{vi_const:.6f}',f'{vi_pub-vi_const:.6f}'])
    f=ROOT/'results/scenario_sensitivity.csv'
    with f.open('w',newline='',encoding='utf-8') as h:
        w=csv.writer(h);w.writerow(['assumed_D_m2_s','assumed_transfer_m_s',
          'predicted_published_viability','predicted_uniform_viability',
          'published_minus_uniform_viability']);w.writerows(vals)
    print('Uncertainty sweep: positive final column favours uniform schedule (model only)')
    for row in vals:print(' '.join(row))

if __name__=='__main__':main()
