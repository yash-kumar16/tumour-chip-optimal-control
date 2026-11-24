"""Reproduce published SN38 input schedule and solve an illustrative CTO problem.

Run: python run_project.py
Provenance and limitations: README.md, data/SOURCES.md
"""
from pathlib import Path
import csv
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from tumour_chip.model import Scenario, SphericalTransport
from tumour_chip.control import (load_published_schedule,
                                  calibrate_kill_rate, ControlProblem)


def write_csv(path, header, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def stairs(edges, vals):
    return np.r_[edges[0], edges[1:]], np.r_[vals, vals[-1]]


def main():
    scenario = Scenario()
    schedule = load_published_schedule(ROOT/'data/petreus2021_microfluidic_24h.csv')
    transport = SphericalTransport(scenario)
    kappa = calibrate_kill_rate(transport, schedule, 0.60)
    control = ControlProblem(transport, schedule, kappa)
    reference = schedule.sn38_nM
    auc = schedule.auc_nM_hours
    maximum = max(reference)
    constant = np.full(len(reference), auc / 24.0)
    optimal_result = control.optimise(maximum_nM=maximum, target_auc_nM_hours=auc)
    optimal = optimal_result.x
    solutions = {'Published schedule': reference, 'Uniform AUC': constant, 'Model-optimised': optimal}

    records=[]
    for name, values in solutions.items():
        records.append([name, control.viability(values), float(np.dot(schedule.durations_hours, values)),
                        float(values.max())])
    write_csv(ROOT/'results/model_comparison.csv',
              ['schedule','predicted_relative_viability','inlet_auc_nM_h','peak_nM'],records)
    write_csv(ROOT/'results/schedules.csv',
              ['segment','start_h','end_h','published_nM','uniform_nM','optimised_nM'],
              [(i+1, schedule.edges_hours[i], schedule.edges_hours[i+1],
                reference[i], constant[i], optimal[i]) for i in range(len(reference))])

    print('SOURCE: Petreus et al. (2021), Table 1; SN38, SW620 microfluidic setup')
    print('WARNING: diffusion, permeability, EC50 and response law are ASSUMPTIONS.')
    print(f'Published 24h inlet AUC: {auc:.4f} nM*h; flow: 20 uL/min.')
    print(f'Conditional calibration to approximate 0.60 relative metabolic viability: kappa={kappa:.5f} day^-1')
    print('Model-predicted survival (NOT measured outcomes):')
    for row in records:
        print(f'  {row[0]:19s}  viability={row[1]:.5f}, AUC={row[2]:.2f}, peak={row[3]:.3f}')

    fig, axes = plt.subplots(2, 1, figsize=(8.4, 7.5), layout='constrained')
    for name, schedule_conc in solutions.items():
        axes[0].stairs(schedule_conc, schedule.edges_hours, label=name, linewidth=2)
    axes[0].set(xlabel='Time after inlet exposure starts (h)', ylabel='Inlet SN38 (nM)',
                xlim=(0,24), ylim=(0,6.1), title='24 h inlet schedules at identical SN38 AUC')
    axes[0].legend(frameon=False, ncol=2, fontsize=9)
    bars=axes[1].bar([x[0] for x in records], [x[1] for x in records])
    axes[1].set(ylabel='Predicted viability relative to control',
                title='Six-day model endpoint (conditional, not experimentally validated)', ylim=(0,1.0))
    for rect,v in zip(bars,[x[1] for x in records]):
        axes[1].text(rect.get_x()+rect.get_width()/2, v+0.015, f'{v:.3f}',
                     ha='center', fontsize=10)
    fig.savefig(ROOT/'figures/inlet_control_comparison.png',dpi=160)
    plt.close(fig)

    # Tissue-side modelled drug levels for the published profile.
    state = control.fields(reference)
    hmax = min(48, transport.times_hours[-1])
    nmax = int(round(hmax * 60.0 / scenario.timestep_minutes))
    fig, ax=plt.subplots(figsize=(8.4,4.2),layout='constrained')
    im=ax.pcolormesh(transport.times_hours[:nmax+1],
                     transport.centres_m * 1e6, state[:nmax+1].T,
                     shading='auto',cmap='viridis', vmin=0, vmax=5.5)
    ax.set(xlabel='Time (hours)',ylabel='Distance from spheroid centre (micrometres)',
           title='Modelled SN38 diffusion into spheroid: published inlet schedule')
    fig.colorbar(im,ax=ax,label='Predicted internal SN38 (nM)')
    fig.savefig(ROOT/'figures/spheroid_concentration_heatmap.png',dpi=160)
    plt.close(fig)

    # Distinguish the 24 hourly values, which include model-derived mouse PK,
    # from in-vitro concentrations actually programmed into the reservoirs.
    invivo = np.genfromtxt(ROOT/'data/petreus2021_invivo_free_concentration_24h.csv',delimiter=',',
                           names=True,encoding='utf-8')
    fig, ax=plt.subplots(figsize=(8.4,3.5),layout='constrained')
    ax.plot(invivo['time_hour'],invivo['sn38_nM'], 'o-',label='Paper: model-supported mouse free-plasma PK')
    ax.stairs(reference,schedule.edges_hours,linewidth=2,label='Paper: microfluidic reservoir settings')
    ax.set(xlabel='Time (h)',ylabel='SN38 (nM)',xlim=(0,24),ylim=(0,6.0),
           title='Reproduction of Petreus et al. (2021), Table 1 values')
    ax.legend(frameon=False,fontsize=8)
    fig.savefig(ROOT/'figures/published_petreus_data.png',dpi=160)
    plt.close(fig)

    with open(ROOT/'results/model_assumptions.txt','w',encoding='utf-8') as f:
        f.write('All values below are scenario assumptions unless noted as sourced.\n')
        f.write(f'R={scenario.radius_m} m: chosen midpoint of 290-370 um reported diameter.\n')
        f.write(f'D={scenario.diffusion_m2_s} m^2/s: UNMEASURED ASSUMPTION.\n')
        f.write(f'P={scenario.transfer_m_s} m/s: UNMEASURED ASSUMPTION for Matrigel/tissue boundary, NOT Pearce membrane.\n')
        f.write(f'EC50={scenario.ec50_nM} nM: UNMEASURED ASSUMPTION.\n')
        f.write(f'clearance={scenario.clearance_per_day} /day: ASSUMED zero.\n')
        f.write(f'conditional fitted kappa={kappa:.10f} /day from approximate 0.60 viability target (Fig.2d).\n')
        f.write('Endpoint is a proxy for relative CellTiter-Glo signal, not cell counts.\n')
    return records

if __name__=='__main__':
    main()
