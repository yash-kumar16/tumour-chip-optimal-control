# Tumour-on-a-Chip: PDE-Constrained Drug Delivery Optimisation

**Status:** Reproducible research prototype (v0.1). **Not** an experimentally validated predictive model and not clinical dosing advice.

A deliberately limited mathematical/computational project joining **drug transport**, **mathematical biology** and **optimal control theory**. It uses the actual SN38 drug-inlet schedule published by Petreus et al. (2021) for an SW620 colorectal tumour spheroid-on-chip experiment, then asks a mathematical counterfactual:

> Under the *same 24-hour integrated inlet concentration and maximum inlet concentration*, can another time-dependent dosing profile produce greater modelled suppression of a spheroid's relative metabolic-viability proxy?

The first version is a **radial diffusion PDE in the tumour spheroid**, coupled to an assumed concentration-response model and constrained dosing optimisation. A future extension will model **flow in an explicit microchannel and transport across an explicit porous membrane**. Importantly, the **Petreus device is a single-channel Ibidi spheroid chip, not the two-channel porous-membrane geometry in Pearce's advertised PhD**; this first model must not be misrepresented as a reproduction of that PhD or of Petreus' physical device.

## 1. Data: what is actually from a published paper?

**Primary source:** Petreus et al., *Tumour-on-chip microfluidic platform for assessment of drug pharmacokinetics and treatment response*, **Communications Biology 4, 1001 (2021)**, DOI: [10.1038/s42003-021-02526-y](https://doi.org/10.1038/s42003-021-02526-y).

Our transcriptions of **Table 1**:

- `data/petreus2021_microfluidic_24h.csv`: actual programmed microfluidic reservoir concentrations (SN38 and AZD0156), seven nonzero schedule segments totalling 24 hours, reservoirs 2-8. Reservoir 1 has no nonzero drug-concentration row in Table 1. The study describes an eight-reservoir system; Table 1 provides **seven drug-containing segments**. Only **SN38** is used in the present single-drug model.
- `data/petreus2021_invivo_free_concentration_24h.csv`: the reported 24 hourly drug-concentration values from the study's **in vivo free-drug plasma profile**. **These values are a mixture of measured mouse PK and fitted/model-derived estimates, not 24 independent raw LC-MS samples.**
- `data/petreus2021_reported_approximate_outcomes.csv`: approximately 60% relative CellTiter-Glo *metabolic viability* for the six-day **microfluidic SN38** monotherapy experiment (approximately 40% reduction vs untreated), 40% for a **different static plate** experiment, and 55% control-normalised spheroid **volume** for a distinct seven-day experiment. They are different endpoints and setups; no pooling or false precision.

Additional directly reported setup facts: **20 microlitres/min** flow; initial spheroid diameter about **290-370 micrometres**. Reference AUC is `47.2 nM*h`, obtained exactly from the published microfluidic schedule. At 20 uL/min, this is proportional to approximately **56.64 pmol** SN38 passing the inlet during the 24-hour schedule (not the mass absorbed by cells).

The paper explicitly provides **Supplementary Data 1 (ZIP)** containing source data for figures. This repository does **not** claim to contain that file: access to the download was unsuccessful during preparation. The CSVs above are transparent **Table 1 transcriptions** instead. Link to the publication and source package in `data/SOURCES.md`.

## 2. What we model (not directly measured in the paper)

For radial coordinate `r` inside a spherical spheroid of radius `R`, with drug concentration `c(r,t)` in nM:

```math
\frac{\partial c}{\partial t} = D\frac1{r^2}\frac{\partial}{\partial r}\!\left(r^2\frac{\partial c}{\partial r}\right)-k_{\mathrm{clear}}c,\quad 0<r<R.
```

Center symmetry and boundary mass-transfer conditions:

```math
\partial_r c(0,t)=0,\qquad D\partial_r c(R,t)=P\,[u(t)-c(R,t)],\qquad c(r,0)=0.
```

Here `u(t)` is the inlet/bath concentration. The **effective mass-transfer coefficient `P` describes a tissue/Matrigel interface**, not a fitted or derived permeability of a periodically perforated membrane. Likewise `D` and `EC50` are not empirical estimates for SN38 in the exact study; all are explicitly assumed sensitivity parameters (see below).

An *illustrative* relative-viability model is:

```math
V[u] = \frac{3}{R^3}\int_0^R \exp\left[-\kappa\int_0^T\frac{c(r,t)}{EC_{50}+c(r,t)}\,dt_{\rm days}\right]r^2\,dr.
```

This is normalised to the hypothetical untreated exponential-growth control; no actual cell population is tracked or calibrated. In particular **CellTiter-Glo detects metabolic activity rather than counting viable cells**. The model outcome is a **proxy**, not a validated prediction of a live-cell count.

### Explicit scenario assumptions

| Quantity | Value in prototype | Status |
|---|---:|---|
| Spheroid radius `R` | `165 um` | Midpoint of reported 290-370 um starting diameters |
| SN38 tissue diffusivity `D` | `1e-12 m^2/s` | **Unmeasured assumption**, explored in sensitivity sweep |
| External transfer `P` | `1e-7 m/s` | **Unmeasured assumption**, explored in sensitivity sweep |
| SN38 `EC50` | `2.5 nM` | **Unmeasured assumption**; not measured dose-response |
| Tissue clearance | `0 / day` | **Simplifying assumption** |
| Duration | `6 days` | Published six-day comparison window |
| Spherical volume cells | `30` | Numerical resolution, not a biological parameter |
| Time step | `15 min` | Numerical resolution |
| `kappa` | Fitted to ~`V=0.60` in the microfluidic reference arm | **ONE-point conditional calibration**, *not* external validation |

**Identifiability:** one approximate endpoint cannot determine diffusivity, permeability, EC50, drug killing, and growth simultaneously. We fix D/P/EC50 as scenarios and calibrate only an effective `kappa`. We do not infer a unique tumour response mechanism.

## 3. Optimisation problem (CTO)

We choose seven piecewise-constant concentrations over exactly the durations published in Table 1:

```math
\min_{u_1,\ldots,u_7}\;V[u]
\quad\text{subject to}\quad
0\leq u_j\leq5.5\,\mathrm{nM},\qquad
\sum_{j=1}^7\Delta t_j u_j=47.2\,\mathrm{nM\,h}.
```

This is a **direct transcription** of a continuous-time PDE-constrained control problem: the actual decision is a time-varying function, here restricted to seven experimentally relevant constant-in-time segments. The equality holds the integral inlet exposure fixed. Because the paper reports constant flow, **total supplied drug mass is equal** between these schedules (but tissue uptake need not be).

The solver is a conservative **finite-volume method in spherical coordinates** plus backward Euler time integration. Because the transport PDE is linear in `u`, we precompute each schedule segment's unit-input response. This yields exact derivatives of the **discretised** nonlinear objective with respect to the seven concentrations (validated against finite differences). The constrained numerical optimiser is SLSQP.

**NB:** The optimisation is only within this model and this seven-dimensional control family; it does not optimise a cancer patient's therapy. There are no toxicity or treatment-resistance data and no clinical constraints.

## 4. Reproduce everything

```bash
python -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_project.py
python robustness.py
pytest -q
```

Outputs:

- `figures/published_petreus_data.png`: published model-supported mouse PK values versus the paper's stepped microfluidic reservoir settings.
- `figures/spheroid_concentration_heatmap.png`: **modelled**, not measured, radial concentration field.
- `figures/inlet_control_comparison.png`: source inlet schedule, equal-AUC uniform input, and model-optimised input, with conditional viability predictions.
- `results/model_comparison.csv`, `results/schedules.csv`, `results/model_assumptions.txt`: numerical output with traceable units.
- `results/scenario_sensitivity.csv`: illustrative grid of **assumed** diffusivity and transfer values, holding the baseline fitted kill parameter fixed.
- `tests/test_model.py`: source-data checks, conservation, bounds/positivity, calibration, analytical-gradient check, equal-AUC constraints and timestep refinement.

## 5. What we find in the first *conditional simulation*

For the base assumptions above (run `python run_project.py`):

| Inlet profile (all 24 h, same AUC) | Modelled viability relative to control |
|---|---:|
| Published descending SN38 schedule | `0.6000` **by calibration** |
| Uniform `1.9667 nM` schedule | `0.5634` |
| Optimised seven-segment schedule | `0.5628` |

The uniform schedule already captures almost all the objective improvement in **this particular assumed model**; the model-optimised schedule improves on it only by about `0.0006` on the relative-viability scale. The model does **not** reproduce measured responses to either counterfactual schedule (none were provided). Importantly, numerical optimisation is not a discovery of a better experimental/clinical dosing regimen: its conclusions depend on *unmeasured* D/P/EC50 and the assumed Hill-response law. The sweep varying D/P is a start at testing that dependence, not a substitute for biological measurement.

Full numerical derivation: [`docs/DERIVATION.md`](docs/DERIVATION.md).

## 6. Connections and non-connections to the intended PhD

Related: diffusion PDEs, effective interface conditions, biological drug exposure, mathematical optimisation, finite-volume discretisation, sensitivity analysis, transparent mathematical derivations.

Not yet achieved: explicit channel advection, multiscale membrane-hole derivation, experimental fit of tissue diffusivity, true PK-PD model, toxicity and resistance, independently validated prediction, or adjoint PDE derivation. **Do not put these in the CV as completed work.**

Most directly relevant to Philip Pearce's published mathematical paper:

Brennan, Yeo, Pearce & Dalwadi (2026), *Effective permeability conditions for diffusive transport through impermeable membranes with gaps*, Proc. R. Soc. A, DOI: [10.1098/rspa.2025.0703](https://doi.org/10.1098/rspa.2025.0703).

That paper is **not** a source for the numeric values of `P` or `D` used here. It considers a different membrane geometry. A rigorous follow-on project would compare fully resolved flux through a periodically perforated 2D membrane with an asymptotically derived macroscopic interface condition, then use its derived parameter in a separate two-channel PDE-constrained optimisation problem.

## Research questions, papers, and proposed mathematical theory

- **[Research agenda: cell-state switching, extinction versus expected burden, optimal-control structure, and fast-slow reduction](docs/RESEARCH_AGENDA.md)** — proposed work only; explicitly distinguishes implemented and unimplemented components.
- **[Related papers and foundational mathematical references](REFERENCES.md)** — direct DOI, arXiv and publisher links; citation details in [`references.bib`](references.bib).
- **[Mathematical derivation of the implemented finite-volume and discrete-gradient methods](docs/DERIVATION.md)** — present implementation, with no unimplemented adjoint claims.

## 7. Next milestones

1. Secure the actual Petreus **Supplementary Data 1 ZIP**, check raw viabilities and uncertainty/replicate structure and redo parameter fit without relying on rounded text claims.
2. Seek independently measured SN38 transport/response parameters **for a compatible cell line and matrix**. If not available, fit parameters to new data with confidence intervals instead of fabricating them.
3. Add explicit 1D channel **advection-diffusion** coupled to radial tissue transport; estimate flow velocity only if actual channel cross-sectional geometry is documented.
4. Build the separate two-channel **perforated membrane** model and verify reduced/interface flux laws against a resolved calculation.
5. Derive a true adjoint optimality system and verify gradients, then introduce toxicity penalties only with appropriate empirical data.

## Research integrity

*Source paper data* = exact transcribed numerical Table 1 values. *Observed biological endpoint* = approximate value from paper text (one calibration target only). *Model output* = numerical prediction from an intentionally simplified PDE/viability model. **Never label model outputs as measurements or claim model validation on the fitted endpoint.**

No clinical, individual, or patient data used. Code is for academic demonstration.
