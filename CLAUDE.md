# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Rate-equation models of Yb/Er upconversion nanoparticles (UCNPs), built on Dyscolab's [jablonski](https://github.com/dyscolab/jablonski) (spectroscopic systems) and [poincare](https://github.com/dyscolab/poincare) (ODE system definition/simulation). Both are pinned to specific git commits in `pyproject.toml`. Parameter values come mostly from Anderson et al. 2013 (`sources/anderson2013.pdf`). This is thesis research code: there is no test suite or lint config; models are exercised from Jupyter notebooks in `notebooks/`.

## Environment

Managed with pixi (Python 3.13, conda-forge); the package is installed editable into the pixi env.

- `pixi install` — create/update the environment
- `pixi run python ...` / `pixi shell` — run inside the env
- Notebooks use the pixi env's `ipykernel`
- Code has been formatted with `ruff format`; commits follow conventional-commit style (`feat(scope): ...`, `fix:`, `chore:`, `refactor:`)

`build/` is a stale, git-ignored setuptools artifact — ignore it.

## Architecture

### Models are declarative `SpectroscopicSystem` classes

Every model is a jablonski `SpectroscopicSystem` subclass whose class body declares:
- `Parameter`s via `assign(default=<pint quantity>)` — can be expressions of other parameters
- energy levels as `SingletState`s via `initial(energy=..., default=<initial population>)`
- transitions (`Absorption`, `Fluorescence`, `InternalConversion`, `MassAction`) or nested sub-systems that wire states together by keyword (e.g. `EnergyTransferUpconversion(sensitizer_high=Yb2, ...)`)

Sub-systems are composed by instantiating them as class attributes with states/parameters of the parent passed in; poincare flattens this into one ODE system. Parameters are referenced/overridden in simulations by the class attribute (e.g. `AndersonModel.energy_flux`).

### Three generations of the same physics

- `anderson.py` — `AndersonModel`: monolithic reproduction of Anderson's Yb/Er model (level Er6 as one state).
- `anderson_unfolded.py` — `AndersonModelUnfolded`: same, but Er6 split into `Er6s` (S3/2) and `Er6h` (H11/2) linked by `ManifoldThermalization`, with temperature-dependent non-radiative rates. Used by `notebooks/paper.ipynb`.
- `ucnp.py` + `ions/` — `UCNP`: the modular refactor in progress. `ions/erbium.py` (`ErIon`) and `ions/yterbium.py` (`YbIon`) define single-ion systems; `UCNP` composes them and adds the Yb↔Er energy-transfer couplings, referencing states as `Yb.Yb2`, `Er.Er3`, etc.

`EnergyTransferUpconversion` and `ManifoldThermalization` are duplicated across these files; the canonical versions for new code live in `transitions.py` / `ions/erbium.py`.

### `transitions.py` — reusable transition building blocks

Temperature-dependent multiphonon relaxation/excitation sub-systems (energy-gap law with Bose–Einstein occupation). `Base.__init__` auto-computes `_gap_wavenumber` from the `source`/`target` state energies (using the pint `"sp"` spectroscopy context), and the phonon number `p = ceil(gap / phonon_wavenumber)`. `*Ref` variants scale a rate given at reference temperature `T_A`; the non-`Ref` variants use an absolute prefactor.

### `jablonski_patch.py` — simulation layer (use this, not raw poincare)

Local overrides/extensions of poincare/jablonski simulation APIs:
- `Simulator` — immutable builder: `Simulator(model).with_solver(...).with_values({...}).with_transform(..., append=True).with_post(...)`. `with_post` stores metadata that is attached as `ds.attrs` on results.
- `EmissionTransform(kind)` / `LineEnergyPost()` — generators that produce, per radiative transition, an output `line_<transition>` (emission rate) and its energy in attrs. `graph_spectra(ds)` consumes exactly this pairing to draw a Gaussian-broadened spectrum.
- `SteadyState.sweep(sim, variable=..., values=...)` — steady-state sweep returning an `xr.Dataset` with units stripped.
- `PulsedExcitation` — square-pulse excitation via jablonski's `piecewise`.
- `to_dict` / `to_values` — serialize parameter values (pint quantities as tuples) and back.

`spectral_steady_sweep` references an undefined `system` and is broken; prefer `SteadyState.sweep` with an `EmissionTransform`.

### Units

Everything uses pint's application registry (`u = get_application_registry()`), and poincare checks dimensionality. Constants `h`, `c`, `k`, `hc` live in `utils.py`. Level energies are written as wavenumbers times `hc` (e.g. `10200 * hc / u.cm`); `AndersonModel` uses an older `energy_factor` convention instead. Excitation in the Anderson models is `energy_flux` (W/cm²) converted to photon flux at 980 nm.

### Level-numbering conventions

Er levels 1–9 follow Anderson's numbering (Er1 = ⁴I15/2 ground … Er9). Radiative transitions are named `rad<from><to>` with color comments (blue from 8, green from 6, red from 5). Per `notes.md`, the H9/2(8)→I13/2(2) line overlaps the green H11/2(6)→I15/2(1) emission.
