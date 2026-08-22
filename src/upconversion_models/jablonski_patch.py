from itertools import chain, pairwise
from typing import Iterable, Mapping

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pint
import scipy.constants as constants
import xarray as xr
from jablonski._typing import Pumper, Time
from jablonski.simulation import lines_to_energies
from jablonski.states import SpectroscopicSystem
from jablonski.util import SpectraKind, emission_transitions
from poincare import Parameter, Simulator, SteadyState, solvers
from poincare.simulator import Components, Initial
from symbolite import Real

from upconversion_models.utils import c, h


def spectral_steady_state_emission(
    system: SpectroscopicSystem,
    excitation: Mapping[Components, Initial],
    solver=solvers.LSODA(),
    kind: SpectraKind = "emission",
):
    lines = {
        f"line_{transition}": transition
        for transition in emission_transitions(system, kind=kind)
    }

    transform = {k: v.radiative_decay.rate_law for k, v in lines.items()}

    sim = Simulator(system, transform=transform, append_transform=True)

    steady = SteadyState(solver=solver)

    ds = steady.solve(sim, values=excitation)

    for line in lines:
        ds.attrs[line] = lines[line].energy_difference

    return ds

def wavelength_to_rgb(wavelength, gamma=0.8):
    """Aproxima el color RGB percibido para una longitud de onda en nm."""
    wavelength = float(wavelength)
    if 380 <= wavelength < 440:
        R, G, B = -(wavelength - 440) / (440 - 380), 0.0, 1.0
    elif 440 <= wavelength < 490:
        R, G, B = 0.0, (wavelength - 440) / (490 - 440), 1.0
    elif 490 <= wavelength < 510:
        R, G, B = 0.0, 1.0, -(wavelength - 510) / (510 - 490)
    elif 510 <= wavelength < 580:
        R, G, B = (wavelength - 510) / (580 - 510), 1.0, 0.0
    elif 580 <= wavelength < 645:
        R, G, B = 1.0, -(wavelength - 645) / (645 - 580), 0.0
    elif 645 <= wavelength <= 780:
        R, G, B = 1.0, 0.0, 0.0
    else:
        R, G, B = 0.0, 0.0, 0.0

    if 380 <= wavelength < 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif 420 <= wavelength < 700:
        factor = 1.0
    elif 700 <= wavelength <= 780:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
    else:
        factor = 0.0

    R = (R * factor) ** gamma if R > 0 else 0.0
    G = (G * factor) ** gamma if G > 0 else 0.0
    B = (B * factor) ** gamma if B > 0 else 0.0
    return (R, G, B)

def graph_spectra(
    system: SpectroscopicSystem,
    excitation: Mapping[Components, Initial],
    width=5,
    solver=solvers.LSODA(),
    kind: SpectraKind = "emission",
):
    def gaussian(x, mu, A, sigma):
        return A * np.exp(-(((x - mu) / sigma) ** 2))

    ds = spectral_steady_state_emission(system, excitation, solver, kind)

    x = np.arange(350, 750, 1)
    y = np.zeros(x.size)

    for line, energy in ds.attrs.items():
        wavelength = 1 / (energy.magnitude / (h * c)).magnitude * 1e7
        y += gaussian(x, wavelength, ds[line].values.item(), width)
        if wavelength < x.max() and wavelength > x.min(): 
            plt.axvline(wavelength, ls="--", color="gray")
            plt.text(x=wavelength -5, y=1.07, s=line.removeprefix("line_"), rotation=90)

    for i in range(len(x) - 1):
        color = wavelength_to_rgb(x[i])
        plt.fill_between(x[i:i + 2], (y/y.max())[i:i + 2], color=color, alpha=0.1)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Normalized intensity [u.a.]")
    return plt.plot(x, y/y.max(), color="black", lw=1)


def spectral_time_resolved_emission(
    system: SpectroscopicSystem,
    excitation: dict[Time, Mapping[Components, Initial | Real | None]],
    save_at: npt.NDArray[np.float64],
    kind: SpectraKind = "emission",
    join_by_energy: bool = False,
    solver=solvers.LSODA(),
) -> xr.Dataset:
    """Single transition square excitation."""

    lines = {
        f"line_{transition}": transition
        for transition in emission_transitions(system, kind=kind)
    }

    transform = {k: v.radiative_decay.rate_law for k, v in lines.items()}

    sim = Simulator(system, transform=transform, append_transform=True)
    ds = piecewise(sim, events=excitation, save_at=save_at, solver=solver)
    if not join_by_energy:
        for line in lines:
            ds.attrs[line] = lines[line].energy_difference
        return ds[list(lines.keys())]
    else:
        return lines_to_energies(lines, ds)


def spectral_steady_sweep(
    system: SpectroscopicSystem,
    variable: Components,
    values: Iterable[Initial],
    kind: SpectraKind = "emission",
    solver=solvers.LSODA(),
):
    lines = {
        f"line_{transition}": transition
        for transition in emission_transitions(system, kind=kind)
    }

    transform = {k: v.radiative_decay.rate_law for k, v in lines.items()}

    sim = Simulator(system, transform=transform, append_transform=True)

    steady = SteadyState(solver=solvers.LSODA())

    results = {v: steady.solve(sim, values={variable: v}) for v in values}

    return xr.Dataset(
        {
            line: xr.DataArray(
                np.array(
                    [results[val][line] for val in values],
                ),
                dims=str(variable),
                coords={str(variable): values},
            )
            for line in lines
        }
    )


def piecewise(
    sim: Simulator,
    *,
    events: dict[Time, Mapping[Components, Initial | Real | None]],
    save_at: npt.NDArray[np.float64],
    solver=solvers.LSODA(),
) -> xr.Dataset:
    adimensionalized_events = {
        k.m_as("s") if isinstance(k, pint.Quantity) else k: v for k, v in events.items()
    }
    event_keys = list(adimensionalized_events.keys())
    t_events = np.sort(event_keys)
    save_at = np.union1d(save_at, t_events)
    pos = np.searchsorted(save_at, t_events)
    save_ats = np.split(save_at, pos + 1)
    t_spans = pairwise(chain((0,), t_events, (save_at[-1],)))

    dss = []
    state = {}
    for t_span, save_at in zip(t_spans, save_ats):
        ds = sim.solve(t_span=t_span, save_at=save_at, values=state, solver=solver)
        for k, v in adimensionalized_events.get(save_at[-1], {}).items():
            if v is None and k in state:
                del state[k]
            else:
                state[k] = v
            # str(k) porque en el output no usamos el objeto Variable aun
            as_str = str(k)
            if as_str in ds:
                ds[as_str].values[-1] = v

        state.update({k: ds[str(k)].values[-1] for k in sim.compiled.variables})
        dss.append(ds)

    ds = xr.concat(dss, dim="time")
    return ds


class InspectableSpectroscopySystem(SpectroscopicSystem):
    def to_dict(self) -> dict[str, pint.Quantity | float | int]:
        cfg_dict = {}
        for val in self._yield(Parameter):
            if isinstance(val.default, pint.Quantity | float | int):
                cfg_dict[val.name] = val.default
        return cfg_dict
