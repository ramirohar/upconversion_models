from itertools import chain, pairwise
from typing import Iterable, Mapping

import numpy as np
import numpy.typing as npt
import pint
import scipy.constants as constants
import xarray as xr
from poincare import Parameter, Simulator, SteadyState, solvers
from poincare.simulator import Components, Initial
from symbolite import Real

from jablonski._typing import Pumper, Time
from jablonski.states import SpectroscopicSystem
from jablonski.util import SpectraKind, emission_transitions
from jablonski.simulation import lines_to_energies


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
