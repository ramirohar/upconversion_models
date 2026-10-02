from __future__ import annotations

from itertools import chain, pairwise
from typing import (
    Generator,
    Hashable,
    Iterable,
    Mapping,
    Sequence,
    override,
)

import matplotlib.pyplot as plt
import numpy as np
import pint
import pint_xarray
import poincare
import xarray as xr
from jablonski._typing import Time
from jablonski.states import SpectroscopicSystem
from jablonski.util import (
    Fluorescence,
    Phosphorescence,
    RadiativeDecay,
    SpectraKind,
)
import numpy.typing as npt
from numpy.typing import ArrayLike
from poincare import Parameter, System, solvers
from poincare.compile import SystemCompiler
from poincare.simulator import Components
from poincare.types import Initial
from scipy_events import Events
from symbolite import Real
from symbolite.core.value import Value

from upconversion_models.utils import c, h
from poincare.helpers import get_from

u = pint.get_application_registry()


class TransformGen:
    def get(self, system: type[System] | System) -> Mapping[Hashable, Value]: ...


class PostGen:
    def get(self, system: type[System] | System) -> dict: ...


class LineEnergyPost(PostGen):
    def get(self, system) -> dict:
        return {
            f"line_{transition}": transition.energy_difference
            for transition in emission_transitions(system, "emission")
        }


class EmissionTransform(TransformGen):
    def __init__(self, kind: SpectraKind):
        self.kind: SpectraKind = kind

    def get(self, system: type[System] | System):

        lines = {
            f"line_{transition}": transition
            for transition in emission_transitions(system, kind=self.kind)
        }

        transform: dict[Hashable, Value] = {
            k: v.radiative_decay.rate_law for k, v in lines.items()
        }
        return transform


class Simulator(poincare.Simulator):
    def __init__(
        self,
        system: System | type[System],
        /,
        *,
        backend: Backend = "numpy",
        solver: solvers.Solver = solvers.LSODA(),
        append_transform: bool = False,
    ):
        self.model = system
        compiler = SystemCompiler(self.model, backend=backend)
        self.compiled = compiler.compiled
        self.transform = self._compile_transform(None)
        self.values = {}
        self.solver = solver
        self.append_transform = append_transform
        self.post = {}

    def with_values(
        self, values: Mapping[Components, Initial], /, *, append: bool = True
    ) -> Simulator:
        sim = self.__class__.__new__(self.__class__)
        sim.model = self.model
        sim.compiled = self.compiled
        sim.transform = self.transform
        sim.values = values if not append else self.values | values
        sim.solver = self.solver
        sim.post = self.post
        return sim

    def with_solver(self, solver: solvers.Solver, /) -> Simulator:
        sim = self.__class__.__new__(self.__class__)
        sim.model = self.model
        sim.compiled = self.compiled
        sim.transform = self.transform
        sim.values = self.values
        sim.solver = solver
        sim.post = self.post
        return sim

    def with_transform(
        self,
        transform: Sequence[Value]
        | Mapping[Hashable, Value]
        | None
        | TransformGen = None,
        /,
        *,
        append: bool = False,
    ) -> Simulator:
        if isinstance(transform, TransformGen):
            transform = transform.get(self.model)
        sim = self.__class__.__new__(self.__class__)
        sim.model = self.model
        sim.compiled = self.compiled
        sim.transform = sim._compile_transform(transform, append)
        sim.values = self.values
        sim.solver = self.solver
        sim.post = self.post
        return sim

    def with_post(self, post: PostGen | dict) -> Simulator:
        if isinstance(post, PostGen):
            post = post.get(self.model)
        sim = self.__class__.__new__(self.__class__)
        sim.model = self.model
        sim.compiled = self.compiled
        sim.transform = self.transform
        sim.values = self.values
        sim.solver = self.solver
        sim.post = post
        return sim

    def solve(
        self,
        # values: Mapping[Components, Initial | Value] = {},
        *,
        t_span: tuple[float, float] | None = None,
        save_at: ArrayLike | None = None,
        # solver: solvers.Solver = solvers.LSODA(),
        events: Sequence[Events] = (),
        check_dimensionality: bool = True,
    ) -> xr.DataTree:
        ds = super().solve(
            t_span=t_span,
            save_at=save_at,
            events=events,
            check_dimensionality=check_dimensionality,
        )
        ds.attrs = self.post
        return ds


class SteadyState(poincare.SteadyState):
    @override
    def sweep(
        self,
        sim: Simulator,
        /,
        *,
        variable: Components,
        values: Iterable[Initial],
    ):
        results = {v: self.solve(sim, values={variable: v}) for v in values}
        ds = xr.Dataset(
            {
                str(var): xr.DataArray(
                    np.array(
                        [
                            results[v][var].item().magnitude
                            if isinstance(results[v][var].item(), pint.Quantity)
                            else results[v][var].item()
                            for v in values
                        ]
                    ),
                    dims=str(variable),
                    coords={str(variable): values},
                )
                for var in [str(var) for var in sim.transform.output.keys()]
            }
            | {
                "time": xr.DataArray(
                    np.array([results[v]["time"].item() for v in values]),
                    dims=str(variable),
                    coords={str(variable): values},
                ),
                "event": xr.DataArray(
                    np.array([results[v]["event"].item() for v in values]),
                    dims=str(variable),
                    coords={str(variable): values},
                ),
            }
        )
        ds.attrs = sim.post
        return ds


def emission_transitions(
    system: System | type[System],
    kind: SpectraKind = "emission",
) -> Generator[RadiativeDecay, None, None]:
    if kind == "emission":
        include = (Fluorescence, Phosphorescence)
    elif kind == "fluorescence":
        include = Fluorescence
    elif kind == "phosphorescence":
        include = Phosphorescence
    else:
        raise ValueError(f"kind must be {SpectraKind}")

    for transition in system._yield(include):
        if isinstance(transition, RadiativeDecay):
            yield transition


def emission_transform(
    system: SpectroscopicSystem, kind: SpectraKind = "emission"
) -> Mapping[Hashable, Value]:
    lines = {
        f"line_{transition}": transition
        for transition in emission_transitions(system, kind=kind)
    }

    transform = {k: v.radiative_decay.rate_law for k, v in lines.items()}
    return transform


def spectral_steady_state(
    sim: poincare.Simulator,
    excitation: Mapping[Components, Initial],
    kind: SpectraKind = "emission",
) -> xr.DataTree:
    lines = {
        f"line_{transition}": transition
        for transition in emission_transitions(sim.model, kind=kind)
    }

    transform = {k: v.radiative_decay.rate_law for k, v in lines.items()}

    sim = sim.with_transform(transform)

    steady = SteadyState()

    ds = steady.solve(sim, values=excitation)

    for line in lines:
        ds.attrs[line] = lines[line].energy_difference

    return ds


# Update this to work with update
def spectral_steady_sweep(
    sim: Simulator,
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

    steady = SteadyState(solver=solver)

    results = {v: steady.solve(sim, values={variable: v}) for v in values}

    ds = xr.Dataset(
        {
            line: xr.DataArray(
                np.array(
                    [results[val][line].values.item() for val in values],
                ),
                dims=str(variable),
                coords={str(variable): values},
            )
            for line in lines
        }
    )

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
    ds: xr.DataTree,
    width=5,
):
    def gaussian(x, mu, A, sigma):
        return A * np.exp(-(((x - mu) / sigma) ** 2))

    x = np.arange(350, 750, 1)
    y = np.zeros(x.size)

    for line, energy in ds.attrs.items():
        wavelength = 1 / (energy.magnitude / (h * c)).magnitude * 1e7
        y += gaussian(x, wavelength, ds[line].values.item(), width)
        if wavelength < x.max() and wavelength > x.min():
            plt.axvline(wavelength, ls="--", color="gray")
            plt.text(
                x=wavelength - 5, y=1.07, s=line.removeprefix("line_"), rotation=90
            )

    for i in range(len(x) - 1):
        color = wavelength_to_rgb(x[i])
        plt.fill_between(x[i : i + 2], (y / y.max())[i : i + 2], color=color, alpha=0.1)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Normalized intensity [u.a.]")
    return plt.plot(x, y / y.max(), color="black", lw=1)


def piecewise(
    sim: Simulator,
    *,
    events: dict[Time, Mapping[Components, Initial | Real | None]],
    save_at: npt.NDArray[np.float64],
) -> xr.Dataset:
    """jablonski's piecewise, with each segment's start time prepended to its save_at.

    poincare's solve starts integrating at save_at[0] and ignores t_span[0], and
    jablonski's segments start at the first save point after the event, so the
    stretch between the event and that point was never integrated.
    """
    try:
        event_keys = np.array([key.to(u.s).magnitude for key in events.keys()])
    except (AttributeError, pint.DimensionalityError):
        raise pint.PintError(
            "events keys must be pint Quantities and have time dimensionality."
        )

    try:
        adimensional_save_at = save_at.to(u.s).magnitude
    except (AttributeError, pint.DimensionalityError):
        raise pint.PintError(
            "save_at must be pint Quantity and have time dimensionality."
        )
    t_events = np.sort(event_keys)
    adimensional_save_at = np.union1d(adimensional_save_at, t_events)
    pos = np.searchsorted(adimensional_save_at, t_events)
    adimensional_save_ats = np.split(adimensional_save_at, pos + 1)
    adimensional_t_spans = pairwise(chain((0,), t_events, (adimensional_save_at[-1],)))
    dss = []
    state = {}
    for t_span, save_at in zip(adimensional_t_spans, adimensional_save_ats):
        # Integrate from t_span[0], then drop the added point.
        prepend = save_at[0] > t_span[0]
        if prepend:
            save_at = np.r_[t_span[0], save_at]
        ds = sim.with_values(state).solve(
            t_span=np.array(t_span) * u.s, save_at=save_at * u.s
        )
        if prepend:
            ds = ds.isel(time=slice(1, None))
            save_at = save_at[1:]

        for k, v in events.get(save_at[-1] * u.s, {}).items():
            if v is None and k in state:
                del state[k]
            else:
                state[k] = v
            # str(k) porque en el output no usamos el objeto Variable aun
            as_str = str(k)
            if as_str in ds:
                ds[as_str][-1] = v

        state.update({k: ds[str(k)][-1].item() for k in sim.compiled.variables})
        dss.append(ds.pint.dequantify())

    ds = xr.concat(dss, dim="time")

    pint_xarray.setup_registry(u)
    ds = ds.pint.quantify()
    u.force_ndarray_like = False
    return ds


def spectral_time_resolved_emission(
    sim: Simulator,
    excitation: dict[Time, Mapping[Components, Initial | Real | None]],
    save_at: npt.NDArray[np.float64],
    kind: SpectraKind = "emission",
    join_by_energy: bool = False,
) -> xr.Dataset:
    """Single transition square excitation."""
    ds = piecewise(sim, events=excitation, save_at=save_at)
    return ds
    # if not join_by_energy:
    #     for line in lines:
    #         ds.attrs[line] = lines[line].energy_difference
    #     return ds[list(lines.keys())]
    # else:
    #     return lines_to_energies(lines, ds)


class PulsedExcitation:
    def __init__(
        self, pulse_width: pint.Quantity, pulse_height: pint.Quantity, pump: Parameter
    ) -> None:
        self.pulse_width = pulse_width
        self.pulse_height = pulse_height
        self.pump = pump

    def solve(
        self,
        sim: Simulator,
        save_at: npt.NDArray[np.float64],
    ) -> xr.Dataset:

        if isinstance(self.pump, str):
            self.pump = get_from(self.pump, sim.model)

        excitation = {
            0 * u.s: {self.pump: self.pulse_height},
            self.pulse_width: {self.pump: 0 * self.pulse_height.units},
        }
        return piecewise(sim, events=excitation, save_at=save_at)
        

def to_dict(values: Mapping[Node, Initial]):
    d = {}
    for k, v in values.items():
        d = d | {k.__repr__(): v.to_tuple() if isinstance(v, pint.Quantity) else v}
    return d


def to_values(dict: dict[str, Any], model):
    return {
        get_from(param=param, model=model): u.Quantity.from_tuple(value)
        if isinstance(value, tuple)
        else value
        for param, value in dict.items()
    }
