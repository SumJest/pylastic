import math
from typing import Callable

import numpy as np
from scipy import optimize

from pylastic.elater import Elastic
from pylastic.schemas import VariationsElasticModuli, Variation, Extremum, Value, Axis
from pylastic.utils.vectors import *


async def brute_force(func, ranges, resolution: 25) -> list[tuple]:
    step = [(_range[1] - _range[0]) / resolution for _range in ranges]

    result = []

    for i in np.arange(ranges[0][0], ranges[0][1], step[0]):
        for j in np.arange(ranges[1][0], ranges[1][1], step[1]):
            result.append(([i, j], func([i, j])))
    return result


# Functions to minimize/maximize
async def minimize(func, dim):
    if dim == 2:
        r = ((-2 * np.pi, 2 * np.pi), (-2 * np.pi, 2 * np.pi))
        n = 25
    elif dim == 3:
        r = ((0, np.pi), (0, np.pi), (0, np.pi))
        n = 10
        return optimize.brute(func, r, Ns=n, full_output=True, finish=optimize.fmin)[0:2]

    global_minimum = None
    axises = []
    for a in await brute_force(func, r, 15):
        local_minimum = optimize.minimize(func, a[0])
        if global_minimum is not None:
            if global_minimum - local_minimum.fun > 1e-3:
                global_minimum = local_minimum.fun
                axises.clear()
            elif abs(global_minimum - local_minimum.fun) < 1e-3:
                axises.append(local_minimum.x)
        else:
            global_minimum = local_minimum.fun
    # clear

    new_axises = [dirVec(*a) for a in axises]
    for i in range(len(new_axises) - 1):
        if new_axises[i] is None:
            continue
        for j in range(i + 1, len(new_axises)):
            if new_axises[j] is None:
                continue
            result = np.cross(new_axises[i], new_axises[j])
            non_zeros = [a for a in result if abs(a) > 1e-1]

            if len(non_zeros) == 0:
                new_axises[j] = None
    new_axises = [list(map(lambda x: x, a)) for a in new_axises if a is not None]
    # TODO: Фильтровать всякие окружности и тд
    return new_axises, global_minimum


async def maximize(func, dim):
    res = await minimize(lambda x: -func(x), dim)
    return res[0], -res[1]


async def young(Young: Callable) -> Variation:
    maxE = await maximize(Young, 2)
    minE = await minimize(Young, 2)

    anisE = (maxE[1] / minE[1])

    minimum = Extremum(symbol="E", value=Value(value=minE[1]), axes=[Axis(vector=axis) for axis in minE[0]])
    maximum = Extremum(symbol="E", value=Value(value=maxE[1]), axes=[Axis(vector=axis) for axis in maxE[0]])

    return Variation(name="Young's modulus", anisotropy=anisE, maximum=maximum, minimum=minimum)


async def linear(LC: Callable) -> Variation:
    minLC = await minimize(LC, 2)
    maxLC = await maximize(LC, 2)
    is_isotropy = bool(abs(minLC[1] - maxLC[1]) < 1e-3)

    anisLC = (maxLC[1] / minLC[1]) if minLC[1] > 0 else "undefined"

    min_axes = [Axis(vector=axis) for axis in minLC[0]]
    max_axes = [Axis(vector=axis) for axis in maxLC[0]]

    if is_isotropy:
        min_axes = max_axes = []

    minimum = Extremum(symbol="β", value=Value(value=minLC[1], unit="TPa-1"),
                       axes=min_axes)
    maximum = Extremum(symbol="β", value=Value(value=maxLC[1], unit="TPa-1"),
                       axes=max_axes)

    variation = Variation(name="Linear compressibility", anisotropy=anisLC, maximum=maximum, minimum=minimum)
    if is_isotropy:
        setattr(variation, 'is_isotropy', is_isotropy)

    return variation


async def shear(shear: Callable) -> Variation:
    minG = await minimize(shear, 3)
    maxG = await maximize(shear, 3)

    anisG = (maxG[1] / minG[1])

    minimum = Extremum(symbol="G", value=Value(value=minG[1]), axes=[Axis(vector=dirVec1(*minG[0]))],
                       second_axis=Axis(vector=dirVec2(*minG[0])))
    maximum = Extremum(symbol="G", value=Value(value=maxG[1]), axes=[Axis(vector=dirVec1(*maxG[0]))],
                       second_axis=Axis(vector=dirVec2(*maxG[0])))

    return Variation(name="Shear modulus", anisotropy=anisG, maximum=maximum, minimum=minimum)


async def poisson(Poisson: Callable) -> Variation:
    minNu = await minimize(Poisson, 3)
    maxNu = await maximize(Poisson, 3)

    anisNu = (maxNu[1] / minNu[1]) if minNu[1] * maxNu[1] > 0 else "undefined"

    minimum = Extremum(symbol="ν", value=Value(value=minNu[1], unit=''), axes=[Axis(vector=dirVec1(*minNu[0]))],
                       second_axis=Axis(vector=dirVec2(*minNu[0])))
    maximum = Extremum(symbol="ν", value=Value(value=maxNu[1], unit=''), axes=[Axis(vector=dirVec1(*maxNu[0]))],
                       second_axis=Axis(vector=dirVec2(*maxNu[0])))

    return Variation(name="Poisson's ratio", anisotropy=anisNu, maximum=maximum, minimum=minimum)


async def get_variations(elas: Elastic) -> VariationsElasticModuli:
    _young = await young(elas.Young)
    _young_reverse = await young(lambda *args: 1000 / elas.Young(*args))

    _young_reverse.name = 'Uniaxial Compressibility'
    _young_reverse.minimum.value.unit = 'TPa-1'
    _young_reverse.maximum.value.unit = 'TPa-1'

    _linear_compressibility = await linear(elas.LC)

    if _linear_compressibility.anisotropy == 'undefined':
        _linear_compressibility_reverse = _linear_compressibility.model_copy(deep=True)
        _linear_compressibility_reverse.minimum, _linear_compressibility_reverse.maximum = _linear_compressibility.maximum.model_copy(deep=True), _linear_compressibility.minimum.model_copy(deep=True)

        _linear_compressibility_reverse.minimum.value.value = 1000 / _linear_compressibility_reverse.minimum.value.value
        _linear_compressibility_reverse.maximum.value.value = 1000 / _linear_compressibility_reverse.maximum.value.value
    else:
        _linear_compressibility_reverse = await linear(lambda *args: 1000 / elas.LC(*args))

    _linear_compressibility_reverse.name = 'Linear stiffness'
    _linear_compressibility_reverse.minimum.value.unit = 'GPa'
    _linear_compressibility_reverse.maximum.value.unit = 'GPa'

    _shear = await shear(elas.shear)
    _poisson_ratio = await poisson(elas.Poisson)

    return VariationsElasticModuli(young=_young, young_reverse=_young_reverse,
                                   linear_compressibility=_linear_compressibility,
                                   linear_compressibility_reverse=_linear_compressibility_reverse,
                                   shear=_shear, poisson_ratio=_poisson_ratio)
