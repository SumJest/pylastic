from math import log, sqrt
from typing import Tuple, Any

from pylastic.elater import Elastic


async def __get_properties(elas: Elastic) -> tuple[Any, Any, Any, Any]:
    averages = elas.averages()
    kv = averages[0][0]
    kr = averages[1][0]
    gv = averages[0][2]
    gr = averages[1][2]
    return kv, kr, gv, gr


async def elastic_anisotropy_index(elas: Elastic) -> float:
    kv, kr, gv, gr = await __get_properties(elas)
    return kv / kr + 5 * gv / gr - 6


async def log_elastic_anisotropy_index(elas: Elastic) -> float:
    kv, kr, gv, gr = await __get_properties(elas)
    return sqrt(log(kv/kr) ** 2 + (5 * log(gv/gr)) ** 2)
