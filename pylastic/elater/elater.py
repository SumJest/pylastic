from pylastic.elater.elastic import Elastic
from pylastic.parser import Tensor
from pylastic.schemas import Analyzed
from pylastic.utils.eigenvalues import get_eigenvalues
from pylastic.utils.formulas import elastic_anisotropy_index, log_elastic_anisotropy_index
from pylastic.utils.properties import get_properties
from pylastic.utils.variations import get_variations


class Analyzer:
    elas: Elastic
    tensor: Tensor
    analyzed: Analyzed = None

    def __init__(self, tensor: Tensor):
        self.tensor = tensor
        self.elas = Elastic(tensor.matrix)

    async def analyze(self):
        properties = await get_properties(self.elas)
        eigenvals = await get_eigenvalues(self.elas)

        is_stable = min(eigenvals.lambdas, key=lambda x: x.value).value > 0
        if is_stable:
            variations = await get_variations(self.elas)
        else:
            variations = None

        au = await elastic_anisotropy_index(self.elas)
        al = await log_elastic_anisotropy_index(self.elas)
        self.analyzed = Analyzed(average_properties=properties,
                                 eigenvalues_of_the_stiffness_matrix=eigenvals,
                                 variations_elastic_moduli=variations,
                                 au=au, al=al, is_stable=is_stable)

    def get_result(self):
        return self.analyzed.model_dump()
