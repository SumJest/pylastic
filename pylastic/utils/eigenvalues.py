import numpy as np

from pylastic.elater import Elastic
from pylastic.schemas import Value, EigenvaluesStiffnessMatrix


async def get_eigenvalues(elas: Elastic) -> EigenvaluesStiffnessMatrix:
    eigenvals = sorted(np.linalg.eig(elas.CVoigt)[0])
    lambdas = [Value(value=eigenval) for eigenval in eigenvals]
    return EigenvaluesStiffnessMatrix(lambdas=lambdas)
