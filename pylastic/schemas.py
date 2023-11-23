from typing import Optional

from pydantic import BaseModel, Field


class Value(BaseModel):
    value: float
    unit: str = Field('GPa')


# TODO: Change value+unit to Value object
class Modulus(BaseModel):
    name: str = Field(...)
    symbol: str = Field(...)
    value: float = Field(...)
    unit: str = Field('')


class Property(BaseModel):
    schema: str = Field(...)
    modules: list[Modulus] = Field(...)


class EigenvaluesStiffnessMatrix(BaseModel):
    lambdas: list[Value] = Field(...)


class Axis(BaseModel):
    vector: list[float] = Field(...)


class Extremum(BaseModel):
    symbol: str = Field(...)
    value: Value = Field(...)
    axes: list[Axis] = Field(...)
    second_axis: Axis | None = Field(default=None)


class Variation(BaseModel):
    name: str = Field(...)
    anisotropy: float | str = Field(...)
    maximum: Extremum = Field(...)
    minimum: Extremum = Field(...)
    is_isotropy: bool = Field(default=False)


class VariationsElasticModuli(BaseModel):
    young: Variation = Field(...)
    young_reverse: Variation = Field(...)
    linear_compressibility: Variation = Field(...)
    linear_compressibility_reverse: Variation = Field(...)
    shear: Variation = Field(...)
    poisson_ratio: Variation = Field(...)


class Analyzed(BaseModel):
    average_properties: list[Property] = Field(...)
    eigenvalues_of_the_stiffness_matrix: EigenvaluesStiffnessMatrix = Field(...)
    variations_elastic_moduli: VariationsElasticModuli | None = Field(...)
    au: float = Field(...)
    al: float = Field(...)
    is_stable: bool = Field(default=True)