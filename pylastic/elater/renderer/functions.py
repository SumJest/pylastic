import enum


class ElasticFunction(enum.Enum):
    YOUNG = 0, "Young's modulus"
    YOUNG_INVERSE = 1, "Uniaxial compressibility"
    LINEAR = 2, "Linear compressibility"
    LINEAR_INVERSE = 3, "Linear stiffness"

    def __init__(self, num, title):
        self.num = num
        self.title = title
