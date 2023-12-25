import enum


class ElasticFunction(enum.Enum):
    YOUNG = 0, "Young's modulus"
    YOUNG_INVERSE = 1, "Uniaxial compressibility"
    LINEAR = 2, "Linear compressibility"
    LINEAR_INVERSE = 3, "Linear stiffness"

    def __init__(self, num, title):
        self.num = num
        self.title = title

    @staticmethod
    def by_num(num: int):
        for member in ElasticFunction:
            if member.num == num:
                return member
        raise ValueError(f"{num} is not a valid ElasticFunction")

    @staticmethod
    def by_name(name: str):
        return getattr(ElasticFunction, name)