import asyncio

from pylastic.parser import Tensor
from .exceptions import TensorNotFound
from pylastic.elater import Analyzer
from pylastic.schemas import Variation
from pylastic.OUTs import OUT


class ElastPars(OUT):
    analyzer: Analyzer

    def run_analyze(self):
        asyncio.run(self.analyzer.analyze())

    def __init__(self, raw: str):
        super().__init__(raw)
        self.out_type = 'elastcon'
        tensor = Tensor(raw.encode())
        if not tensor.found:
            raise TensorNotFound('Tensor matrix not found')
        self.analyzer = Analyzer(tensor)
        self.run_analyze()

    @staticmethod
    def get_columns():
        result = OUT.get_columns()
        elastcon_parameters = {
            'bulk_modulus': {"type": "real", "alias": "Kh", "visible": "True"},
            'AU': {"type": "real", "alias": "Универсальный индекс упругой анизотропии", "visible": "true"},
            'mechanical_stability': {"type": "text", "alias": "Механическая стабильность", "visible": "true"},
            'young_modulus': {"type": "text", "alias": "Модуль Юнга", "visible": "true"},
            'linear_compressibility': {"type": "text", "alias": "Линейная сжимаемость", "visible": "true"},
            'shear_modulus': {"type": "text", "alias": "Модуль сдвига", "visible": "true"},
            'poisson_ratio': {"type": "text", "alias": "Коэффициент Пуассона", "visible": "true"},
            'dump': {"type": "jsonb", "alias": "Full data", "visible": "false"}
        }
        result.update(elastcon_parameters)
        return result

    def __format_variation(self, variation: Variation):
        return f'{variation.minimum.value.value}\n{variation.maximum.value.value}\n{variation.anisotropy}'

    def __get_result(self):
        hill_properties = self.analyzer.analyzed.average_properties[-1]
        result = {'dump': self.analyzer.analyzed.model_dump(), 'bulk_modulus': hill_properties.modules[0].value,
                  'au': self.analyzer.analyzed.au, 'mechanical_stability': self.analyzer.analyzed.is_stable}
        if self.analyzer.analyzed.is_stable:
            result['young_modulus'] = self.__format_variation(self.analyzer.analyzed.variations_elastic_moduli.young)
            result['linear_compressibility'] = self.__format_variation(self.analyzer.analyzed.variations_elastic_moduli
                                                                       .linear_compressibility)
            result['shear_modulus'] = self.__format_variation(self.analyzer.analyzed.variations_elastic_moduli.shear)
            result['poisson_ratio'] = self.__format_variation(self.analyzer.analyzed.variations_elastic_moduli
                                                              .poisson_ratio)
        else:
            result['young_modulus'] = ''
            result['linear_compressibility'] = ''
            result['shear_modulus'] = ''
            result['poisson_ratio'] = ''
        return result

    def get_result(self):
        result = super().get_result()
        result.update(self.__get_result())
        return result


if __name__ == "__main__":
    with open(r'C:\Users\Roman\PycharmProjects\Site\files\BERQEU_631-plast_EL.out', 'r') as fio:
        parse = ElastPars(fio.read())
        print(parse.get_result())
