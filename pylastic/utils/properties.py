from pylastic.elater import Elastic
from pylastic.schemas import Modulus, Property


async def get_properties(elas: Elastic) -> list[Property]:
    schemas = ["Voigt", "Reuss", "Hill"]
    names = ["Bulk modulus", "Young's modulus", "Shear modulus", "Poisson's ratio"]
    symbols = [["K<sub>V</sub>", "E<sub>V</sub>", "G<sub>V</sub>", "v<sub>V</sub>"],
               ["K<sub>R</sub>", "E<sub>R</sub>", "G<sub>R</sub>", "v<sub>R</sub>"],
               ["K<sub>H</sub>", "E<sub>H</sub>", "G<sub>H</sub>", "v<sub>H</sub>"]]
    units = [["GPa" if i < 3 else '' for i in range(4)] for i in range(3)]
    averages = elas.averages()
    properties: list[Property] = []
    for i in range(len(averages)):
        schema = schemas[i]
        modules = []
        for j in range(len(averages[i])):
            modules.append(Modulus(name=names[j], symbol=symbols[i][j], value=averages[i][j], unit=units[i][j]))
        properties.append(Property(schema=schema, modules=modules))
    return properties
