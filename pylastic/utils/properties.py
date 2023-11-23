from pylastic.elater import Elastic
from pylastic.schemas import Modulus, Property


async def get_properties(elas: Elastic) -> list[Property]:
    schemas = ["Voigt", "Reuss", "Hill"]
    names = ["Bulk modulus", "Young's modulus", "Shear modulus", "Poisson's ratio"]
    symbols = [["Kv", "Ev", "Gv", "vV"], ["Kr", "Er", "Gr", "vR"], ["Kh", "Eh", "Gh", "vH"]]
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
