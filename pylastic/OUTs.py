import os
import re
import subprocess
from datetime import datetime

months_ru = {
    "янв": "01",
    "фев": "02",
    "мар": "03",
    "апр": "04",
    "мая": "05",
    "июн": "06",
    "июл": "07",
    "авг": "08",
    "сен": "09",
    "окт": "10",
    "ноя": "11",
    "дек": "12"
}


class OUT:
    def __init__(self, Raw, smilesCode=''):
        self.out_type = 'single_point'
        self.Symmetry = ""
        self.SpaceGroup = 0
        self.ParamsCellBase = []
        self.NumberAtoms = 0
        self.TableFreeAtoms = []
        self.CellPrimitive = False
        self.Functional = []
        self.TableBasisSet = []
        self.DateResearch = ''
        self.Elements = {}
        self.typeOpt = ''
        self.id_chains = ''

        self.ParseHead(Raw)

        self.visual_data = {'smiles': smilesCode}

    def ParseHead(self, text):
        textHead = re.search(
            "(CRYSTAL\s|SLAB(\s)?|POLYMER\s|HELIX\s|NANOTUBE\s|MOLECULE\s)(.*\n)*?"
            "(END)?.*"
            "((\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4})|"
            "([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s\d{4})|"
            "([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d{4})\s\d{2}:\d{2}:\d{2}\s+.*))",
            text)[0]
        lines = textHead.split('\n')
        i = 0
        self.Symmetry = lines[i].split(' ')[0]
        i += 1
        if re.search('^0 0 0', lines[i]):
            i += 1
        if re.search('1 0 0', lines[i]):
            self.SpaceGroup = lines[i + 1].strip()
        else:
            self.SpaceGroup = lines[i].strip()
        i += 1
        pattern = r'^\d+$'  # строка с цифрами
        while i < len(lines) and not (re.search(pattern, lines[i])):
            if re.search('^\s?(\d+[.]?\d+(\s*)?)+', lines[i]):
                self.ParamsCellBase = lines[i] + '; '
                i += 1
                break
            i += 1
        self.NumberAtoms = int(lines[i])
        i += 1
        pattern = r'^\D*$'  # строка цифр
        pattern_coords = r'\s?+(\d+)\s+|([-]?\d+.\d+(\w[-]\d+)?)'
        while i < len(lines) and not (re.search(pattern, lines[i])):
            if re.search(pattern_coords, lines[i]):
                matches = re.findall(pattern_coords, lines[i])
                string = ''
                for match in matches:
                    for tup in match:
                        if tup.__len__() >= 1:
                            if re.search('^E', tup):
                                pass
                            else:
                                string += tup.strip() + ' '
                string = string[:-1]
                # string = matches[1][1]
                # string.strip()
                self.TableFreeAtoms.append(string)
            i += 1
        for row in self.TableFreeAtoms:
            atom = re.search("\d+", row)[0]
            count = self.Elements.get(atom)
            if count == None:
                count = 0
            self.Elements.update({atom: count + 1})
        if self.Elements.items().__len__() >= 1:
            key_list = []
            for key in self.Elements:
                key_list.append(key)
            # убрать кавычки, отсортировать и потом вернуть кавычки
            # for key in key_list:
            #     self.Elements.update({int(key):self.Elements[key]})
            #     self.Elements.pop(key)
            # self.Elements = dict(sorted(self.Elements.items()))
            # for key in key_list:
            #     self.Elements.update({str(key):self.Elements[key]})
            #     self.Elements.pop(key)
        if lines[i] == 'PRIMITIV':
            self.CellPrimitive = True
        # typeOptBlock
        typeOptBlock = []
        pattern = r'^\w+[\s]?$'  # строка из символов без пробелов и т.п.
        while i < len(lines) and (re.search(pattern, lines[i])):
            if lines[i] == 'END':
                i += 1
                if lines[i + 1] == 'END':
                    i += 1
                break
            typeOptBlock.append(lines[i])
            i += 1
        # typeOptBlock = str(typeOptBlock)
        # if re.search('ATOMONLY', typeOptBlock):
        #     self.typeOpt = 'ATOMONLY'
        # elif re.search('CVOLOPT', typeOptBlock):
        #     self.typeOpt = 'CVOLOPT'
        # else:
        #     self.typeOpt = 'Full'
        # if re.search('OPTGEOM', typeOptBlock):
        #     self.out_type = 'optimization'
        # elif re.search('FREQCALC', typeOptBlock):
        #     self.out_type = 'hessian'
        #     if re.search('INTRAMAN', typeOptBlock):
        #         self.out_type = 'raman'
        # else:
        #     self.out_type = 'single_point'
        i += 1
        while not re.search('^\d+\s+\d+', lines[i]):
            i += 1
        while (i <= lines.__len__()):
            self.TableBasisSet.append(lines[i])
            i += 1
            if re.search('.*END.*', lines[i]):
                i += 1
                break
        # Functional
        pattern = r'^\w+[-]?\w+?$'  # строка из символов без пробелов и т.п.
        while i < len(lines):
            if re.search('.*END.*', lines[i]):
                i += 1
                continue
            if re.search(pattern, lines[i]):
                self.Functional.append(lines[i])
                if re.search('.*END.*', lines[i + 1]):
                    break
            i += 1
        removable_vars = ['XLGRID',
                          'xlgrid'
                          'numerica',
                          'NUMERICA']
        for var in removable_vars:
            if var in self.Functional:
                self.Functional.remove(var)

        if re.search("(\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4})", textHead):
            date_string = re.search("(END)?.*(\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4})", textHead)[0]
            date_string = date_string.replace("END", '')
            date_string = date_string.replace("YEKT", '')
            date_string = date_string.replace("SCF", '').strip()
            date_string = re.sub(' +', ' ', date_string)
            date_object = datetime.strptime(date_string, "%a %b %d %H:%M:%S %Y")
            self.DateResearch = date_object.strftime("%Y-%m-%d")

        elif re.search("([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s\d{4})", textHead):
            date_string = \
            re.search("(END)?.*([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s\d{4})",
                      textHead)[0]
            date_string = date_string.replace("END", '')
            for ru_month, num_month in months_ru.items():
                date_string = date_string.replace(ru_month, num_month)
            date_object = parser.parse(date_string, dayfirst=True, fuzzy=True)
            self.DateResearch = date_object.strftime("%Y-%m-%d")
        elif re.search("([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d{4})\s\d{2}:\d{2}:\d{2}\s+.*)", textHead):
            date_string = \
            re.search("(END)?.*([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d{4})\s\d{2}:\d{2}:\d{2}\s+.*)", textHead)[0]
            for ru_month, num_month in months_ru.items():
                date_string = date_string.replace(ru_month, num_month)
            self.DateResearch = '03-03-03'
        else:
            self.DateResearch = '09-09-09'

        # Thu May 23 12:12:12 2020
        # if re.search("(END)?\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4}", textHead):
        #     date_string = re.search("(END)?\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4}", textHead)[0]
        #     date_string = date_string.replace("END", '')
        #     date_string = date_string.replace("YEKT", '')
        #     date_object = datetime.strptime(date_string, "%a %b %d %H:%M:%S %Y")
        #     self.DateResearch = date_object.strftime("%Y-%m-%d")
        # elif re.search("(END)?[а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s+\d{4}", textHead):
        #     re.search("(END)?[а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s+\d{4}",
        #               textHead)[0]
        #
        #
        # try:
        #     date_string = re.search("(END)?\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4}", textHead)[0]
        #     date_string = date_string.replace("END", '')
        #     date_string = date_string.replace("YEKT", '')
        #     date_object = datetime.strptime(date_string, "%a %b %d %H:%M:%S %Y")
        #     self.DateResearch = date_object.strftime("%Y-%m-%d")
        # except Exception:
        #     try:
        #         date_string = re.search("(END)?[а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s+\d{4}", textHead)[0]
        #     except:
        #         date_string = re.search("(END)?[а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s*(\d{4})\s\d{2}:\d{2}:\d{2}\s+.*", textHead)[0]
        #     date_string = date_string.replace("END", '')
        #     try:
        #         date_object = parser.parse(date_string, dayfirst=True, fuzzy=True)
        #         self.DateResearch = date_object.strftime("%Y-%m-%d")
        #     except Exception as e:
        #         for ru_month, num_month in months_ru.items():
        #             date_string = date_string.replace(ru_month, num_month)
        #         date_string = re.sub("^\w{2}\s", "", date_string)
        #         date_string = re.sub("[+]\d{2}\s", "", date_string)
        #         try:
        #             date_object = datetime.strptime(date_string, "%m %d %H:%M:%S %Y")
        #         except:
        #             date_object = datetime.strptime(date_string, "%d %m %Y %H:%M:%S %Z")
        #     self.DateResearch = date_object.strftime("%Y-%m-%d")

    def get_result(self):
        result = {
            'out_type': self.out_type,
            'Symmetry': self.Symmetry,
            'SpaceGroup': self.SpaceGroup,
            'ParamsCellBase': self.ParamsCellBase,
            'NumberAtoms': self.NumberAtoms,
            # 'TableFreeAtoms': self.TableFreeAtoms,
            # 'typeOpt': self.typeOpt,
            'CellPrimitive': self.CellPrimitive,
            # 'TableBasisSet': self.TableBasisSet,
            'Functional': self.Functional,
            'Elements': self.Elements,
            'DateResearch': self.DateResearch,
            'visual_data': self.visual_data,
            'id_chains': self.id_chains
        }
        return result

    @staticmethod
    def get_columns():
        result = {
            'id': {"type": 'text', "alias": "Идентификатор", "visible": "True"},
            'filename': {"type": 'text', "alias": "Имя файла", "visible": "True"},
            'out_type': {"type": 'text', "alias": "Тип аут файла", "visible": "True"},
            'Symmetry': {"type": 'text', "alias": "Симметрия", "visible": "True"},
            'SpaceGroup': {"type": 'text', "alias": "Пространственная группа", "visible": "True"},
            'ParamsCellBase': {"type": 'text', "alias": "Базовые параметры ячейки", "visible": "True"},
            'NumberAtoms': {"type": 'int', "alias": "Кол-во атомов", "visible": "True"},
            # 'TableFreeAtoms',
            # 'typeOpt',
            'CellPrimitive': {"type": 'bool', "alias": "Примитивная ячейка", "visible": "True"},
            # 'TableBasisSet',
            'Functional': {"type": 'jsonb', "alias": "Функционал", "visible": "True"},
            'Elements': {"type": 'jsonb', "alias": "Элементы", "visible": "True"},
            'DateResearch': {"type": 'date', "alias": "Дата исследования", "visible": "True"},
            'visual_data': {"type": 'text', "alias": "Данные для визуализации", "visible": "True"},
            'ID_chains': {"type": 'text', "alias": "Связанные файлы", "visible": "True"}
        }
        return result

    @staticmethod
    def get_smiles_code(source_folder, file):
        try:
            outobject = OUT_pdb(filename=source_folder + file)
            outobject.pbc_convert(pbc='mol')
            outobject.write_pdb(outfile=source_folder + file + '.pdb',
                                infile=source_folder + file,
                                write_conect=False)
        except Exception as e:
            return 'Error in pbc_convert:\n' + str(e)
        if os.name == 'nt':
            command_pdb_to_mol2 = [
                'C:\\Users\\Serejka\\PycharmProjects\\out_pdb\\OpenBabel\\obabel.exe',
                "-ipdb", f"{source_folder}{file}.pdb",
                "-omol2", "-O", f"{source_folder}{file}.mol2"
            ]
            command_mol2_to_smi = [
                'C:\\Users\\Serejka\\PycharmProjects\\out_pdb\\OpenBabel\\obabel.exe',
                f"{source_folder}{file}.mol2",
                "-osmi", "-O", f"{source_folder}{file}.smi"
            ]
            # subprocess.call("C:\\Users\\Serejka\\PycharmProjects\\out_pdb\\OpenBabel\\obabel.exe -ipdb " + source_folder+file+'.pdb' + " -omol2 -O "+ source_folder+file+".mol2", stdout='', shell=True)
            # subprocess.call("C:\\Users\\Serejka\\PycharmProjects\\out_pdb\\OpenBabel\\obabel.exe -imol2 " + source_folder+file+'.mol2' + " -osmi -O "+ source_folder+file+".smi", shell=True)
        else:
            command_pdb_to_mol2 = [
                'obabel',
                "-ipdb", f"{source_folder}{file}.pdb",
                "-omol2", "-O", f"{source_folder}{file}.mol2"
            ]
            command_mol2_to_smi = [
                'obabel',
                f"{source_folder}{file}.mol2",
                "-osmi", "-O", f"{source_folder}{file}.smi"
            ]
            # os.system(
            #     "obabel -ipdb " + source_folder + file + '.pdb' + " -omol2 -O " + source_folder + file + ".mol2")
            # os.system(
            #     "obabel -imol2 " + source_folder + file + '.mol2' + " -osmi -O " + source_folder + file + ".smi")
        process0 = subprocess.Popen(command_pdb_to_mol2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process0.communicate()
        process1 = subprocess.Popen(command_mol2_to_smi, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process1.communicate()

        fileSmi = open(source_folder + file + ".smi")
        smilesCode = fileSmi.read().split('\t')[0]
        if smilesCode.__len__() < 1:
            smilesCode = 'No smiles code'
        fileSmi.close()
        os.remove(source_folder + file + '.pdb')
        os.remove(source_folder + file + '.mol2')
        os.remove(source_folder + file + '.smi')
        return smilesCode

    @staticmethod
    def define_type(Raw):
        if re.search('SUPERCEL', Raw):
            return 'SUPERCEL'
        elif re.search('ELASTCON', Raw):
            return 'ELASTCON'
        elif re.search('Firefly Project', Raw):
            return 'Firefly Project'
        elif re.search('OPTGEOM', Raw):
            return 'optimization'
        elif re.search('FREQCALC', Raw):
            if re.search('INTRAMAN', Raw):
                return 'raman'
            else:
                return 'hessian'
        elif re.search(
                "(CRYSTAL\s|SLAB(\s)?|POLYMER\s|HELIX\s|NANOTUBE\s|MOLECULE\s)(.*\n)*?"
                "(END)?.*"
                "((\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\w+\s+\d{4})|"
                "([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s(\d+)\s\d{2}:\d{2}:\d{2}\s+[+]\d{2}\s\d{4})|"
                "([а-яА-ЯёЁ]+\s(\d+\s)?[а-яА-ЯёЁ]{3}\s(\d{4})\s\d{2}:\d{2}:\d{2}\s+.*))",
                Raw):
            return 'base_out'
        else:
            return 'unknown_type'
