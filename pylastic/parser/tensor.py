import logging
import re


class Tensor:
    type: str
    matrix: str
    found: bool = False

    def __init__(self, file: bytes):
        self.__parse(file)

    def __parse(self, file: bytes):
        lines = file.decode().split('\n')
        extra_fields = {}
        found = False
        tensor_index = -1
        logging.info("PARSER: Start parsing file for tensor matrix")
        for i in range(len(lines)):
            if "SYMMETRIZED ELASTIC CONSTANTS FOR" not in lines[i]:
                continue

            results = re.match("^\s*SYMMETRIZED\s+ELASTIC\s+CONSTANTS\sFOR\s*(?P<type>[\w][\w\s]+[\w])\s*CASE,\s+IN\s+(?P<unit>[\w]+)[\s]*$",
                               lines[i])
            extra_fields = results.groupdict()
            logging.info(f"PARSER: Found matrix with parameters {extra_fields} on line {i}")
            found = True
            tensor_index = i + 1
            break
        self.found = found
        if not found:
            logging.warning("PARSER: Matrix not found")
            return
        self.matrix = "\n".join(lines[tensor_index+1:tensor_index+7])
        self.type = extra_fields['type']
        logging.info("PARSER: Parsing finished")
