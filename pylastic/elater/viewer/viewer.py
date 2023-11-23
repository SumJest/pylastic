from mako.template import Template

from pylastic.schemas import Analyzed
import pathlib


class Viewer:
    analyzed: Analyzed
    template: str

    def __init__(self, data: Analyzed):
        self.analyzed = data
        self.__read_template()

    def __read_template(self):
        with open(pathlib.Path(__file__).parent.resolve() / 'templates' / 'result.html', 'r') as fio:
            self.template = fio.read()

    async def render(self):
        template = Template(self.template)
        return template.render(**self.analyzed.model_dump())
