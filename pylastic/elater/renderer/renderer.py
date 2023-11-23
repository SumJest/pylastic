import logging
import pathlib
from typing import Callable

import numpy as np
from mako.template import Template
from plotly.subplots import make_subplots

from pylastic.elater import Analyzer
from pylastic.elater.renderer.functions import ElasticFunction
from plotly import graph_objects


class ElasticRender:
    analyzer: Analyzer
    template: str
    height: int
    width: int
    full_html: bool
    accuracy: int
    path_to_plotly_js: str

    def __init__(self, analyzer: Analyzer, height: int = 500, width: int = 500, full_html: bool = False,
                 accuracy: int = 500, path_to_plotly_js: str = 'plotly.min.js'):
        self.analyzer = analyzer
        self.height = height
        self.width = width
        self.full_html = full_html
        self.accuracy = accuracy
        self.path_to_plotly_js = path_to_plotly_js

    async def __render_projection(self, func: Callable, function: ElasticFunction):
        phi = np.linspace(0, 2 * np.pi, self.accuracy)
        v_func = np.vectorize(lambda *args: func(args))
        fig = make_subplots(rows=1, cols=3)

        hovertemplate = "r: %{customdata[0]:.2f}<br>&#966;: %{customdata[1]:.2f}<extra></extra>"

        max_point = 0

        r = v_func(np.pi / 2, phi)
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        customdata = np.stack((r.T, phi.T), axis=-1)

        max_point = max(max_point, max(max(x), max(y)))

        fig.add_trace(graph_objects.Scatter(x=x, y=y, fill=None, hovertemplate=hovertemplate,
                                            customdata=customdata, name='xy'), row=1, col=1)
        r = v_func(phi, 0)
        x = r * np.cos(phi)
        z = r * np.sin(phi)
        customdata = np.stack((r.T, phi.T), axis=-1)
        max_point = max(max_point, max(max(x), max(z)))

        fig.add_trace(graph_objects.Scatter(x=z, y=x, fill=None, hovertemplate=hovertemplate,
                                            customdata=customdata, name='xz'), row=1, col=2)

        r = v_func(phi, np.pi / 2)
        y = r * np.cos(phi)
        z = r * np.sin(phi)
        customdata = np.stack((r.T, phi.T), axis=-1)
        max_point = max(max_point, max(max(y), max(y)))

        fig.add_trace(graph_objects.Scatter(x=z, y=y, fill=None, hovertemplate=hovertemplate,
                                            customdata=customdata, name='yz'), row=1, col=3)

        axis_range = round(1.1 * max_point)

        fig.update_xaxes(range=[-axis_range, axis_range])
        fig.update_yaxes(range=[-axis_range, axis_range])

        fig.update_layout(title=f'Projections of {function.name}', autosize=False,
                          width=self.width * 3, height=self.height)
        return fig.to_html(include_plotlyjs=self.path_to_plotly_js, full_html=self.full_html)

    async def __ploty_render(self, func: Callable, function: ElasticFunction):

        phi = np.linspace(0, np.pi, self.accuracy // 2)

        theta = np.linspace(0, 2 * np.pi, self.accuracy)

        theta, phi = np.meshgrid(theta, phi)

        v_func = np.vectorize(lambda *args: func(args))
        r = v_func(theta, phi)
        x = r * np.cos(phi) * np.sin(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(theta)

        customdata = np.stack((r.T, theta.T, phi.T), axis=-1)
        hovertemplate = "r: %{customdata[0]:.2f}<br>&#952;: %{customdata[1]:.2f}" + \
                        "<br>&#966;: %{customdata[2]:.2f}<extra></extra>"

        fig = graph_objects.Figure(data=[graph_objects.Surface(x=x, y=y, z=z, surfacecolor=r,
                                                               customdata=customdata, hovertemplate=hovertemplate)])
        fig.update_layout(title=f'Surface of {function.name}', autosize=False,
                          width=self.width, height=self.height)
        return fig.to_html(include_plotlyjs=self.path_to_plotly_js, full_html=self.full_html)

    def __get_func_callable(self, function: ElasticFunction) -> Callable:
        match function:
            case ElasticFunction.YOUNG:
                return self.analyzer.elas.Young
            case ElasticFunction.YOUNG_INVERSE:
                return lambda *args: (1000 / self.analyzer.elas.Young(*args))
            case ElasticFunction.LINEAR:
                return self.analyzer.elas.LC
            case ElasticFunction.LINEAR_INVERSE:
                return lambda *args: (1000 / self.analyzer.elas.LC(*args))
            case _:
                raise Exception("Wrong function")

    async def html_surface(self, function: ElasticFunction):
        f = self.__get_func_callable(function)

        return await self.__ploty_render(f, function)

    async def html_projections(self, function: ElasticFunction):
        f = self.__get_func_callable(function)

        return await self.__render_projection(f, function)
