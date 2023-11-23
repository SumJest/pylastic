import math

__all__ = ['dirVec', 'dirVec1', 'dirVec2']


def dirVec(theta, phi):
    return [math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)]


def dirVec1(theta, phi, chi):
    return [math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)]


def dirVec2(theta, phi, chi):
    return [math.cos(theta) * math.cos(phi) * math.cos(chi) - math.sin(phi) * math.sin(chi),
            math.cos(theta) * math.sin(phi) * math.cos(chi) + math.cos(phi) * math.sin(chi),
            - math.sin(theta) * math.cos(chi)]
