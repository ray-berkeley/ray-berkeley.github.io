import re
import unittest
from xml.etree import ElementTree as ET

import numpy as np

from overlay import compose
from spectrum import SpectrumRenderer


class SvgContourTests(unittest.TestCase):
    def assert_contour(self, svg, identity: str, vertices: np.ndarray, xlim: list, ylim: list) -> None:
        group = svg.find(f".//{{*}}g[@id='{identity}']")
        self.assertIsNotNone(group)
        paths = group.findall('{*}path')
        self.assertEqual(len(paths), 1)
        path = paths[0]
        commands = re.findall('[A-Za-z]', path.attrib['d'])
        closed = np.array_equal(vertices[0], vertices[-1])
        self.assertEqual(commands.count('M'), 1)
        self.assertEqual(commands[-1] == 'z', closed)
        self.assertTrue(set(commands) <= {'M', 'L', 'z'})
        actual = np.array(re.findall(r'-?\d+(?:\.\d+)?', path.attrib['d']), dtype=float).reshape(-1, 2)
        if closed and np.array_equal(actual[0], actual[-1]):
            actual = actual[:-1]
        clip_id = path.attrib['clip-path'].removeprefix('url(#').removesuffix(')')
        clip = svg.find(f".//{{*}}clipPath[@id='{clip_id}']/{{*}}rect")
        left, top, width, height = [float(clip.attrib[key]) for key in ['x', 'y', 'width', 'height']]
        points = vertices[:-1] if closed else vertices
        expected = np.column_stack((
            left + (points[:, 0] - xlim[0]) / (xlim[1] - xlim[0]) * width,
            top + (points[:, 1] - ylim[1]) / (ylim[0] - ylim[1]) * height,
        ))
        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-6)

    def test_long_open_and_closed_contours_keep_every_vertex(self) -> None:
        theta = np.linspace(0, 2 * np.pi, 513)
        closed = np.column_stack((np.cos(theta), np.sin(theta)))
        closed[-1] = closed[0]
        opened = closed[:400] / 2
        result = dict(positive=[closed.tolist()], negative=[opened.tolist()], positive_levels=[1], negative_levels=[-1], x_range=[1, -1], y_range=[1, -1], settings=dict(contour_levels=[-1, 1], xlabel='C (ppm)', ylabel='C (ppm)'))
        style = dict(visible=True, color='#005998', negative_color='#b75442', colormap='solid', linewidth='.6', xoffset='.25', yoffset='-.4')
        settings = dict(xmin='-2', xmax='2', ymin='-2', ymax='2', width='6', height='6', title='', fontsize='10', equal=True)
        exported = compose([dict(id=1, kind='spectrum', result=result, settings=style)], settings)
        svg = ET.fromstring(exported['svg'])
        for sign, vertices in [('positive', closed), ('negative', opened)]:
            self.assert_contour(svg, f'dataset-1-{sign}-0', vertices + [.25, -.4], exported['settings']['xlim'], exported['settings']['ylim'])

    def test_single_spectrum_export_keeps_individual_contours(self) -> None:
        result = SpectrumRenderer(dict(source='demo', negative=True)).render()
        svg = ET.fromstring(result['svg'])
        groups = [element for element in svg.iter() if element.attrib.get('id', '').startswith('contour-')]
        self.assertEqual(len(groups), result['paths'])
        for index, group in enumerate(groups):
            # The legacy renderer writes negative levels before positive levels.
            vertices = np.asarray((result['negative'] + result['positive'])[index])
            self.assert_contour(svg, group.attrib['id'], vertices, result['settings']['xlim'], result['settings']['ylim'])


if __name__ == '__main__':
    unittest.main()
