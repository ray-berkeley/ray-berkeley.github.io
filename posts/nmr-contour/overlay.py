import io

from matplotlib.figure import Figure
from matplotlib import colormaps
from matplotlib.colors import to_hex
import numpy as np
from spectrum import require


SEQUENTIAL_MAPS = {'Blues', 'Greens', 'Purples', 'Oranges', 'BuPu', 'PuBu', 'GnBu', 'BuGn', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'PuBuGn', 'viridis', 'plasma', 'inferno', 'magma', 'cividis'}


def reference_offset(item: dict, axis: str) -> float:
    offset = float(item['settings'].get(f'{axis}offset') or 0)
    require(np.isfinite(offset) and abs(offset) <= 1e5, 'Reference offset must be between -100000 and 100000 ppm.')
    return offset - float(item['result']['settings'].get(f'{axis}offset') or 0)


def contour_colors(item: dict) -> dict:
    style, result = item['settings'], item['result']
    name = style.get('colormap', 'solid')
    require(name == 'solid' or name in SEQUENTIAL_MAPS, 'Choose a sequential colormap.')
    colors = dict(positive=[style['color']] * len(result['positive']), negative=[style['negative_color']] * len(result['negative']), gradient=[])
    if name != 'solid':
        cmap = colormaps[name + ('_r' if style.get('reverse_colormap') else '')]
        levels = sorted(value for value in result['settings']['contour_levels'] if value > 0)
        palette = {level: to_hex(cmap(index / max(1, len(levels) - 1))) for index, level in enumerate(levels)}
        colors['positive'] = [palette[level] for level in result['positive_levels']]
        colors['gradient'] = [to_hex(cmap(index / 15)) for index in range(16)]
    return colors


def compose(datasets: list[dict], settings: dict) -> dict:
    visible = [item for item in datasets if item['settings'].get('visible')]
    require(bool(visible), 'Show at least one dataset.')
    spectra = [item for item in visible if item['kind'] == 'spectrum']
    markers = []
    for item in visible:
        if item['kind'] != 'assignments':
            continue
        style = item['settings']
        points = []
        for row in item['result']['rows']:
            label = f'{row[3]}{row[2]} ({row[0]}/{row[1]})'
            pairs = [(row[5], row[7])]
            if style.get('mirror'):
                pairs.append((row[7], row[5]))
            for x, y in pairs:
                points.append([x + float(style['xoffset']), y + float(style['yoffset']), label])
        markers.append(dict(id=item['id'], points=points, settings=style))
    ranges = {}
    for axis, coordinate in [('x', 0), ('y', 1)]:
        values = [value + reference_offset(item, axis) for item in spectra for value in item['result'][f'{axis}_range']]
        values += [point[coordinate] for marker in markers for point in marker['points']]
        low, high = min(values), max(values)
        if low == high:
            low, high = low - 1, high + 1
        lower = float(settings.get(f'{axis}min') or low)
        upper = float(settings.get(f'{axis}max') or high)
        require(lower < upper, f'{axis}: minimum must be below maximum.')
        ranges[f'{axis}lim'] = [upper, lower]
    fig = Figure(figsize=(float(settings['width']), float(settings['height'])), layout='constrained')
    ax = fig.subplots()
    colors = {item['id']: contour_colors(item) for item in spectra}
    for layer, item in enumerate(reversed(visible), start=2):
        style = item['settings']
        if item['kind'] == 'spectrum':
            offset = [reference_offset(item, axis) for axis in ['x', 'y']]
            for sign in ['positive', 'negative']:
                for index, (segment, color) in enumerate(zip(item['result'][sign], colors[item['id']][sign], strict=True)):
                    ax.plot(*(np.asarray(segment) + offset).T, color=color, linewidth=float(style['linewidth']), zorder=layer, gid=f"dataset-{item['id']}-{sign}-{index}")
        else:
            marker = next(marker for marker in markers if marker['id'] == item['id'])
            for index, (x, y, label) in enumerate(marker['points']):
                ax.plot(x, y, marker='+', color=style['color'], markersize=6, linestyle='none', zorder=layer, gid=f"dataset-{item['id']}-marker-{index}")
                if style.get('labels'):
                    ax.annotate(label, (x, y), xytext=(4, 4), textcoords='offset points', fontsize=float(settings['fontsize']), color=style['color'], clip_on=True, zorder=layer, gid=f"dataset-{item['id']}-label-{index}")
    labels = {f'{axis}label': settings.get(f'{axis}label') or (spectra[0]['result']['settings'][f'{axis}label'] if spectra else f'{axis.upper()} (ppm)') for axis in ['x', 'y']}
    ax.set(**ranges, **labels, title=settings['title'])
    ax.tick_params(labelsize=float(settings['fontsize']))
    for text in [ax.title, ax.xaxis.label, ax.yaxis.label]:
        text.set_fontsize(float(settings['fontsize']))
    if settings.get('equal'):
        ax.set_aspect('equal')
    output = io.StringIO()
    fig.savefig(output, format='svg')
    return dict(svg=output.getvalue(), colors=colors, markers=markers, settings=dict(settings, **ranges, **labels))
