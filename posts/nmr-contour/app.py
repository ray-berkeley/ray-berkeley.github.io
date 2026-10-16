import base64
import csv
import io
import json

import matplotlib
matplotlib.use('Agg')
import holoviews as hv
from bokeh.embed import components
from bokeh.events import RangesUpdate
from bokeh.models import CustomJS
from bokeh.resources import INLINE
import numpy as np
import pynmrstar
from flask import Flask, abort, render_template, request
from spectrum import SpectrumRenderer
from werkzeug.exceptions import HTTPException

hv.extension('bokeh')
app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 64 * 1024 * 1024
matplotlib.rcParams.update({'svg.fonttype': 'none', 'font.family': 'DejaVu Sans'})


def require(condition: bool, message: str) -> None:
    if not condition:
        abort(422, message)


def data_uri(content: bytes, mime: str) -> str:
    return f'data:{mime};base64,{base64.b64encode(content).decode()}'


def plot_result() -> str:
    upload = request.files.get('spectrum')
    settings = dict(request.form) if request.method == 'POST' else {'source': 'demo', 'negative': 'on'}
    result = SpectrumRenderer(settings, upload.read() if upload else b'', upload.filename if upload else '').render()
    positive, negative = result['positive'], result['negative']
    color = settings.get('color', '#214f55')
    linewidth = float(settings.get('linewidth') or .6)
    fontsize = float(settings.get('fontsize') or 10)
    xlim, ylim = result['settings']['xlim'], result['settings']['ylim']
    xlabel, ylabel = result['settings']['xlabel'], result['settings']['ylabel']
    title = result['title']
    interactive = (
        hv.Path(positive).opts(color=color, line_width=linewidth * 96 / 72)
        * hv.Path(negative).opts(color='#b75442', line_width=linewidth * 96 / 72)
    ).opts(
        responsive=True, height=560, xlabel=xlabel, ylabel=ylabel,
        xlim=tuple(sorted(xlim)), ylim=tuple(sorted(ylim)),
        invert_xaxis=True, invert_yaxis=True, title=title,
        tools=['pan', 'wheel_zoom', 'box_zoom', 'reset'],
        active_tools=['pan', 'wheel_zoom'], toolbar='above', show_legend=False,
        fontsize={'labels': fontsize, 'xticks': fontsize, 'yticks': fontsize},
    )
    if request.form.get('equal'):
        interactive = interactive.opts(data_aspect=1)
    model = hv.render(interactive, backend='bokeh')
    model.js_on_event(RangesUpdate, CustomJS(code="""
        for (const axis of ['x', 'y']) {
            const values = [cb_obj[axis + '0'], cb_obj[axis + '1']];
            document.querySelector('[name=' + axis + 'min]').value = Math.min(...values);
            document.querySelector('[name=' + axis + 'max]').value = Math.max(...values);
        }
        document.getElementById('export-notice').hidden = false;
    """))
    plot_script, plot_div = components(model)
    return render_template('plot.html', plot_script=plot_script, plot_div=plot_div, svg=data_uri(result['svg'].encode(), 'image/svg+xml'), settings=data_uri(json.dumps(result['settings'], indent=2).encode(), 'application/json'), name=result['name'], shape=result['shape'], minimum=f"{result['minimum']:.4g}", maximum=f"{result['maximum']:.4g}", base=f"{result['base']:.4g}", paths=result['paths'], x_range=f"{result['x_range'][0]:.3f} → {result['x_range'][1]:.3f}", y_range=f"{result['y_range'][0]:.3f} → {result['y_range'][1]:.3f}")


@app.get('/')
def index():
    return render_template('index.html', initial=plot_result(), bokeh_resources=INLINE.render())


@app.post('/plot')
def plot():
    return plot_result()


@app.post('/crosspeaks')
def crosspeaks():
    upload = request.files.get('star')
    require(upload is not None and bool(upload.filename), 'Choose an NMR-STAR file.')
    entry = pynmrstar.Entry.from_string(upload.read().decode())
    loops = entry.get_loops_by_category('Atom_chem_shift')
    require(bool(loops), 'No Atom_chem_shift loop found.')
    loop_number = SpectrumRenderer(dict(request.form)).number('loop', 1, 1, len(loops))
    require(loop_number.is_integer(), 'Shift list number must be an integer.')
    loop = loops[int(loop_number)-1]
    tags = ['Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Val']
    require(all(tag in loop.tags for tag in tags), 'Shift loop requires entity, residue, atom, and value fields.')
    direct = request.form.get('direct', 'CA').strip().upper()
    indirect = request.form.get('indirect', 'CB').strip().upper()
    require(bool(direct) and bool(indirect), 'Enter both atom names.')
    residues = {}
    for assembly, entity, index, residue, atom, value in loop.get_tag(tags):
        if atom not in (direct, indirect) or value in ('.', '?'):
            continue
        key = (assembly, entity, index, residue)
        shifts = residues.setdefault(key, {})
        require(atom not in shifts, f'Duplicate assignment for {key}, {atom}; resolve ambiguity before exporting.')
        shifts[atom] = float(value)
        require(np.isfinite(shifts[atom]), 'Nonfinite chemical shift.')
    rows = [(*key, direct, values[direct], indirect, values[indirect]) for key, values in residues.items() if direct in values and indirect in values]
    require(bool(rows), 'No residues contain both selected atoms in this shift list.')
    columns = ['Assembly', 'Entity', 'Residue index', 'Residue', 'Direct atom', 'Direct ppm', 'Indirect atom', 'Indirect ppm']
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(columns)
    writer.writerows(rows)
    return render_template('crosspeaks.html', rows=rows[:100], count=len(rows), columns=columns, loops=len(loops), csv=data_uri(output.getvalue().encode(), 'text/csv'))


@app.errorhandler(ValueError)
def validation_error(error):
    return render_template('error.html', message=str(error)), 422


@app.errorhandler(HTTPException)
def request_error(error):
    return render_template('error.html', message=error.description), error.code


if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5050, threaded=False)
