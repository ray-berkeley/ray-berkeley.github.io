import hashlib
import io
import re
from importlib.metadata import version
from functools import cached_property
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile, is_zipfile

import matplotlib
matplotlib.use('Agg')
import nmrglue as ng
import numpy as np
from matplotlib.figure import Figure

matplotlib.rcParams.update({'svg.fonttype': 'none', 'font.family': 'DejaVu Sans'})


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def demo() -> tuple:
    x, y = np.linspace(80, 10, 512), np.linspace(80, 10, 384)
    xx, yy = np.meshgrid(x, y)
    z = np.zeros_like(xx)
    for a, b, amplitude in [(55, 32, 1), (62, 28, .8), (48, 20, .65), (68, 42, .5)]:
        z += amplitude * np.exp(-((xx-a)/.65)**2-((yy-b)/.85)**2)
        z += amplitude * np.exp(-((xx-b)/.65)**2-((yy-a)/.85)**2)
    z -= .3 * np.exp(-((xx-40)/.8)**2-((yy-60)/.8)**2)
    return z, x, y, '¹³C', '¹³C', [], 'Synthetic demonstration', None


class SpectrumRenderer:
    def __init__(self, settings: dict, payload: bytes = b'', filename: str = ''):
        self.settings = settings
        self.payload = payload
        self.filename = filename

    @cached_property
    def data(self) -> tuple:
        data = self.read_spectrum()
        self.payload = b''
        return data

    @cached_property
    def automatic_base(self) -> float:
        z = self.data[0]
        return float(np.median(z) + 1.5 * np.std(z))

    @cached_property
    def intensity_range(self) -> tuple[float, float]:
        z = self.data[0]
        return float(z.min()), float(z.max())

    def number(self, name: str, default: float, low: float, high: float) -> float:
        raw = self.settings.get(name, '').strip()
        require(not raw or bool(re.fullmatch(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?', raw)), f'{name}: enter a number.')
        value = float(raw) if raw else float(default)
        require(np.isfinite(value) and low <= value <= high, f'{name}: enter a value from {low} to {high}.')
        return value


    def read_spectrum(self) -> tuple:
        source = self.settings.get('source', 'demo')
        if source == 'demo':
            return demo()
        require(bool(self.payload) and bool(self.filename), 'Choose a spectrum file.')
        payload = self.payload
        require(len(payload) <= 64 * 1024**2, 'Spectrum upload exceeds 64 MB.')
        if source == 'auto':
            if payload.startswith(b'UCSF NMR'):
                source = 'ucsf'
            elif is_zipfile(io.BytesIO(payload)):
                source = 'bruker'
            else:
                raise ValueError('Choose a UCSF spectrum or a Bruker ZIP containing 2rr, procs, and proc2s. A standalone 2rr file is not sufficient.')
            self.settings['source'] = source
        require(source in ('ucsf', 'bruker'), 'Choose a UCSF spectrum or a Bruker ZIP.')
        digest = hashlib.sha256(payload).hexdigest()
        with TemporaryDirectory() as folder:
            root = Path(folder)
            levels = []
            if source == 'ucsf':
                path = root / 'spectrum.ucsf'
                path.write_bytes(payload)
                metadata, z = ng.sparky.read(str(path))
                require(z.ndim == 2, 'Only processed 2D spectra are supported.')
                x = ng.sparky.make_uc(metadata, z, dim=1).ppm_scale()
                y = ng.sparky.make_uc(metadata, z, dim=0).ppm_scale()
                xlab, ylab = metadata['w2']['nucleus'], metadata['w1']['nucleus']
            else:
                with ZipFile(io.BytesIO(payload)) as archive:
                    members = [m for m in archive.infolist() if not m.is_dir() and '__MACOSX' not in m.filename]
                    require(sum(m.file_size for m in members) <= 256 * 1024**2, 'Unpacked ZIP exceeds 256 MB.')
                    spectra = [m for m in members if Path(m.filename).name == '2rr']
                    require(len(spectra) == 1, 'ZIP must contain exactly one processed dataset (one 2rr file).')
                    parent = Path(spectra[0].filename).parent
                    for name in ('2rr', 'procs', 'proc2s', 'clevels'):
                        matches = [m for m in members if Path(m.filename) == parent / name]
                        require(len(matches) <= 1, f'Duplicate {name} in ZIP.')
                        if matches:
                            (root / name).write_bytes(archive.read(matches[0]))
                    require(all((root / name).exists() for name in ('2rr', 'procs', 'proc2s')), 'ZIP requires 2rr, procs, and proc2s in the same directory.')
                metadata, z = ng.bruker.read_pdata(str(root), bin_files=['2rr'], read_acqus=False, scale_data=True)
                require(z.ndim == 2, 'Only processed 2D spectra are supported.')
                axes = []
                for key, size in [('procs', z.shape[1]), ('proc2s', z.shape[0])]:
                    params = metadata[key]
                    require(all(k in params for k in ('OFFSET', 'SW_p', 'SF', 'SI')), f'{key} requires OFFSET, SW_p, SF, and SI.')
                    require(params['SF'] > 0 and params['SW_p'] > 0 and params['SI'] == size, 'Unsupported Bruker axis metadata or cropped data; export UCSF instead.')
                    axes.append(params['OFFSET'] - np.arange(size) * params['SW_p'] / params['SF'] / size)
                x, y = axes
                xlab, ylab = metadata['procs'].get('AXNUC', 'F2'), metadata['proc2s'].get('AXNUC', 'F1')
                if (root / 'clevels').exists():
                    for line in (root / 'clevels').read_text().splitlines():
                        if re.match(r'^\s*[+-]?(?:\d|\.\d)', line):
                            levels.extend(float(value) for value in line.split())
        require(min(z.shape) >= 2 and z.size <= 16_000_000, 'Spectrum must have at most 16 million points and at least two points on each axis.')
        require(np.isfinite(z).all(), 'Spectrum contains nonfinite intensities.')
        return z, x, y, xlab, ylab, levels, Path(self.filename).name, digest


    def limits(self, axis: np.ndarray, prefix: str) -> tuple[float, float]:
        low = self.number(f'{prefix}min', float(axis.min()), -1e6, 1e6)
        high = self.number(f'{prefix}max', float(axis.max()), -1e6, 1e6)
        require(low < high, f'{prefix}: minimum must be below maximum.')
        require(low < axis.max() and high > axis.min(), f'{prefix}: range does not overlap the spectrum.')
        return high, low


    def render(self, export_svg: bool = True) -> dict:
        z, x, y, xlab, ylab, saved, name, digest = self.data
        x = x + self.number('xoffset', 0, -1e5, 1e5)
        y = y + self.number('yoffset', 0, -1e5, 1e5)
        mode = self.settings.get('levels', 'auto')
        require(mode in ('auto', 'manual', 'saved'), 'Unknown contour mode.')
        require(mode != 'manual' or bool(self.settings.get('base', '').strip()), 'Enter a manual base intensity.')
        factor = self.number('factor', 1.4, 1.001, 10)
        count = self.number('count', 12, 1, 40)
        require(count.is_integer(), 'Contour count must be an integer.')
        if mode == 'saved':
            require(bool(saved), 'No clevels file found. Use automatic or manual levels.')
            levels = np.array(sorted(set(v for v in saved if v != 0)))
            require(len(levels) <= 80 and np.isfinite(levels).all(), 'Invalid or excessive saved contour levels.')
        else:
            if mode == 'auto':
                base = self.automatic_base
                require(base > 0, 'Automatic base is not positive. Set a manual base contour.')
            else:
                base = self.number('base', 1, 1e-30, 1e30)
            levels = base * factor ** np.arange(int(count))
            if self.settings.get('negative'):
                levels = np.concatenate((-levels[::-1], levels))
        xlim, ylim = self.limits(x, 'x'), self.limits(y, 'y')
        width = self.number('width', 6, 2, 16)
        height = self.number('height', 6, 2, 16)
        linewidth = self.number('linewidth', .6, .1, 3)
        fontsize = self.number('fontsize', 10, 5, 24)
        color = self.settings.get('color', '#214f55')
        require(bool(re.fullmatch(r'#[0-9a-fA-F]{6}', color)), 'Choose a valid contour color.')
        xlabel = self.settings.get('xlabel', '').strip() or f'{xlab} (ppm)'
        ylabel = self.settings.get('ylabel', '').strip() or f'{ylab} (ppm)'
        title = self.settings.get('title', '').strip() or name
        fig = Figure(figsize=(width, height), layout='constrained')
        ax = fig.subplots()
        contour = ax.contour(x, y, z, levels=levels, colors=[('#b75442' if v < 0 else color) for v in levels], linewidths=linewidth, negative_linestyles='solid')
        # One artist per continuous contour gives separately editable SVG paths.
        contour.remove()
        paths = 0
        positive, negative = [], []
        positive_levels, negative_levels = [], []
        for level, segments in zip(contour.levels, contour.allsegs, strict=True):
            for segment in segments:
                if len(segment) >= 2:
                    if export_svg:
                        ax.plot(segment[:, 0], segment[:, 1], color='#b75442' if level < 0 else color, linewidth=linewidth)
                    (negative if level < 0 else positive).append(segment)
                    (negative_levels if level < 0 else positive_levels).append(float(level))
                    paths += 1
        ax.set(xlim=xlim, ylim=ylim)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_title(title, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        if self.settings.get('equal'):
            ax.set_aspect('equal')
        settings = dict(self.settings)
        settings.update(source_name=name, source_sha256=digest, contour_levels=levels.tolist(), xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel, shape=list(z.shape), axes='X = direct (columns); Y = indirect (rows)', software={package: version(package) for package in ('nmrglue', 'numpy', 'matplotlib')}, automatic_rule='median + 1.5 × standard deviation; display heuristic, not a noise estimate')
        result = dict(positive_levels=positive_levels, negative_levels=negative_levels, positive=[s.tolist() for s in positive], negative=[s.tolist() for s in negative], settings=settings, name=name, shape=list(z.shape), minimum=self.intensity_range[0], maximum=self.intensity_range[1], base=float(min(abs(levels))), paths=paths, x_range=[float(x.max()), float(x.min())], y_range=[float(y.max()), float(y.min())], title=title)
        if export_svg:
            output = io.BytesIO()
            fig.savefig(output, format='svg')
            result['svg'] = output.getvalue().decode()
        return result
