import {assignmentPoints, contourColors, offset} from './plot-data.js';

export class OverlayPlot {
  constructor(selector) {
    this.selector = selector;
    this.element = document.querySelector(selector);
    this.layers = new Map();
  }
  initialize() {
    if (!this.initializing) this.initializing = this.create();
    return this.initializing;
  }
  async create() {
    this.figure = Bokeh.Plotting.figure({name: 'nmr-overlay', x_range: [80, 10], y_range: [80, 10], sizing_mode: 'stretch_width', height: 560, tools: 'pan,wheel_zoom,box_zoom,reset,save', active_scroll: 'wheel_zoom', toolbar_location: 'above'});
    this.figure.js_event_callbacks = {rangesupdate: [new Bokeh.CustomJS({code: "document.querySelector('#spectrum-form').dispatchEvent(new CustomEvent('viewchange'));"})]};
    this.view = await Bokeh.Plotting.show(this.figure, this.selector);
    this.observer = new ResizeObserver(() => this.resize());
    this.observer.observe(this.element);
  }
  resize() {
    if (!this.settings) return;
    const figure = this.figure;
    this.element.style.minHeight = this.settings.equal ? '0' : '';
    if (this.settings.equal) {
      const ratio = Math.abs((this.settings.ylim[0] - this.settings.ylim[1]) / (this.settings.xlim[0] - this.settings.xlim[1]));
      figure.frame_width = Math.max(100, Math.floor(this.element.clientWidth - 100));
      figure.frame_height = Math.round(figure.frame_width * ratio);
      figure.height = figure.frame_height + 100;
    } else {
      figure.frame_width = null; figure.frame_height = null; figure.height = 560;
    }
  }
  add(item) {
    const name = `dataset-${item.id}`;
    const renderers = item.kind === 'spectrum'
      ? ['positive', 'negative'].map(() => {
        const renderer = this.figure.multi_line({name, xs: [], ys: [], line_color: '#005998'});
        renderer.glyph.line_color = {field: 'colors'};
        return renderer;
      })
      : [this.figure.scatter({name, x: [], y: [], marker: 'cross', size: 9}), this.figure.text({name, x: [], y: [], text: [], x_offset: 5, y_offset: -5})];
    const layer = {renderers};
    this.layers.set(item.id, layer);
    return layer;
  }
  update(items, settings) {
    this.element.hidden = false;
    const figure = this.figure;
    for (const id of this.layers.keys()) if (!items.some(item => item.id === id)) this.layers.delete(id);
    for (const item of items.filter(item => item.result)) {
      const layer = this.layers.get(item.id) || this.add(item);
      const offsets = [offset(item, 'x'), offset(item, 'y')];
      const geometryKey = JSON.stringify([offsets, item.settings.mirror]);
      const geometryChanged = layer.result !== item.result || layer.geometryKey !== geometryKey;
      const colorKey = JSON.stringify([item.settings.color, item.settings.negative_color, item.settings.colormap, item.settings.reverse_colormap]);
      if (item.kind === 'spectrum') {
        for (const [index, sign] of ['positive', 'negative'].entries()) {
          const renderer = layer.renderers[index];
          let data = renderer.data_source.data;
          if (geometryChanged) data = {xs: item.result[sign].map(segment => segment.map(point => point[0] + offsets[0])), ys: item.result[sign].map(segment => segment.map(point => point[1] + offsets[1]))};
          if (geometryChanged || layer.colorKey !== colorKey) renderer.data_source.data = {...data, colors: contourColors(item, sign)};
          renderer.glyph.line_width = Number(item.settings.linewidth) * 96 / 72;
          renderer.visible = item.settings.visible;
        }
      } else {
        const [marker, label] = layer.renderers;
        if (geometryChanged) {
          const points = assignmentPoints(item);
          const data = {x: points.map(point => point[0]), y: points.map(point => point[1]), text: points.map(point => point[2])};
          marker.data_source.data = data; label.data_source.data = data;
        }
        marker.glyph.line_color = item.settings.color;
        label.glyph.text_color = item.settings.color;
        label.glyph.text_font_size = `${settings.fontsize}pt`;
        marker.visible = item.settings.visible;
        label.visible = item.settings.visible && item.settings.labels;
      }
      Object.assign(layer, {result: item.result, geometryKey, colorKey});
    }
    const renderers = [...items].reverse().flatMap(item => this.layers.get(item.id)?.renderers || []);
    if (figure.renderers.length !== renderers.length || figure.renderers.some((renderer, index) => renderer !== renderers[index])) figure.renderers = renderers;
    for (const axis of ['x', 'y']) {
      const limits = settings[`${axis}lim`], range = figure[`${axis}_range`];
      const previous = this.settings?.[`${axis}lim`];
      if (!previous || limits.some((value, index) => value !== previous[index])) {
        range.setv({start: limits[0], end: limits[1], reset_start: limits[0], reset_end: limits[1]});
      }
    }
    figure.title.text = settings.title;
    figure.title.text_font_size = `${settings.fontsize}pt`;
    for (const [side, axis] of [['below', 'x'], ['left', 'y']]) for (const model of figure[side]) {
      model.axis_label = settings[`${axis}label`];
      model.axis_label_text_font_size = `${settings.fontsize}pt`;
      model.major_label_text_font_size = `${settings.fontsize}pt`;
    }
    this.settings = settings;
    this.resize();
  }
  currentView() {
    const settings = {...this.settings};
    for (const axis of ['x', 'y']) {
      const range = this.figure[`${axis}_range`];
      settings[`${axis}lim`] = [range.start, range.end];
      settings[`${axis}min`] = String(Math.min(range.start, range.end));
      settings[`${axis}max`] = String(Math.max(range.start, range.end));
    }
    return settings;
  }
  clear() {
    if (this.figure) this.figure.renderers = [];
    this.layers.clear();
    this.element.hidden = true;
  }
  reset() {
    this.clear();
    this.settings = null;
  }
}
