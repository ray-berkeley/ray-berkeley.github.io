import {contourPalettes} from './color-palettes.js';

export function offset(item, axis) {
  const value = Number(item.settings[`${axis}offset`] || 0);
  if (!Number.isFinite(value) || Math.abs(value) > 1e5) throw Error('Reference offset must be between -100000 and 100000 ppm.');
  return value;
}
export function contourColors(item, sign) {
  const {settings, result} = item;
  if (sign === 'negative' || settings.colormap === 'solid') return Array(result[sign].length).fill(sign === 'negative' ? settings.negative_color : settings.color);
  const palette = contourPalettes[settings.colormap + (settings.reverse_colormap ? '_r' : '')];
  if (!palette) throw Error('Choose a sequential colormap.');
  const levels = result.settings.contour_levels.filter(value => value > 0).sort((a, b) => a - b);
  const colors = new Map(levels.map((level, index) => [level, palette[Math.min(255, Math.floor(index / Math.max(1, levels.length - 1) * 256))]]));
  return result.positive_levels.map(level => colors.get(level));
}
export function assignmentPoints(item) {
  const xoffset = offset(item, 'x'), yoffset = offset(item, 'y');
  return item.result.rows.flatMap(row => {
    const label = `${row[3]}${row[2]} (${row[0]}/${row[1]})`;
    const pairs = item.settings.mirror ? [[row[5], row[7]], [row[7], row[5]]] : [[row[5], row[7]]];
    return pairs.map(([x, y]) => [x + xoffset, y + yoffset, label]);
  });
}
export function viewSettings(items, settings) {
  const visible = items.filter(item => item.settings.visible && item.result);
  const spectra = visible.filter(item => item.kind === 'spectrum');
  const markers = visible.filter(item => item.kind === 'assignments').flatMap(assignmentPoints);
  const result = {...settings};
  for (const [coordinate, axis] of ['x', 'y'].entries()) {
    const values = spectra.flatMap(item => item.result[`${axis}_range`].map(value => value + offset(item, axis)));
    values.push(...markers.map(point => point[coordinate]));
    let low = Math.min(...values), high = Math.max(...values);
    if (low === high) { low--; high++; }
    const lower = Number(settings[`${axis}min`] || low), upper = Number(settings[`${axis}max`] || high);
    if (!Number.isFinite(lower) || !Number.isFinite(upper) || lower >= upper) throw Error(`${axis}: minimum must be below maximum.`);
    result[`${axis}lim`] = [upper, lower];
    result[`${axis}label`] = settings[`${axis}label`] || (spectra.length ? spectra[0].result.settings[`${axis}label`] : `${axis.toUpperCase()} (ppm)`);
  }
  return result;
}
export function provenance(item) {
  if (!item.result) return null;
  if (item.kind === 'assignments') return {source_sha256: item.result.source_sha256, parser: item.result.parser};
  const settings = {...item.result.settings};
  for (const axis of ['x', 'y']) {
    settings[`${axis}offset`] = item.settings[`${axis}offset`];
    settings[`${axis}lim`] = settings[`${axis}lim`].map(value => value + offset(item, axis));
  }
  return settings;
}
