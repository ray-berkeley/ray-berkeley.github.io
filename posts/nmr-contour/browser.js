import {ColorPicker, colorBackground} from './color-picker.js';
import {OverlayPlot} from './plot.js';
import {provenance, viewSettings} from './plot-data.js';

const form = document.querySelector('#spectrum-form');
const datasetList = document.querySelector('#dataset-list');
let selectedId = null;
let draggedId = null;
const status = document.querySelector('#status');
const errors = document.querySelector('#error');
const downloads = document.querySelector('#downloads');
const entityFields = ['name', 'levels', 'base', 'count', 'factor', 'negative', 'color', 'colormap', 'reverse_colormap', 'negative_color', 'linewidth', 'xoffset', 'yoffset', 'loop', 'direct', 'indirect', 'labels', 'mirror'];
const globalFields = ['xmin', 'xmax', 'ymin', 'ymax', 'xlabel', 'ylabel', 'title', 'width', 'height', 'fontsize', 'equal'];
const calculationFields = {spectrum: ['levels', 'base', 'count', 'factor', 'negative'], assignments: ['loop', 'direct', 'indirect']};
const palette = ['#005998', '#c45a28', '#6a4794', '#23845d'];
const exampleView = {xmin: '10', xmax: '80', ymin: '10', ymax: '80', xlabel: '', ylabel: '', title: 'FUS-LC: uniformly labeled + N112 (Murray et al.)', width: '6', height: '6', fontsize: '10', equal: true};
let datasets = [], nextId = 1, worker, ready = false, revision = 0, request = 0, pending = new Map();
let workerEpoch = 0, calculatingEpoch = null, timer = null, frame = null, exporting = false, settingsUrl, svgUrl, csvUrl, csvResult;
const plot = new OverlayPlot('#plot');
const svgButton = document.querySelector('#svg-download');
const reloadExample = document.querySelector('#reload-example');
let loadingExamples = false;
const selected = () => datasets.find(item => item.id === selectedId);
const colorPicker = new ColorPicker(form, settings => {
  if (!selected()) return;
  Object.assign(selected().settings, settings);
  showControls(); updatePills(); schedule();
});
function fieldValue(name) {
  const input = form.elements[name];
  return input.type === 'checkbox' ? input.checked : input.value;
}
function setField(name, value) {
  const input = form.elements[name];
  if (input.type === 'checkbox') input.checked = Boolean(value);
  else input.value = value;
}
function showError(error) {
  console.error(error);
  errors.textContent = String(error).trim().split('\n').at(-1);
  errors.hidden = false;
  downloads.hidden = true;
  status.textContent = 'Correct the input to continue.';
}
function updatePills() {
  if (draggedId !== null) return;
  const focus = document.activeElement;
  const focusId = focus?.closest('.dataset-pill')?.dataset.id;
  const focusAction = focus?.dataset.action;
  datasetList.replaceChildren(...datasets.map((item, index) => {
    const pill = document.createElement('li');
    pill.className = 'dataset-pill';
    pill.dataset.id = item.id;
    pill.classList.toggle('selected', item.id === selectedId);
    pill.classList.toggle('muted', !item.settings.visible);
    const button = (action, text, label) => {
      const element = document.createElement('button');
      element.type = 'button'; element.dataset.action = action;
      element.textContent = text; element.setAttribute('aria-label', label); element.title = label;
      return element;
    };
    const handle = button('drag', '⠿', `Drag ${item.settings.name} to reorder`);
    handle.className = 'drag-handle'; handle.draggable = true;
    const edit = button('select', '', `Edit ${item.settings.name}`);
    edit.className = 'pill-select'; edit.setAttribute('aria-pressed', String(item.id === selectedId));
    const swatch = document.createElement('span'); swatch.className = 'pill-swatch';
    swatch.style.background = colorBackground(item.settings);
    swatch.setAttribute('aria-hidden', 'true');
    const name = document.createElement('span'); name.className = 'pill-name'; name.textContent = item.settings.name;
    edit.append(swatch, name);
    const actions = document.createElement('div'); actions.className = 'pill-actions';
    const visible = button('visible', '', `${item.settings.visible ? 'Hide' : 'Show'} ${item.settings.name}`);
    visible.setAttribute('aria-pressed', String(item.settings.visible));
    visible.innerHTML = `<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M2 12s3.5-7 10-7 10 7 10 7-3.5 7-10 7S2 12 2 12Z"/><circle cx="12" cy="12" r="3"/>${item.settings.visible ? '' : '<path d="m3 3 18 18"/>'}</svg>`;
    const up = button('up', '↑', `Move ${item.settings.name} up`); up.disabled = index === 0;
    const down = button('down', '↓', `Move ${item.settings.name} down`); down.disabled = index === datasets.length - 1;
    const remove = button('remove', '×', `Remove ${item.settings.name}`);
    actions.append(visible, up, down, remove); pill.append(handle, edit, actions);
    return pill;
  }));
  document.querySelector('#layer-hint').hidden = datasets.length < 2;
  if (focusId && focusAction) (datasetList.querySelector(`[data-id="${focusId}"] [data-action="${focusAction}"]:not(:disabled)`) || datasetList.querySelector(`[data-id="${focusId}"] [data-action="select"]`))?.focus();
}
function announce(message) { document.querySelector('#dataset-announcement').textContent = message; }
function selectDataset(id) {
  selectedId = id; updatePills(); showControls();
}
function moveDataset(id, targetIndex) {
  const index = datasets.findIndex(item => item.id === id);
  if (index < 0 || index === targetIndex || targetIndex < 0 || targetIndex >= datasets.length) return;
  const [item] = datasets.splice(index, 1); datasets.splice(targetIndex, 0, item);
  updatePills(); schedule(); announce(`${item.settings.name}, layer ${targetIndex + 1} of ${datasets.length}.`);
}
datasetList.addEventListener('click', event => {
  const button = event.target.closest('button'); if (!button) return;
  const id = Number(button.closest('.dataset-pill').dataset.id);
  const index = datasets.findIndex(item => item.id === id), item = datasets[index];
  switch (button.dataset.action) {
    case 'select': selectDataset(id); break;
    case 'visible': item.settings.visible = !item.settings.visible; updatePills(); schedule(); break;
    case 'up': moveDataset(id, index - 1); break;
    case 'down': moveDataset(id, index + 1); break;
    case 'remove':
      datasets.splice(index, 1);
      removeFromWorker(id);
      if (selectedId === id) selectedId = datasets[Math.min(index, datasets.length - 1)]?.id ?? null;
      updatePills(); showControls(); schedule(); announce(`Removed ${item.settings.name}.`);
      datasetList.querySelector(`[data-id="${selectedId}"] [data-action="select"]`)?.focus();
      break;
  }
});
datasetList.addEventListener('dragstart', event => {
  if (!event.target.matches('.drag-handle')) return;
  draggedId = Number(event.target.closest('.dataset-pill').dataset.id);
  event.dataTransfer.setData('application/x-nmr-dataset', String(draggedId));
  event.dataTransfer.effectAllowed = 'move';
});
function clearDropMarkers() {
  for (const pill of datasetList.children) pill.classList.remove('drop-before', 'drop-after');
}
datasetList.addEventListener('dragover', event => {
  if (draggedId === null) return;
  event.preventDefault(); event.dataTransfer.dropEffect = 'move'; clearDropMarkers();
  const target = event.target.closest('.dataset-pill');
  if (target) target.classList.add(event.clientY < target.getBoundingClientRect().top + target.offsetHeight / 2 ? 'drop-before' : 'drop-after');
});
datasetList.addEventListener('drop', event => {
  if (draggedId === null) return;
  event.preventDefault(); event.stopPropagation();
  const target = event.target.closest('.dataset-pill');
  const id = draggedId; draggedId = null;
  if (target) {
    const from = datasets.findIndex(item => item.id === id);
    let to = datasets.findIndex(item => String(item.id) === target.dataset.id);
    if (event.clientY >= target.getBoundingClientRect().top + target.offsetHeight / 2) to++;
    if (from < to) to--;
    moveDataset(id, to);
  }
  clearDropMarkers(); updatePills();
});
datasetList.addEventListener('dragend', () => { draggedId = null; clearDropMarkers(); updatePills(); });
function showControls() {
  const item = selected();
  for (const name of entityFields) form.elements[name].disabled = !item;
  for (const name of ['negative_color', 'linewidth']) form.elements[name].closest('label').hidden = item?.kind !== 'spectrum';
  document.querySelector('#assignment-controls').hidden = item?.kind !== 'assignments';
  document.querySelector('#contour-controls').hidden = item?.kind !== 'spectrum';
  for (const name of ['loop', 'direct', 'indirect', 'labels', 'mirror']) form.elements[name].disabled = item?.kind !== 'assignments';
  for (const name of ['levels', 'base', 'count', 'factor', 'negative', 'negative_color', 'linewidth']) form.elements[name].disabled = item?.kind !== 'spectrum';
  if (item) {
    for (const name of entityFields) setField(name, item.settings[name]);
    if (item.kind === 'spectrum') {
      form.elements.base.disabled = item.settings.levels !== 'manual';
      form.elements.base.required = item.settings.levels === 'manual';
      for (const name of ['count', 'factor', 'negative']) form.elements[name].disabled = item.settings.levels === 'saved';
    }
  }
  colorPicker.update(item);
  updateAssignmentDownload();
}
function updateAssignmentDownload() {
  const item = selected();
  const csv = document.querySelector('#star-download');
  csv.hidden = !item?.result?.csv || item.key !== calculationKey(item);
  if (csvResult !== item?.result) {
    if (csvUrl) URL.revokeObjectURL(csvUrl);
    csvUrl = null; csvResult = item?.result;
  }
  if (!csv.hidden && !csvUrl) {
    csvUrl = URL.createObjectURL(new Blob([item.result.csv], {type: 'text/csv'}));
    csv.href = csvUrl;
    csv.download = `${item.file.name}.csv`;
  }
}
function scientificSettings(item) {
  return Object.fromEntries(calculationFields[item.kind].map(name => [name, item.settings[name]]));
}
function calculationKey(item) { return JSON.stringify(scientificSettings(item)); }
function needsCalculation(item) { return item.settings.visible && (!item.result || item.key !== calculationKey(item)); }
function queueDisplay() {
  if (frame !== null) return;
  frame = requestAnimationFrame(() => { frame = null; refreshDisplay().catch(showError); });
}
function schedule(delay = 0) {
  revision++;
  downloads.hidden = true;
  errors.hidden = true;
  if (ready) status.textContent = 'Updating plot…';
  queueDisplay();
  if (delay) {
    clearTimeout(timer);
    timer = setTimeout(() => { timer = null; calculate(); }, delay);
  } else if (timer === null) calculate();
}
function rpc(data) {
  return new Promise((resolve, reject) => {
    const id = ++request;
    pending.set(id, {resolve, reject});
    worker.postMessage({...data, request: id}, data.buffer ? [data.buffer] : []);
  });
}
function removeFromWorker(id) {
  const epoch = workerEpoch;
  if (ready) rpc({type: 'remove', id}).catch(error => { if (epoch === workerEpoch) showError(error); });
}
function startWorker() {
  const epoch = ++workerEpoch;
  worker?.terminate();
  for (const job of pending.values()) job.reject(Error('Calculation restarted.'));
  pending.clear();
  clearTimeout(timer); timer = null;
  for (const item of datasets) { item.loaded = false; item.key = null; }
  ready = false;
  revision++;
  errors.hidden = true;
  downloads.hidden = true;
  status.textContent = 'Loading plot…';
  worker = new Worker(new URL('./worker.js', import.meta.url), {type: 'module'});
  worker.onmessage = ({data}) => {
    if (epoch !== workerEpoch) return;
    if (data.type === 'ready') { ready = true; schedule(); }
    else if (data.type === 'fatal') { ready = false; showError(data.message); }
    else {
      const job = pending.get(data.request);
      pending.delete(data.request);
      if (data.type === 'error') job?.reject(Error(data.message));
      else job?.resolve(data.result);
    }
  };
  worker.onerror = event => {
    if (epoch !== workerEpoch) return;
    ready = false;
    for (const job of pending.values()) job.reject(Error(event.message));
    pending.clear();
    showError(event.message);
  };
}
async function addFiles(files) {
  for (const file of files) {
    if (file.size > 64 * 1024 ** 2) { showError(`${file.name}: file exceeds 64 MB.`); continue; }
    const head = new Uint8Array(await file.slice(0, 4096).arrayBuffer());
    const text = new TextDecoder().decode(head);
    const kind = text.startsWith('UCSF NMR') || (head[0] === 80 && head[1] === 75) ? 'spectrum' : 'assignments';
    if (kind === 'assignments' && (!/(?:^|\n)\s*data_/.test(text) || file.size > 16 * 1024 ** 2)) {
      showError(`${file.name}: choose UCSF, a Bruker ZIP with 2rr/procs/proc2s, or NMR-STAR (up to 16 MB).`);
      continue;
    }
    const item = {id: nextId++, file, kind, settings: {name: file.name, visible: true, levels: 'auto', base: '', count: '12', factor: '1.4', negative: false, colormap: 'solid', reverse_colormap: false, color: palette[datasets.length % palette.length], negative_color: '#b75442', linewidth: '.6', xoffset: '0', yoffset: '0', loop: '1', direct: 'CA', indirect: 'CB', labels: true, mirror: false}};
    datasets.push(item);
    selectedId = item.id;
    updatePills();
    showControls();
    schedule();
  }
}
async function exampleFile(url, name) {
  const response = await fetch(url);
  if (!response.ok) throw Error(`Example download failed: HTTP ${response.status}`);
  return new File([await response.arrayBuffer()], name);
}
async function calculate() {
  const epoch = workerEpoch;
  if (!ready || calculatingEpoch === epoch || !form.checkValidity()) return;
  calculatingEpoch = epoch;
  try {
    while (ready && epoch === workerEpoch && timer === null) {
      const item = datasets.find(needsCalculation);
      if (!item) break;
      const settings = scientificSettings(item), key = JSON.stringify(settings);
      try {
        if (!item.loaded) {
          status.textContent = `Reading ${item.settings.name}…`;
          const buffer = await item.file.arrayBuffer();
          if (epoch !== workerEpoch) return;
          if (!datasets.includes(item)) continue;
          await rpc({type: 'load', id: item.id, kind: item.kind, buffer, filename: item.file.name});
          if (epoch !== workerEpoch) return;
          item.loaded = true;
        }
        if (!datasets.includes(item)) continue;
        status.textContent = `Calculating ${item.settings.name}…`;
        const result = await rpc({type: 'calculate', id: item.id, settings, key});
        if (epoch !== workerEpoch) return;
        if (datasets.includes(item)) { item.result = result; item.key = key; queueDisplay(); }
      } catch (error) {
        if (epoch !== workerEpoch) return;
        if (datasets.includes(item) && item.settings.visible && key === calculationKey(item)) { showError(error); break; }
      }
    }
  } finally {
    if (epoch === workerEpoch) { calculatingEpoch = null; queueDisplay(); }
  }
}
function configuredView() { return Object.fromEntries(globalFields.map(name => [name, fieldValue(name)])); }
function manifest() {
  return {layerOrder: 'top to bottom', configuredView: configuredView(), view: plot.currentView(), datasets: datasets.map(item => ({name: item.file.name, kind: item.kind, settings: {...item.settings}, provenance: provenance(item)})), bokehjs: Bokeh.version};
}
function updateSettingsDownload() {
  if (settingsUrl) URL.revokeObjectURL(settingsUrl);
  settingsUrl = URL.createObjectURL(new Blob([JSON.stringify(manifest(), null, 2)], {type: 'application/json'}));
  Object.assign(document.querySelector('#settings-download'), {href: settingsUrl, download: 'nmr-overlay.json'});
}
async function refreshDisplay() {
  if (!form.checkValidity()) { downloads.hidden = true; status.textContent = 'Complete valid settings for the selected dataset.'; return; }
  if (!datasets.some(item => item.settings.visible)) {
    plot.clear(); downloads.hidden = true;
    document.querySelector('#stats').textContent = '';
    status.textContent = 'Add or show a dataset.';
    return;
  }
  if (!datasets.some(item => item.settings.visible && item.result)) return;
  await plot.initialize();
  // Read current state after initialization, which may overlap user edits.
  if (!datasets.some(item => item.settings.visible && item.result)) return;
  const settings = viewSettings(datasets, configuredView());
  plot.update(datasets, settings);
  document.querySelector('#spectrum-name').textContent = settings.title || 'NMR overlay';
  document.querySelector('#stats').textContent = datasets.filter(item => item.settings.visible && item.result).map(item => `${item.settings.name}: ${item.kind === 'spectrum' ? `${item.result.paths} contour paths` : `${item.result.count} assigned pairs`}`).join(' · ');
  const complete = ready && !datasets.some(needsCalculation) && errors.hidden;
  downloads.hidden = !complete;
  if (complete) {
    updateSettingsDownload();
    status.textContent = exporting ? 'Preparing SVG…' : 'Ready.';
  } else if (errors.hidden) status.textContent = ready ? 'Updating contours…' : 'Loading plot…';
  updateAssignmentDownload();
}
async function exportSvg() {
  if (exporting || downloads.hidden || !form.checkValidity()) return;
  const generation = revision, epoch = workerEpoch;
  const settings = plot.currentView();
  const items = datasets.filter(item => item.settings.visible).map(({id, key, settings}) => ({id, key, settings: {...settings}}));
  exporting = true; svgButton.disabled = true; status.textContent = 'Preparing SVG…';
  try {
    const result = await rpc({type: 'export', settings, datasets: items});
    if (epoch !== workerEpoch) return;
    if (generation !== revision || JSON.stringify(settings) !== JSON.stringify(plot.currentView())) { status.textContent = 'View changed. Select Download SVG again.'; return; }
    if (svgUrl) URL.revokeObjectURL(svgUrl);
    svgUrl = URL.createObjectURL(new Blob([result.svg], {type: 'image/svg+xml'}));
    const link = document.createElement('a');
    Object.assign(link, {href: svgUrl, download: 'nmr-overlay.svg'});
    document.body.append(link); link.click(); link.remove();
    status.textContent = 'Ready.';
  } catch (error) { if (epoch === workerEpoch && generation === revision) showError(error); }
  finally { exporting = false; svgButton.disabled = false; }
}
function settingsChanged(event) {
  if (entityFields.includes(event.target.name) && selected()) {
    const value = fieldValue(event.target.name);
    if (selected().settings[event.target.name] === value) return;
    selected().settings[event.target.name] = value;
    if (event.target.name === 'color') { selected().settings.colormap = 'solid'; setField('colormap', 'solid'); }
    colorPicker.update(selected());
    updatePills();
    if (['levels', 'colormap'].includes(event.target.name)) showControls();
  }
  if (entityFields.includes(event.target.name) || globalFields.includes(event.target.name)) schedule(selected() && calculationFields[selected().kind].includes(event.target.name) ? 350 : 0);
}
form.addEventListener('input', settingsChanged);
form.addEventListener('change', settingsChanged);
form.addEventListener('submit', event => event.preventDefault());
form.addEventListener('viewchange', () => { revision++; queueDisplay(); });
svgButton.addEventListener('click', exportSvg);
document.querySelector('#settings-download').addEventListener('click', updateSettingsDownload);
document.querySelector('#files').addEventListener('change', async event => { await addFiles(event.target.files); event.target.value = ''; });
const dropzone = document.querySelector('#dropzone');
for (const event of ['dragenter', 'dragover']) document.addEventListener(event, e => { if (e.dataTransfer.types.includes('Files')) { e.preventDefault(); dropzone.classList.add('dragging'); } });
document.addEventListener('dragleave', e => { if (!e.relatedTarget) dropzone.classList.remove('dragging'); });
document.addEventListener('drop', e => { if (e.dataTransfer.types.includes('Files')) { e.preventDefault(); dropzone.classList.remove('dragging'); addFiles(e.dataTransfer.files).catch(showError); } });
document.querySelector('#restart').addEventListener('click', startWorker);
reloadExample.addEventListener('click', event => {
  if (event.button !== 0 || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
  event.preventDefault(); loadMurrayExample();
});
for (const [name, value] of Object.entries(exampleView)) setField(name, value);
showControls(); startWorker();
loadMurrayExample();

async function loadMurrayExample() {
  if (loadingExamples) return;
  loadingExamples = true;
  reloadExample.setAttribute('aria-busy', 'true');
  reloadExample.setAttribute('aria-disabled', 'true');
  try {
    const files = await Promise.all([
      exampleFile('./examples/CC_2D-U-FUS.ucsf', 'Uniformly labeled FUS-LC (Murray et al.).ucsf'),
      exampleFile('./examples/CC-2D-N112-FUS.ucsf', 'N112 FUS-LC S112C (Murray et al.).ucsf'),
      exampleFile('./examples/bmr30304.str', 'BMRB 30304 assignments.str'),
    ]);
    for (const item of datasets) removeFromWorker(item.id);
    datasets = []; selectedId = null; draggedId = null;
    clearTimeout(timer); timer = null;
    for (const [name, value] of Object.entries(exampleView)) setField(name, value);
    plot.reset(); updatePills(); showControls(); schedule();
    await addFiles(files);
    moveDataset(datasets.find(item => item.kind === 'assignments').id, 0);
  } catch (error) { showError(error); }
  finally {
    loadingExamples = false;
    reloadExample.removeAttribute('aria-busy');
    reloadExample.removeAttribute('aria-disabled');
  }
}
