const form = document.querySelector('#assignment-form');
const status = document.querySelector('#star-status');
const errorBox = document.querySelector('#star-error');
const results = document.querySelector('#star-results');
let worker;
let ready = false;
let busy = false;
let revision = 0;
let timer;
let downloadURL;

function fail(message) {
  console.error(message);
  errorBox.textContent = String(message).trim().split('\n').at(-1);
  errorBox.hidden = false;
  results.hidden = true;
  status.textContent = 'Assignment parsing failed. Correct the input and try again.';
}
function initialize() {
  worker?.terminate();
  ready = false;
  busy = false;
  results.hidden = true;
  errorBox.hidden = true;
  status.textContent = 'Loading the assignment reader…';
  worker = new Worker(new URL('./assignment-worker.js', import.meta.url), {type: 'module'});
  worker.onerror = event => {ready = false; fail(event.message);};
  worker.onmessage = ({data}) => {
    if (data.type === 'ready') {ready = true; schedule();}
    if (data.type === 'fatal') {ready = false; fail(data.message);}
    if (data.type === 'result' || data.type === 'error') {
      busy = false;
      if (data.revision !== revision) {clearTimeout(timer); timer = setTimeout(calculate, 350); return;}
      if (data.type === 'error') fail(data.message);
      else display(data.result);
    }
  };
}
function display(result) {
  const table = document.querySelector('#star-table');
  table.replaceChildren();
  const head = table.createTHead().insertRow();
  for (const label of result.columns) {
    const cell = document.createElement('th');
    cell.scope = 'col';
    cell.textContent = label;
    head.append(cell);
  }
  const body = table.createTBody();
  for (const record of result.rows) {
    const row = body.insertRow();
    for (const value of record) row.insertCell().textContent = value;
  }
  if (downloadURL) URL.revokeObjectURL(downloadURL);
  downloadURL = URL.createObjectURL(new Blob([result.csv], {type: 'text/csv;charset=utf-8'}));
  document.querySelector('#star-download').href = downloadURL;
  document.querySelector('#star-summary').textContent = `${result.count} candidate correlations · ${result.loops} shift lists in file.`;
  document.querySelector('#star-truncated').hidden = result.count <= 100;
  results.hidden = false;
  status.textContent = 'Ready. Assignments processed in this browser.';
}
function schedule() {
  revision++;
  clearTimeout(timer);
  results.hidden = true;
  errorBox.hidden = true;
  if (ready) status.textContent = 'Waiting for assignment settings…';
  timer = setTimeout(calculate, 350);
}
async function calculate() {
  if (busy) return;
  if (!form.checkValidity()) {status.textContent = 'Choose a file and complete the assignment settings.'; return;}
  const file = form.elements.star.files[0];
  if (file.size > 16 * 1024**2) {fail('NMR-STAR file exceeds 16 MB.'); return;}
  if (!worker) {initialize(); return;}
  if (!ready) return;
  busy = true;
  const target = worker;
  const requested = revision;
  status.textContent = 'Reading assignments…';
  try {
    const settings = Object.fromEntries([...new FormData(form)].filter(([key]) => key !== 'star'));
    const buffer = await file.arrayBuffer();
    if (target !== worker) return;
    worker.postMessage({revision: requested, settings, buffer}, [buffer]);
  } catch (error) {
    if (target === worker) {busy = false; fail(String(error));}
  }
}
form.addEventListener('input', schedule);
form.addEventListener('change', schedule);
form.addEventListener('submit', event => {event.preventDefault(); schedule();});
document.querySelector('#star-restart').addEventListener('click', () => {revision++; initialize();});
