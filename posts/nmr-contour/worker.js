import {loadPyodide} from 'https://cdn.jsdelivr.net/pyodide/v314.0.7/full/pyodide.mjs';
let python;
async function fetchFile(filename) {
  const response = await fetch(`./${filename}`);
  if (!response.ok) throw Error(`${filename}: HTTP ${response.status}`);
  return [filename, new Uint8Array(await response.arrayBuffer())];
}
async function initialize() {
  const [runtime, files] = await Promise.all([loadPyodide(), Promise.all(['spectrum.py', 'assignments.py', 'overlay.py', 'dataset_store.py', 'vendor/nmrstarlib-2.1.1.zip'].map(fetchFile))]);
  python = runtime;
  await python.loadPackage(['numpy', 'scipy', 'matplotlib', 'micropip']);
  await python.runPythonAsync("import micropip\nawait micropip.install('nmrglue==0.11')");
  for (const [filename, bytes] of files) {
    if (filename.endsWith('.zip')) python.unpackArchive(bytes, 'zip', {extractDir: '/home/pyodide'});
    else python.FS.writeFile(`/home/pyodide/${filename}`, bytes);
  }
  await python.runPythonAsync('import json\nfrom dataset_store import DatasetStore\nstore = DatasetStore()');
  self.postMessage({type: 'ready'});
}
async function handle(data) {
  let payload;
  try {
    const {buffer, ...request} = data;
    python.globals.set('request_json', JSON.stringify(request));
    if (buffer) {
      payload = python.toPy(new Uint8Array(buffer));
      python.globals.set('payload', payload);
    }
    const result = JSON.parse(await python.runPythonAsync(`json.dumps(store.dispatch(json.loads(request_json)${buffer ? ', bytes(payload)' : ''}), allow_nan=False)`));
    self.postMessage({type: 'result', request: data.request, result});
  } catch (error) {
    self.postMessage({type: 'error', request: data.request, message: String(error)});
  } finally {
    payload?.destroy();
    for (const key of ['payload', 'request_json']) if (python.globals.has(key)) python.globals.delete(key);
  }
}
let queue = Promise.resolve();
self.onmessage = ({data}) => { queue = queue.then(() => handle(data)); };
initialize().catch(error => self.postMessage({type: 'fatal', message: String(error)}));
