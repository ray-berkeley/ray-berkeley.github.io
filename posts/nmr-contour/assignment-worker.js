import {loadPyodide} from 'https://cdn.jsdelivr.net/pyodide/v314.0.7/full/pyodide.mjs';
let python;
async function initialize() {
  python = await loadPyodide();
  const archive = await fetch('./vendor/nmrstarlib-2.1.1.zip');
  if (!archive.ok) throw Error(`Cannot load NMR-STAR parser: HTTP ${archive.status}`);
  python.unpackArchive(await archive.arrayBuffer(), 'zip', {extractDir: '/home/pyodide'});
  const source = await fetch('./assignments.py');
  if (!source.ok) throw Error(`Cannot load assignments.py: HTTP ${source.status}`);
  python.FS.writeFile('/home/pyodide/assignments.py', await source.text());
  await python.runPythonAsync('import json\nfrom assignments import crosspeaks');
  self.postMessage({type: 'ready'});
}
self.onmessage = async ({data}) => {
  let payload;
  try {
    payload = python.toPy(new Uint8Array(data.buffer));
    python.globals.set('payload', payload);
    python.globals.set('settings_json', JSON.stringify(data.settings));
    const result = JSON.parse(await python.runPythonAsync('json.dumps(crosspeaks(bytes(payload), json.loads(settings_json)), allow_nan=False)'));
    self.postMessage({type: 'result', revision: data.revision, result});
  } catch (error) {
    self.postMessage({type: 'error', revision: data.revision, message: String(error)});
  } finally {
    payload?.destroy();
    python.globals.delete('payload');
    python.globals.delete('settings_json');
  }
};
initialize().catch(error => self.postMessage({type: 'fatal', message: String(error)}));
