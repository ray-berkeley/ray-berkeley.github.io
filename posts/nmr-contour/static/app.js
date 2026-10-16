document.addEventListener('htmx:beforeSwap', (event) => {
  if (event.detail.xhr.status >= 400 && event.detail.xhr.status < 500) {
    event.detail.shouldSwap = true;
    event.detail.isError = false;
  }
  if (event.detail.xhr.status >= 500) {
    event.detail.target.innerHTML = '<p class="error" role="alert">The file could not be processed. See the server terminal for the error details, then correct the file and try again.</p>';
  }
});
document.addEventListener('htmx:sendError', (event) => {
  event.detail.target.innerHTML = '<p class="error" role="alert">Cannot reach the local server. Start it and try again.</p>';
});
const source = document.querySelector('[name=source]');
const file = document.querySelector('[name=spectrum]');
function updateSource() {
  file.disabled = source.value === 'demo';
  file.required = source.value !== 'demo';
  file.accept = source.value === 'bruker' ? '.zip' : '.ucsf';
}
source.addEventListener('change', updateSource);
updateSource();
const mode = document.querySelector('[name=levels]');
const base = document.querySelector('[name=base]');
function updateLevels() {
  base.disabled = mode.value !== 'manual';
  base.required = mode.value === 'manual';
  for (const name of ['count', 'factor', 'negative']) {
    document.querySelector(`[name=${name}]`).disabled = mode.value === 'saved';
  }
}
mode.addEventListener('change', updateLevels);
updateLevels();
