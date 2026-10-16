import {fixedPalettes, gradientPalettes} from './color-palettes.js';

export function colorBackground(settings) {
  if (settings.colormap === 'solid') return settings.color;
  const colors = gradientPalettes[settings.colormap];
  return `linear-gradient(to right, ${(settings.reverse_colormap ? [...colors].reverse() : colors).join(',')})`;
}

export class ColorPicker {
  constructor(form, choose) {
    this.form = form;
    this.choose = choose;
    this.element = document.querySelector('#color-picker');
    this.fixed = document.querySelector('#fixed-swatches');
    this.gradients = document.querySelector('#gradient-swatches');
    this.palette = document.querySelector('#fixed-palette');
    this.hex = document.querySelector('#custom-hex');
    for (const name of Object.keys(fixedPalettes)) this.palette.add(new Option(name, name));
    this.palette.addEventListener('change', () => { this.buildFixed(); this.update(this.item); });
    this.buildFixed();
    for (const [name, colors] of Object.entries(gradientPalettes)) {
      const button = this.button(name, `Use ${name} gradient`);
      button.dataset.colormap = name;
      const preview = document.createElement('span'); preview.className = 'gradient-tile';
      preview.style.background = `linear-gradient(to right, ${colors.join(',')})`;
      preview.setAttribute('aria-hidden', 'true');
      button.prepend(preview);
      button.addEventListener('click', () => choose({colormap: name}));
      this.gradients.append(button);
    }
    this.hex.addEventListener('input', () => {
      const value = this.hex.value.trim();
      const valid = /^#[0-9a-f]{6}$/i.test(value);
      this.hex.setAttribute('aria-invalid', String(!valid));
      if (valid && this.item) choose({color: value.toLowerCase(), colormap: 'solid'});
    });
    this.hex.addEventListener('change', () => {
      if (this.hex.getAttribute('aria-invalid') === 'true') this.hex.value = this.item?.settings.color || '';
      this.hex.setAttribute('aria-invalid', 'false');
    });
    this.element.addEventListener('keydown', event => {
      if (event.key === 'Escape') { this.element.open = false; this.element.querySelector('summary').focus(); }
    });
    document.addEventListener('click', event => { if (!this.element.contains(event.target)) this.element.open = false; });
  }
  button(text, label) {
    const button = document.createElement('button'); button.type = 'button';
    button.textContent = text; button.title = label; button.setAttribute('aria-label', label);
    return button;
  }
  buildFixed() {
    this.fixed.replaceChildren(...fixedPalettes[this.palette.value].map((color, index) => {
      const button = this.button('', `${this.palette.value} ${index + 1}: ${color}`);
      button.style.background = color; button.dataset.color = color;
      button.addEventListener('click', () => this.choose({color, colormap: 'solid'}));
      return button;
    }));
  }
  update(item) {
    this.item = item;
    const spectrum = item?.kind === 'spectrum';
    const settings = item?.settings;
    const gradient = spectrum && settings.colormap !== 'solid';
    const title = settings ? (gradient ? settings.colormap + (settings.reverse_colormap ? ' (reversed)' : '') : settings.color.toUpperCase()) : 'Select a file';
    document.querySelector('#color-current-label').textContent = title;
    document.querySelector('#color-current-swatch').style.background = settings ? colorBackground(settings) : '';
    this.palette.disabled = !item; this.hex.disabled = !item;
    if (document.activeElement !== this.hex) this.hex.value = settings?.color || '';
    for (const button of this.fixed.children) {
      button.disabled = !item;
      button.setAttribute('aria-pressed', String(Boolean(settings && !gradient && settings.color === button.dataset.color)));
    }
    for (const button of this.gradients.children) {
      button.disabled = !spectrum;
      button.setAttribute('aria-pressed', String(Boolean(gradient && settings.colormap === button.dataset.colormap)));
    }
    document.querySelector('#gradient-note').textContent = spectrum ? 'Low → high positive contour levels' : 'Select a spectrum to use contour gradients.';
    this.form.elements.reverse_colormap.disabled = !gradient;
    document.querySelector('#reverse-gradient').hidden = !gradient;
  }
}
