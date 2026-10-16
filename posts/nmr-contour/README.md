# NMR Contour

A local browser app for processed 2D NMR spectra, with an optional Flask interface. nmrglue reads spectra; Bokeh provides the interactive view; Matplotlib produces editable SVG from the same contour coordinates. No database or frontend build. Copy this directory to use it as a separate project.

## Run in the browser

From the website repository root, start a static file server:

```sh
python3 -m http.server 8766 --bind 127.0.0.1
```

Open **http://127.0.0.1:8766/posts/nmr-contour/**. Serve the repository root so the page can load the shared site stylesheet and fonts from `css/`.

The page opens with the Murray et al. uniformly labeled and N112 FUS-LC spectra plus BMRB 30304 assignment markers, titled “FUS-LC: uniformly labeled + N112 (Murray et al.),” with 10–80 ppm on both axes and equal ppm scaling. Drop multiple files anywhere on the page, or use the file chooser. Each UCSF, Bruker ZIP, or NMR-STAR file becomes a separate dataset. Click a file pill to edit its name, color, reference offsets, and contours or assignment markers. Each compact single-row pill has visibility and remove controls; full filenames are available in the button tooltip and editable name field. Drag its handle, or use the up/down buttons, to reorder layers; the top pill is drawn on top in both Bokeh and SVG, including assignment markers and labels. Keyboard users can activate the same buttons with Enter or Space. Remove deletes only the selected dataset from the browser session. The **Reload Murray et al. 2017 example data** link replaces the datasets with both spectra and the assignments, and restores the default view. The other local test spectrum is `../../assets/data/FUS-LC-2D-INEPT.ucsf`.

The unified Color picker offers ColorBrewer Paired, Set1, Set2, and Set3 as individual fixed swatches, a custom native color picker plus hex input, and sequential gradients: Blues, Greens, Purples, Oranges, BuPu, PuBu, GnBu, BuGn, YlGn, YlGnBu, YlOrBr, YlOrRd, PuBuGn, viridis, plasma, inferno, magma, and cividis. Categorical ColorBrewer palettes supply solid swatches only; the contour gradients remain sequential. No diverging gradients are included. Positive contour levels sample the full map evenly in ascending level order; reverse is optional. Negative contours retain their separate color. The browser plot and SVG use the same Matplotlib-generated per-path colors, and JSON saves the palette choice and direction. Assignment markers retain solid colors; the gradient tiles stay visible but disabled, with a note to select a spectrum.

Axis limits, axis labels, figure title, export dimensions, font size, and equal scaling apply to the shared figure. Clear limits to fit visible datasets. All visible datasets are included in the SVG; JSON records each dataset's settings and source digest in top-to-bottom layer order. Assignment CSV is available for the selected assignment dataset. Contour geometry is cached per dataset, so edits to colors or the shared view do not recalculate contours. Each file is independent; there is no automatic referencing or pairing with another dataset.

For Bruker ZIPs, `2rr`, `procs`, and `proc2s` must be together in one directory; `clevels` is optional. A nested experiment folder is supported if it contains exactly one `2rr`. A standalone `2rr` or raw FID is not accepted because the processing metadata is needed for shape, scaling, and ppm referencing. This tool supports processed 2D spectra, not 1D `1r`.

The plot updates after a 350 ms pause when a file or setting changes. Pan or zoom updates the exported SVG bounds without changing the axis input fields. Bokeh Reset returns to the entered limits. Editing those limits updates the reset target. JSON records both configuredView and the exported view. Downloads appear only after the current calculation finishes.

This server only serves files. Spectrum parsing, contour calculation, and SVG generation run inside a browser Web Worker. Selected spectra are not uploaded. The first load needs internet access to download Pyodide 314.0.7 and scientific packages from jsDelivr/PyPI. `nmrglue` is pinned to 0.11; NumPy, SciPy, and Matplotlib use the versions in that Pyodide release. This is not a guaranteed offline installation.

The browser version supports the FUS LC DARR example, UCSF, processed Bruker ZIP, contour controls, axis referencing, interactive viewing, and SVG/settings downloads. NMR-STAR parsing and overlay export run in the same worker. **Restart calculation** terminates the worker and loads a new runtime if a calculation must be cancelled. Large inputs can use substantially more memory than their file size. Width, height, and equal-ppm scaling control the exported SVG; the interactive view fits the available screen.

The browser HTML, CSS, JavaScript, Python modules, `static/`, `vendor/`, and `examples/` can be served by GitHub Pages or another static host alongside the website’s shared `css/` directory. No Flask server is needed for this version. BokehJS 3.9.2 is included locally with its license. The app and bundled examples are eligible for Git tracking; `.venv/`, `__pycache__/`, and compiled Python files remain ignored.

## NMR-STAR assignments in the browser

Drop a `.str`, `.star`, or `.txt` NMR-STAR file, click its file pill, choose the shift list by file order, and enter two distinct atom names (for example CA/CB or N/H). Markers, optional residue labels, and CSV update automatically. Files up to 16 MB are accepted. Mirror across the diagonal is optional and intended only for spectra with identical nuclei on both axes; it is off by default.

These are paired assignments within a residue, not detected peaks. Assembly, entity, residue index, and residue identity define each group. Missing shifts are omitted; incomplete identity, duplicate selected-atom assignments, nonfinite shifts, and invalid list numbers produce errors. Labels include residue name/index plus assembly/entity IDs. All rows are plotted and exported.

Startup and **Reload Murray et al. 2017 example data** load the unmodified [BMRB 30304](https://bmrb.io/data_library/summary/?bmrbId=30304) NMR-STAR deposition from Murray et al., Cell (2017), associated with PDB 5W3N. This supplies 35 intraresidue CA/CB pairs, shown above the spectra with residue labels. The assignments are from the same study as the default Murray pair. They are not verified assignments for the separate, optional FUS LC CP-DARR example. Labels preserve deposited Comp_index_ID numbering; do not infer full-length FUS residue numbers. Provenance is recorded in `examples/README.md`.

The browser uses the unmodified Python parser from MIT-licensed `nmrstarlib` 2.1.1, bundled under `vendor/`. Native `pynmrstar` uses a compiled extension that is not directly available in Pyodide. Browser CSV output was checked against native `pynmrstar` 3.6.2 on a two-entity fixture and public BMRB entry 15000. No file or chemical-shift content is sent to a server.

## Run the Flask version

Python 3.11+ and [uv](https://docs.astral.sh/uv/) are required.

```sh
uv sync --locked
uv run python app.py
```

Open **http://127.0.0.1:5050**. Stop with Ctrl+C. This version requires a Python server and cannot run directly on GitHub Pages. It retains the NMR-STAR assignment tools. Both interfaces use `spectrum.py` for spectrum reading, contour calculation, and SVG generation.

## Included

- Processed 2D UCSF/SPARKY upload.
- Bruker ZIP upload containing exactly one `2rr`, with `procs`, `proc2s`, and optional `clevels` in the same folder. Raw FIDs are not needed.
- Automatic, manual geometric, or saved contour levels; positive and negative contours.
- Spectrum size, intensity extrema, and full ppm ranges.
- Axis limits, additive reference offsets, labels, equal ppm scaling, figure size, line width, font size, and color.
- SVG export with one editable path per continuous contour and editable text. All calculated contour vertices are retained without path simplification; closed loops are explicitly closed and strokes use round joins. JSON settings include input SHA-256, actual levels, axes, and library versions. Settings import is not implemented.
- Uploaded NMR-STAR shift lists to intraresidue atom-pair CSV tables. Choose the shift list by file order and atom names (for example CA/CB, C/CA, or N/H). Entities and assemblies stay separate; duplicate atom assignments fail explicitly. These are candidate correlations, not observed peaks.
- Murray uniformly labeled + N112 FUS-LC as the starting pair, with FUS LC DARR available as another example. The legacy Flask interface retains its synthetic demo.

In the Flask version, click **Update plot** after edits. This resubmits the selected file and settings. Uploads are handled in temporary storage, deleted after the request; previews and downloads are returned in the response. There are no accounts or stored sessions.

## Scientific conventions and limits

X is the direct dimension (array columns); Y is the indirect dimension (rows). Both axes decrease from left/bottom to right/top. UCSF axes use nmrglue's `make_uc`. Each Bruker axis uses its own processing parameters:

```text
ppm[i] = OFFSET - i × SW_p / (SF × SI)
```

This draft requires `SI` to match the loaded dimension length. Cropped or unusual Bruker processing metadata may need UCSF export. Check referencing against the source software before publication. Axis offsets are applied before plot limits.

Automatic levels use `median(data) + 1.5 × std(data)` as a display heuristic, following the Bruker helper in nmrlib. This is not a validated noise estimate. Manual levels use `base × factor**arange(count)`. Saved `clevels` retain signs. There is no interpolation, normalization, smoothing, or peak picking. Negative contours default to rust red and can be changed per dataset.

The browser supports multiple spectra and assignment overlays. The legacy Flask view remains single-spectrum. 1D spectra and remote BMRB lookup are not implemented. Pan, wheel zoom, box zoom, and reset are available in the plot toolbar. Axis fields remain fixed during pan and zoom; Reset restores their configured limits. The browser version updates SVG bounds automatically; in Flask, click Update plot. Bokeh resources are embedded locally, so no CDN or Bokeh server is needed. Maximum request size is 64 MB, unpacked ZIP size 256 MB, and loaded spectrum size 16 million points. Invalid input produces an error; parser exceptions remain visible in the server terminal.

## Browser validation

Checked in headless desktop Chrome on 2026-10-15 using the actual Pyodide worker:

| Input | Shape | Contour paths | Browser calculation time |
| --- | --- | --- | --- |
| FUS-LC-2D-INEPT.ucsf | 1024 × 2048 | 132 | 0.77 s |
| FUS-LC-DARR.ucsf | 1024 × 4096 | 357 | 1.44 s |

These runs used automatic levels with negative contours enabled. Times exclude runtime/package download and browser drawing, and depend on the computer. Native and browser results matched source SHA-256, shapes, intensity extrema, ppm endpoints, contour counts, and levels (relative tolerance 1e-12). All corresponding contour coordinates agreed within 2.85e-14 ppm. The shared reader also returned exactly the same arrays and axes as the original Flask implementation for both spectra. Both SVG downloads were exercised in Chrome; the desktop page was visually inspected.

Manual positive-only contours, zoom-to-export bounds, processed Bruker ZIP with automatic/saved levels, mobile layout after resize, and invalid-file error reporting with hidden stale downloads were checked. This verifies implementation agreement for these inputs, not the correctness of their original acquisition or referencing metadata.

## Earlier Flask validation

Checked on 2026-10-15:

- Synthetic preview, negative/manual contours, and rejection of missing files, fractional level counts, missing manual base, and nonoverlapping limits.
- Both real 1024 × 4096 spectra in the blog's `posts/assets/example_data.zip` parsed and returned SVG output.
- A small NMR-STAR fixture kept two entities with the same residue index separate in the crosspeak CSV.
- Local HTTP response returned 200. No browser UI inspection was possible in this environment.

The blog's UCSF example reports Y endpoints approximately −35.600 to −247.740 ppm, while the Bruker example reports 221.176 to 9.011 ppm. This is unresolved source referencing, not evidence of equivalence between these files. No automatic correction is applied. These checks validate the draft's I/O, not scientific accuracy for every acquisition or processing setup.

## Sources

- [Plot 2D NMR Spectra in the Browser](https://ray-berkeley.github.io/posts/Plot_2D_NMR_Spectra_in_Browser.html)
- [nmrlib](https://github.com/ray-berkeley/nmrlib), inspected at `9a90e8acd44313d56dcdb5144fb45d07f81cec18`.
- [nmrglue processed Bruker reader](https://nmrglue.readthedocs.io/en/latest/reference/generated/nmrglue.fileio.bruker.read_pdata.html)

This draft implements the workflows directly with current underlying libraries. It does not import nmrlib, whose package currently pins Matplotlib to <=3.2.2. The Flask assignment tools use pynmrstar; the browser uses vendored nmrstarlib.

## Overlay validation (2026-10-15)

Headless Chrome exercised both local FUS spectra in one plot, changed only the second dataset to five positive/five negative levels, and confirmed the first dataset's levels did not change. Added BMRB 30304, checked labels in the combined SVG, hid and removed it, and then added the STAR file through an actual drag-and-drop event. Mirroring yielded 70 markers from 35 pairs; a +1 ppm X offset placed the first marker at 58.6 ppm. Desktop and 390 px mobile layouts were captured and inspected; no browser errors occurred. Browser CSV matched native parsing, and native pairs matched pynmrstar 3.6.2 independently. These checks do not establish that the reference assignments match the example spectrum.


The **Reload Murray et al. 2017 example data** link loads the unmodified `examples/CC_2D-U-FUS.ucsf` and `examples/CC-2D-N112-FUS.ucsf` spectra from [Mendeley Data](https://data.mendeley.com/datasets/rxbh442x2g/1), with attribution under CC BY 4.0, plus BMRB 30304 assignments. These are also loaded on startup in distinct colors, with 10–80 ppm on both axes. The N112 sample labels residues 1–111 and has the S112C substitution. Assignment markers and labels can be hidden or removed through their file pill. New spectra default to positive contours only; enable negative contours in Contours when needed. Saved Bruker clevels retain their deposited signs.
