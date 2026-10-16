from assignments import crosspeaks
from overlay import compose
from spectrum import SpectrumRenderer, require


class DatasetStore:
    def __init__(self):
        self.entries = {}

    def dispatch(self, request: dict, payload: bytes | None = None) -> dict:
        action = request['type']
        if action == 'load':
            identity, kind = request['id'], request['kind']
            require(identity not in self.entries, 'Dataset already loaded.')
            require(kind in ('spectrum', 'assignments'), 'Unknown dataset type.')
            source = payload
            if kind == 'spectrum':
                source = SpectrumRenderer({'source': 'auto'}, payload, request['filename'])
                source.data
            self.entries[identity] = dict(kind=kind, source=source, key=None, result=None)
            return {}
        if action == 'remove':
            # Removal can arrive while the browser is still reading the file.
            self.entries.pop(request['id'], None)
            return {}
        if action == 'calculate':
            entry = self.entries[request['id']]
            if entry['key'] != request['key']:
                settings = request['settings']
                if entry['kind'] == 'spectrum':
                    spectrum = entry['source']
                    spectrum.settings = dict(settings, source=spectrum.settings['source'])
                    result = spectrum.render(export_svg=False)
                else:
                    result = crosspeaks(entry['source'], settings)
                entry.update(key=request['key'], result=result)
            return entry['result']
        if action == 'export':
            datasets = []
            for item in request['datasets']:
                entry = self.entries[item['id']]
                require(entry['key'] == item['key'] and entry['result'] is not None, 'Dataset changed; wait for the plot, then export again.')
                datasets.append(dict(item, kind=entry['kind'], result=entry['result']))
            return compose(datasets, request['settings'])
        raise ValueError(f'Unknown worker action: {action}')
