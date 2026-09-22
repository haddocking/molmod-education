# JupyterLite site

## Add a notebook

1. Put the `.ipynb` in `Notebooks/`
2. Build picks it up automatically
3. Check deps: pyodide kernel only installs pure-Python / pyodide wheels, maybe you dep is not compatible
4. Build + serve locally, open it, confirm it runs under pyodide

## Build the jupyterlite locally for testing

```bash
cd jupyterlite
# setup the venv for jupyterlite
python -m venv .venv-lite
source .venv-lite/bin/activate
# install the jupyterlite deps
pip install -r requirements.txt
# build the jupyterlite site
jupyter lite build --contents ../Notebooks --output-dir dist --lite-dir .
```

Output will be in `./dist`

## Serve

```bash
python3 -m http.server -d dist 8000
```

Go to <http://localhost:8000> and check

