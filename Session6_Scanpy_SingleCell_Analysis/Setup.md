# Tutorial Setup Guide

Follow these steps to get your Python environment ready for the tutorial.  
This will install the [`uv`](https://docs.astral.sh/uv/) package manager, create a compute environment with
the required tools, and launch Jupyter Lab.


## 1. Install uv

### macOS / Linux
Run this command in your terminal:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Make sure `uv` is available in all future shells by adding it to your PATH:

```bash
printf '\nexport PATH="%s:$PATH"\n' "$(dirname $(which uv))" >> "$HOME/.profile"
```

### Windows (PowerShell)
Run this instead:

```powershell
iwr https://astral.sh/uv/install.ps1 -UseBasicParsing | iex
```

After that, you can follow the same steps as macOS/Linux.


## 2. Create the Python environment

Inside the tutorial folder, create a new virtual environment:

```bash
uv venv
```

Then install the packages we’ll use:

```bash
uv pip install scanpy uv ipykernel palantir igraph leidenalg rpy2 scikit-image anndata2ri
```

Register this environment as a Jupyter kernel:

```bash
uv run python -m ipykernel install --user --name scanpy --display-name "scanpy env"
```


## 3. Launch Jupyter Lab

Start a [Jupyter Lab](https://jupyterlab.readthedocs.io/en/latest/index.html) instance (from the same directory):

```bash
uvx jupyter lab
```

This launches Jupyter Lab in a temporary environment managed by uv. It’s a quick way to get a notebook server without cluttering your project venv. If you prefer a dedicated setup, install Jupyter Lab inside its own virtual environment and run it from there.

### Expected Output

When Jupyter Lab starts, you should see something like this in the terminal:

```
[I 2025-01-15 10:20:44.123 ServerApp] Jupyter Server 2.14.2 is running at:
[I 2025-01-15 10:20:44.123 ServerApp] http://localhost:8888/lab?token=123abc456def...
[I 2025-01-15 10:20:44.123 ServerApp]  or http://127.0.0.1:8888/lab?token=123abc456def...
[I 2025-01-15 10:20:44.123 ServerApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
```

Usually your browser will open automatically.  
If it does not, **copy the full URL (including the `?token=...` part)** into your browser.


Now you are ready: choose the **"scanpy env"** kernel inside Jupyter Lab and start working with the tutorial notebooks!
