"""Input, output, and runtime-environment helpers."""

import os
from datetime import datetime
from pathlib import Path


def configure_runtime_environment(project_dir, include_celloracle_cache=True):
    """Configure writable cache directories before scientific packages are imported."""
    project_dir = Path(project_dir)
    relative_paths = {
        "MPLCONFIGDIR": "results/matplotlib_cache",
        "XDG_CACHE_HOME": "results/xdg_cache",
    }
    if include_celloracle_cache:
        relative_paths.update(
            {
                "XDG_CONFIG_HOME": "results/xdg_config",
                "NUMBA_CACHE_DIR": "results/celloracle_numba_cache",
            }
        )
    for name, relative in relative_paths.items():
        value = str(project_dir / relative)
        os.environ.setdefault(name, value)
        Path(os.environ[name]).mkdir(parents=True, exist_ok=True)


def unique_output_dir(base_dir):
    """Return the requested directory, or a timestamped sibling if it exists."""
    base = Path(base_dir)
    if not base.exists():
        return base
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return base.with_name(f"{base.name}_{stamp}")


def make_timestamped_output_dir(prefix, children):
    """Create a new timestamped output directory and its required children."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = Path(f"{prefix}_{stamp}")
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {out}")
    for child in children:
        (out / child).mkdir(parents=True, exist_ok=False)
    return out
