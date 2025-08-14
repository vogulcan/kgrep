# kgrep/__init__.py
from __future__ import annotations

# Expose a stable version string if available (nice for --version flags, logs)
try:
    from importlib.metadata import version as _pkg_version

    __version__ = _pkg_version("kgrep")
except Exception:
    __version__ = "0.0.0"


def main() -> None:
    """
    Console entrypoint used by [project.scripts].
    Thin wrapper that calls the real CLI in kgrep.kgrep:main.
    """
    from .kgrep import main as _main

    _main()


__all__ = ["main", "__version__"]
