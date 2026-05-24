from . import analysis, data, plot, model, tools, inversion
import importlib.metadata

try:
    __version__ = importlib.metadata.version("pyissm")
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"