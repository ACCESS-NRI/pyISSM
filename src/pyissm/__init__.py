from . import analysis, data, learn, plot, model, tools
import importlib.metadata

try:
    __version__ = importlib.metadata.version("pyissm")
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"