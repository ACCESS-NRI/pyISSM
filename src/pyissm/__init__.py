<<<<<<< HEAD
from . import analysis, data, plot, model, tools
=======
from . import analysis, data, learn, plot, model, tools
>>>>>>> 3d54be9 (Initial documentation Infrastructure (#88))
import importlib.metadata

try:
    __version__ = importlib.metadata.version("pyissm")
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"