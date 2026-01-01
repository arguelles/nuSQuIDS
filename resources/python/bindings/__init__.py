import sys as _sys
import importlib as _importlib

# Import the extension module
_ext_module = _importlib.import_module('.nuSQuIDS', __name__)

# Export all public symbols from compiled extension
for _name in dir(_ext_module):
    if not _name.startswith('_'):
        globals()[_name] = getattr(_ext_module, _name)

# Import tools
from .nuSQUIDSTools import *

# The C++ class is nuSQUIDS (all caps UIDS), matching the Python binding.
# The compiled extension file is nuSQuIDS.so (mixed case).
# Provide an alias so users can use either capitalization:
#   nsq.nuSQUIDS(...) - matches C++ class name (recommended)
#   nsq.nuSQuIDS(...) - matches module file name (also works)

# Store references to the classes
_nuSQUIDS_class = _ext_module.nuSQUIDS
_nuSQUIDSAtm_class = _ext_module.nuSQUIDSAtm

# Set the class aliases
nuSQuIDS = _nuSQUIDS_class
nuSQuIDSAtm = _nuSQUIDSAtm_class

# Rename the submodule in sys.modules so Python doesn't auto-set nuSQuIDS
# (Python auto-adds submodules to parent package namespace, which would override our alias)
# The extension is still accessible as nusquids._nuSQuIDS_ext if needed
if (__name__ + '.nuSQuIDS') in _sys.modules:
    _sys.modules[__name__ + '._nuSQuIDS_ext'] = _sys.modules.pop(__name__ + '.nuSQuIDS')

# Clean up
del _name
