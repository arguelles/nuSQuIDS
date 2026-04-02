import os as _os
import sys as _sys
import importlib as _importlib


def _setup_data_path():
    """Auto-configure NUSQUIDS_DATA_PATH for pip installations.

    When nuSQuIDS is installed via pip, the compiled library contains hardcoded
    paths from the build machine that don't exist on the user's system. This
    function detects available data locations and sets the environment variable
    so the C++ code can find cross-section and model data files.
    """
    if 'NUSQUIDS_DATA_PATH' in _os.environ:
        return

    # 1. Check pooch cache directory (populated by nusquids-fetch-data)
    if 'NUSQUIDS_DATA_HOME' in _os.environ:
        _cache_base = _os.environ['NUSQUIDS_DATA_HOME']
    elif _sys.platform == 'win32':
        _cache_base = _os.path.join(
            _os.environ.get('LOCALAPPDATA',
                            _os.path.join(_os.path.expanduser('~'), 'AppData', 'Local')),
            'nuSQuIDS')
    else:
        _cache_base = _os.path.join(
            _os.environ.get('XDG_DATA_HOME',
                            _os.path.join(_os.path.expanduser('~'), '.local', 'share')),
            'nuSQuIDS')

    if _os.path.isdir(_cache_base):
        # Find the latest versioned directory that has cross-section data
        try:
            for _entry in sorted(_os.listdir(_cache_base), reverse=True):
                _candidate = _os.path.join(_cache_base, _entry)
                if (_os.path.isdir(_candidate) and
                        _os.path.isfile(_os.path.join(_candidate, 'xsections',
                                                      'csms_proton.h5'))):
                    _os.environ['NUSQUIDS_DATA_PATH'] = _candidate
                    return
        except OSError:
            pass

    # 2. Check for data files bundled with the Python package
    _pkg_data = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'data')
    if _os.path.isdir(_pkg_data):
        _os.environ['NUSQUIDS_DATA_PATH'] = _pkg_data
        return


_setup_data_path()

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
