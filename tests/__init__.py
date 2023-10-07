import sys

MIN_PYTHON = (3, 6)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

import pkg_resources
sys.path.append(pkg_resources.resource_filename('seq', ''))
sys.path.append(pkg_resources.resource_filename('tests', ''))
