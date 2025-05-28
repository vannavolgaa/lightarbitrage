import sys
import os

# Add the directory containing toolsetup.py to the Python module search path
sys.path.insert(0, os.path.dirname(__file__))

from toolsetup import build_list, setup

build_type = os.environ.get('arbitrage_build', 'core')
extras = {b.get_build_name(): b.get_setup_kwargs() for b in build_list}
setup(**extras[build_type])

