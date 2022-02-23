import os
from setuptools import setup
import shutil

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    version = get_version(root='..', relative_to=__file__)
except Exception:
    version = '{version}'
""".lstrip()

with open('docs/source/conf.orig.py', 'r') as src:
    with open('docs/source/conf.py', 'w') as dst:
        dst.write("# THIS FILE WAS COPIED FROM conf.orig.py. DO NOT MODIFY\n")
        shutil.copyfileobj(src, dst)

setup(use_scm_version={'write_to': os.path.join('dbsp_drp', 'version.py'),
                       'write_to_template': VERSION_TEMPLATE})
