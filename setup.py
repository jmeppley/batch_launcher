try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

DESCRIPTION = "Python script to fragment tasks for parallel execution"
LONG_DESCRIPTION = open('README.md').read()
NAME = "batch_launcher"
AUTHOR = "John Eppley"
AUTHOR_EMAIL = "jmeppley@gmail.com"
MAINTAINER = "John Eppley"
MAINTAINER_EMAIL = "jmeppley@gmail.com"
URL = 'http://github.com/jmeppley/batch_launcher'
DOWNLOAD_URL = 'http://github.com/jmeppley/batch_launcher'
LICENSE = 'Apache'
VERSION = '1.0.1'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      scripts=['batch_launcher.py',],
      packages=[],
      python_requires="<3.11",
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GPL License',
          'Natural Language :: English',
          'Programming Language :: Python :: 2.7'],
      )
