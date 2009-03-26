from setuptools import setup, find_packages
import sys, os

desc_lines = open('README', 'r').readlines()

PACKAGE = 'mtconvert'

VERSION = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),'oldowan', PACKAGE, 'VERSION')).read().strip()

setup(name='oldowan.%s' % PACKAGE,
      version=VERSION,
      description=desc_lines[0],
      long_description=''.join(desc_lines[2:]),
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Topic :: Scientific/Engineering :: Bio-Informatics"
      ],
      keywords='bioinformatics',
      author='Ryan Raaum',
      author_email='code@raaum.org',
      url='http://www.raaum.org',
      license='MIT',
      platforms = ['Any'],
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=False,
      install_requires=[
          "oldowan.polymorphism >= 1.0.0",
          "oldowan.mtdna >= 1.0.0"
      ],
      namespace_packages = ['oldowan'],
      data_files=[("oldowan/%s" % PACKAGE, ["oldowan/%s/VERSION" % PACKAGE])],
      zip_safe=False,
      test_suite = 'nose.collector',
      )
