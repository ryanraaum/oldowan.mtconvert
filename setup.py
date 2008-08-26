from setuptools import setup, find_packages
import sys, os

desc_lines = open('README', 'r').readlines()

setup(name='oldowan.mitomotifs',
      version='1.0.2',
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
      url='http://mitomotifs.raaum.org',
      license='MIT',
      platforms = ['Any'],
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=False,
      namespace_packages = ['oldowan'],
      zip_safe=True,
      test_suite = 'nose.collector',
      )
