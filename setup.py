#!/usr/bin/env python

from distutils.core import setup

setup(name='intron_retention_utils',
      version='0.3.0beta',
      description='Python tools for extracting intron retention events',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/intron_retention_utils',
      package_dir = {'': 'lib'},
      packages=['intron_retention_utils'],
      scripts=['intron_retention_utils'],
      license='GPL-3'
     )

