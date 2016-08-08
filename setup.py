#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_intron_retention',
      version='0.2.0',
      description='Python tools for extracting intron retention events',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/Genomon-Project/GenomonIntronRetention.git',
      package_dir = {'': 'lib'},
      packages=['genomon_intron_retention'],
      scripts=['genomon_intron_retention'],
      license='GPL-3'
     )

