import io
import os
import re
import setuptools
from setuptools.command.install import install as _install

with open("readme.md") as f:
	ld = f.read()

NAME = "harp"

#pip version
setuptools.setup(
	name=NAME,
	version="0.1.0",
	author="Colin Kinz-Thompson",
	author_email="colin.kinzthompson@rutgers.edu",
	description="H(ierarhical) A(tomic) R(esolution) P(erception) for experimental, biomolecular structures",
	long_description=ld,
	long_description_content_type="text/markdown",
	license="GPLv3",
	package_data={"": ["LICENSE","*.so"]},
	url="https://github.com/bayes-shape-calc/HARP",
	project_urls={
		'Documentation': 'https://github.com/bayes-shape-calc/HARP/docs',
	},
	packages=setuptools.find_packages(where="."),
	# packages=setuptools.find_packages(include=['harp','harp.*']),
	python_requires='==3.10.12',
	install_requires=[
		"numpy==1.23.5",
		"numba==0.56.4",
		"mrcfile==1.4.3",
		# "gemmi>=0.5.5"
	],
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
		"Topic :: Scientific/Engineering :: Chemistry",
		"Topic :: Scientific/Engineering :: Biology",
	],
	entry_points={
			'console_scripts': [
				'openrcsb=harp.bin.openrcsb:main',
				'harpcalc=harp.bin.harpcalc:main'
			],
	},
	include_package_data=True,
	zip_safe=True,
)
