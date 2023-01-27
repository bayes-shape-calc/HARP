import io
import os
import re
import setuptools
from setuptools.command.install import install as _install

with open("readme.md") as f:
	ld = f.read()

NAME = "HARP"

#pip version
setuptools.setup(
	name=NAME,
	version="0.0.1",
	author="Colin Kinz-Thompson",
	author_email="colin.kinzthompson@rutgers.edu",
	description="H(ierarhical) A(tomic) R(esolution) P(erception) for experimental, biomolecular structures",
	long_description=ld,
	long_description_content_type="text/markdown",
	license="GPLv3",
	package_data={"": ["LICENSE","*.so"]},
	url="",
	packages=setuptools.find_packages(where="."),
	python_requires='>=3.7',
	install_requires=["numba>=0.55","numpy>=1.22","mrcfile>=1.3","gemmi>=0.5.5"],
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
