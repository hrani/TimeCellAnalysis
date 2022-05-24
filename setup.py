""" setup.py : Script for FindSim """
__author__      = "HarshaRani"
__copyright__   = "Copyright 2021 HillTau, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"

import os
import sys
sys.path.append(os.path.dirname(__file__))
import setuptools
#from distutils.core import setup, Extension
from setuptools import setup, Extension
setuptools.setup(
        name="TCAnalaysis",
        description="The code base is to use of timeCell python module on Matlab files",
        version="1.0",
        author= "Upinder S. Bhalla",
        author_email="bhalla@ncbs.res.in",
        maintainer= "HarshaRani",
        email= "hrani@ncbs.res.in",
        long_description = open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        packages=["TimeCellAnalaysis","TimeCellAnalaysis.TcPy"],
        package_dir={'TimeCellAnalaysis': "."},
        ext_modules=[Extension('tc', ['TcPy/tcDefaults.cpp','TcPy/tiScore.cpp'], include_dirs=['extern/pybind11/include'])],
        install_requires = ['numpy','matplotlib','h5py'],
        
        #dependency_links=[
        # Make sure to include the `#egg` portion so the `install_requires` recognizes the package
        #'git+ssh://git@github.com/pybind/pybind11.git#egg=ExampleRepo-0.1'],
	url ="http://github.com/bhallalab/timecellanalaysis",
	package_data = {"TCAnalaysis" : ['TcPy/timeCellMatlabExample.py']},
	license="GPLv3",
	entry_points = {
	
		'console_scripts' : [ 'TimeCellAnalaysis = TimeCellAnalaysis.__main__:run'
				   ]
			},
		)



""" setup.py : Script for FindSim """
__author__      = "HarshaRani"
__copyright__   = "Copyright 2021 HillTau, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"

import os
import sys
sys.path.append(os.path.dirname(__file__))
import setuptools
#from distutils.core import setup, Extension
from setuptools import setup, Extension
setuptools.setup(
        name="TCAnalaysis",
        description=" ",
        version="1.0",
        author= "Upinder S. Bhalla",
        author_email="bhalla@ncbs.res.in",
        maintainer= "HarshaRani",
        maintainer_email= "hrani@ncbs.res.in",
        long_description = open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        packages=["TimeCellAnalaysis","TimeCellAnalaysis.TcPy"],
        package_dir={'TCAnalaysis': "."},
        ext_modules=[Extension('tc', ['TcPy/tcDefaults.cpp'], include_dirs=['extern/pybind11/include'])],
        install_requires = ['numpy','matplotlib','h5py'],
        
        #dependency_links=[
        # Make sure to include the `#egg` portion so the `install_requires` recognizes the package
        #'git+ssh://git@github.com/pybind/pybind11.git#egg=ExampleRepo-0.1'],
	url ="http://github.com/bhallalab/TimeCellAnalaysis",
	package_data = {"TCAnalaysis" : ['TcPy/timeCellMatlabExample.py']},
	license="GPLv3",
	entry_points = {
	
		'console_scripts' : [ 'TCAnalaysis = TCAnalaysis.__main__:run'
				   ]
			},
		)

