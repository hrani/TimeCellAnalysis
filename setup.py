""" setup.py : Script for TimeCellAnalysis """
__author__      = "HarshaRani"
__copyright__   = "Copyright 2022 TimeCellAnalysis, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"

import os,sys,subprocess
import setuptools
from setuptools import setup, Extension,command
from setuptools.command.build_ext import build_ext

class git_clone_external(build_ext):
    def run(self):
        if not  os.path.isdir("extern/pybind11"):
        	subprocess.check_call(['git', 'clone', 'https://github.com/pybind/pybind11/','extern/pybind11/'])
        build_ext.run(self)
        
if sys.version_info >= (3,5):
	{
		setuptools.setup(
		name="TCAnalysis",
		description="The code base is to use of timeCell python module on Matlab files",
		version="1.0",
	    	author= "Upinder S. Bhalla",
	    	author_email="bhalla@ncbs.res.in",
	    	maintainer= "HarshaRani",
		maintainer_email= "hrani@ncbs.res.in",
	    	long_description = open('README.md', encoding='utf-8').read(),
	    	long_description_content_type='text/markdown',
	    	cmdclass = {'build_ext': git_clone_external},
	    	packages=['TimeCellAnalysis','TimeCellAnalysis.TcPy','TimeCellAnalysis.rho-matlab'],
	    	package_dir={'TimeCellAnalysis': "."},
	    	ext_modules=[Extension('tc', ['TcPy/tcDefaults.cpp','TcPy/tiScore.cpp','TcPy/r2bScore.cpp','TcPy/peqScore.cpp','TcPy/tcBind.cpp'], include_dirs=['extern/pybind11/include'])],
	    	install_requires = ['numpy','matplotlib','h5py'],
		url ="http://github.com/bhallalab/timecellanalysis",
                package_data = {"TimeCellAnalysis" : ['TcPy/*','rho-matlab/*',"README.md"]},
			include_package_data=True,
		license="GPLv3",
		entry_points = {
				'console_scripts' : [ 
					'TimeCellAnalysis.ti_demo = TimeCellAnalysis.__main__:run',
					'TimeCellAnalysis.ground_truth_check = TimeCellAnalysis.__main__:run_ground_truth',
					'TimeCellAnalysis.r2b_demo = TimeCellAnalysis.__main__:run_r2bdemo',
					'TimeCellAnalysis.run_batchanalysis = TimeCellAnalysis.__main__:run_batchanalysis'
						   			]
							},
					)
	}
else:
	print("Need python >= 3.5")
