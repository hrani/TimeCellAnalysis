"""__main__.py: 
Entry point for this package.
"""
    
""" setup.py : Script for TimeCell Analaysis """
author__      = "HarshaRani"
__copyright__   = "Copyright 2018 FindSim, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"




def run():
    from TCAnalaysis.TcPy import timeCellMatlabExample
    timeCellMatlabExample.main()
 
if __name__ == '__main__':
    run()
    

