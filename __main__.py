"""__main__.py:
Entry point for this package.
"""

""" setup.py : Script for TimeCell Analaysis """
author__      = "HarshaRani"
__copyright__   = "Copyright 2022 TimeCellAnalysis, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"

def run():
    from TimeCellAnalysis.TcPy import ti_demo
    ti_demo.main()

def run_ground_truth():
    from TimeCellAnalysis.TcPy import ground_truth_check
    ground_truth_check.main()

def run_r2bdemo():
    from TimeCellAnalysis.TcPy import r2b_demo
    r2b_demo.main()

def run_batchanalysis():
    from TimeCellAnalysis.TcPy import run_batch_analysis
    run_batch_analysis.main()

if __name__ == '__main__':
    run()
