#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
using namespace std;
namespace py = pybind11;

#include "tcHeader.h"

PYBIND11_MAKE_OPAQUE(std::vector<double>);

AnalysisParams::AnalysisParams()
	: 
		csOnsetFrame( 75 ),
		usOnsetFrame( 190 ),
		circPad( 20 ),
		binFrames( 3 ),
		numShuffle( 1000 ),
		epsilon( 1.0e-6 )
{
	circShuffleFrames = circPad*2 + usOnsetFrame - csOnsetFrame;;
}

TiAnalysisParams::TiAnalysisParams()
	: 
		transientThresh( 2.0 ),
		tiPercentile( 99.0 ),

		// frac of trials where at least 1 peak>transientThresh*sdev + mean
		fracTrialsFiredThresh( 0.25 ),
		frameDt( 0.08 )
{;}

PeqAnalysisParams::PeqAnalysisParams()
	: 
		alpha( 10.0 ),
		beta( 1.0 ),
		gamma( 10.0 ),
		transientThresh( 2.0 ),
		hitWindow( 5.0 ) // +/- 5 frames.
{;}
