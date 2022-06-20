#include <vector>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

PYBIND11_MAKE_OPAQUE(std::vector<double>);
using namespace std;
namespace py = pybind11;
#include <tcHeader.h>

vector< CellScore > r2bScore(py::array_t<double> xs, const AnalysisParams& ap, double r2bThresh, double r2bPercentile );
vector< CellScore > tiScore( py::array_t<double> xs, const AnalysisParams& ap, const TiAnalysisParams& tip );

vector< CellScore > peqScore( py::array_t<double> xs, const AnalysisParams& ap, const PeqAnalysisParams& peq );

PYBIND11_MODULE(tc, m) {
	py::bind_vector< std::vector< double >>( m, "VectorDouble");

	py::class_<AnalysisParams>( m, "AnalysisParams")
		.def( py::init())
		.def_readwrite( "csOnsetFrame", &AnalysisParams::csOnsetFrame)
		.def_readwrite( "usOnsetFrame", &AnalysisParams::usOnsetFrame)
		.def_readwrite( "circPad", &AnalysisParams::circPad)
		.def_readwrite( "circShuffleFrames", &AnalysisParams::circShuffleFrames)
		.def_readwrite( "binFrames", &AnalysisParams::binFrames)
		.def_readwrite( "numShuffle", &AnalysisParams::numShuffle)
		.def_readwrite( "epsilon", &AnalysisParams::epsilon);

	py::class_<TiAnalysisParams>( m, "TiAnalysisParams")
		.def( py::init())
		.def_readwrite( "transientThresh", &TiAnalysisParams::transientThresh)
		.def_readwrite( "tiPercentile", &TiAnalysisParams::tiPercentile)
		.def_readwrite( "fracTrialsFiredThresh", &TiAnalysisParams::fracTrialsFiredThresh )
		.def_readwrite( "frameDt", &TiAnalysisParams::frameDt );

	py::class_<PeqAnalysisParams>( m, "PeqAnalysisParams")
		.def( py::init())
		.def_readwrite( "alpha", &PeqAnalysisParams::alpha)
		.def_readwrite( "beta", &PeqAnalysisParams::beta)
		.def_readwrite( "gamma", &PeqAnalysisParams::gamma )
		.def_readwrite( "hitWindow", &PeqAnalysisParams::hitWindow );

	py::class_<CellScore>( m, "CellScore")
		.def( py::init())
		.def_readonly( "meanScore", &CellScore::meanScore)
		.def_readonly( "baseScore", &CellScore::baseScore)
		.def_readonly( "percentileScore", &CellScore::percentileScore)
		.def_readonly( "sdev", &CellScore::sdev)
		.def_readonly( "eventWidthMean", &CellScore::eventWidthMean)
		.def_readonly( "eventWidthSdev", &CellScore::eventWidthSdev)
		.def_readonly( "imprecision", &CellScore::imprecision)
		.def_readonly( "sigMean", &CellScore::sigMean)
		.def_readonly( "sigBootstrap", &CellScore::sigBootstrap)
		.def_readonly( "fracTrialsFired", &CellScore::fracTrialsFired)
		.def_readonly( "meanTrace", &CellScore::meanTrace)
		.def_readonly( "meanPkIdx", &CellScore::meanPkIdx);

	py::bind_vector< std::vector< CellScore > >( m, "VectorCellScore");


    m.doc() = "pybind11 for time cells. Take a 2d vector of (cells, trials*frames) and return a vector of time-cell scores structs"; // module docstring

    m.def("r2bScore", &r2bScore, "A function which computes r2b time-cell score for this block of data", py::arg( "data" ), py::arg( "params" ), py::arg( "r2bThresh"), py::arg( "r2bPercentile" ) );
    m.def("tiScore", &tiScore, "A function which computes temporal information time-cell score for this block of data", py::arg( "data" ), py::arg( "params" ), py::arg( "tiParams" ) );
    m.def("peqScore", &peqScore, "A function which computes parametric equation score for this block of data", py::arg( "data" ), py::arg( "params" ), py::arg( "peqParams" ) );
}

