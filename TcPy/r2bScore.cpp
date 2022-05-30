#include <vector>
#include <random>
#include <numeric>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Python.h>
using namespace std;
namespace py = pybind11;
#include <tcHeader.h>
PYBIND11_MAKE_OPAQUE(std::vector<double>);


//////////////////////////////////////////////////////////////////////////
// stuff for r2b calculations.
//////////////////////////////////////////////////////////////////////////

/**
 * trialAve: Computes mean over odd trials for each frame, each cell
 * Passes back vector ret[cell][frame]
 * Returns vector of peak frame number for each cell, for even trials.
 */
vector< unsigned int> trialAve( const double* data, 
		vector< vector< double > >& ret,
		unsigned int numTrials, unsigned int numCells, 
		const AnalysisParams& ap )
{
		// unsigned int offset, std::uniform_int_distribution<std::mt19937::result_type>* shuffler, std::mt19937* rng ) {
	// offset is for odd/even trials
	unsigned int nf = ap.circShuffleFrames;
	ret.resize( numCells );
	for ( auto rr = ret.begin(); rr != ret.end(); ++rr )
		rr->resize( nf, 0.0 );
	vector< unsigned int > pk;
	
	for ( unsigned int cell = 0; cell < numCells; cell++ ) {
		vector< double >& rr = ret[cell];
		vector< double > pksum( nf, 0.0 );
		for( unsigned int ff = 0; ff < nf; ++ff ) {
			const double* d = data + cell + (ff + ap.csOnsetFrame-ap.circPad)*numCells*numTrials;
			for( unsigned int tt = 0; tt < numTrials; tt+=2 ) {
				pksum[ff] += d[ tt * numCells ];
			}
			for( unsigned int tt = 1; tt < numTrials; tt+=2 ) {
				rr[ff] += d[ tt * numCells ];
			}
		}
		pk.push_back( max_element( pksum.begin(), pksum.end() ) - pksum.begin() );
		for( unsigned int ff = 0; ff < nf; ++ff ) {
			rr[ff] /= numTrials/2;
		}
	}
	return pk;
}

vector< unsigned int> shuffleTrialAve( const double* data, 
		vector< vector< double > >& ret,
		unsigned int numTrials, unsigned int numCells,
		std::uniform_int_distribution<std::mt19937::result_type>& shuffler,
	   	std::mt19937& rng, const AnalysisParams& ap ) 
{
	unsigned int nf = ap.circShuffleFrames;
	ret.resize( numCells );
	for ( auto rr = ret.begin(); rr != ret.end(); ++rr )
		rr->resize( nf, 0.0 );
	vector< unsigned int > pk;
	unsigned int ncxnt = numCells*numTrials;

	vector< int > shuff;
	for( unsigned int tt = 0; tt < numTrials; ++tt ) {
		shuff.push_back( shuffler(rng) );
	}
	
	for ( unsigned int cell = 0; cell < numCells; cell++ ) {
		vector< double >& rr = ret[cell];
		vector< double > pksum( nf, 0.0 );
		for( unsigned int tt = 0; tt < numTrials; tt+=2 ) {
			for( unsigned int ff = 0; ff < nf; ++ff ) {
				unsigned int circff = ap.csOnsetFrame - ap.circPad + (ff+shuff[tt]) % nf;
				unsigned int offset = cell + circff * ncxnt;
				pksum[ff] += data[ offset + tt * numCells ];
				rr[ff] += data[ offset + (tt+1) * numCells ];
			}
		}
		pk.push_back( max_element( pksum.begin(), pksum.end() ) - pksum.begin() );
		for( unsigned int ff = 0; ff < nf; ++ff ) {
			rr[ff] /= numTrials/2;
		}
	}
	return pk;
}

double r2b( const vector< double >& ave, unsigned int pkfr) {
	unsigned int f0 = pkfr > 0 ? pkfr-1 : 0;
	unsigned int f2 = pkfr < ave.size()-1 ? pkfr+1 : pkfr;
	double ridge = ave[f0] + ave[pkfr] + ave[f2];
	double background = std::accumulate( ave.begin(), ave.end(), 0.0 ) - ridge;
	if ( ridge < 0.0 || background <= 0.0 ) {
		return 0.0;
	}
	return ridge/background;
}

vector< CellScore > r2bScore( py::array_t<double> xs, const AnalysisParams& ap, double r2bThresh, double r2bPercentile)
{
	py::buffer_info info = xs.request();
	auto data = static_cast< double* >( info.ptr);
	unsigned int numCells = 1;
	unsigned int numTrials = 1;
	assert( info.ndim == 3 );
	numCells = info.shape[2];
	numTrials = info.shape[1];
	// numFrames = info.shape[0];
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> shuffler(0, ap.circShuffleFrames);

	vector< CellScore > ret( numCells );
	vector< vector< double > > tave;
	auto pkfr = trialAve( data, tave, numTrials, numCells, ap );
	vector< double > r2bOriginal;;

	for(unsigned int cc = 0; cc < numCells; cc++){
		r2bOriginal.push_back( r2b( tave[cc], pkfr[cc] ) );
		ret[cc].meanTrace = tave[cc];
		ret[cc].meanPkIdx = pkfr[cc];
    }
	
	vector< double > sumShuff( numCells, 0.0 );
	vector< double > r2bBootstrap( numCells, 0.0 );
	for( unsigned int ii = 0; ii < ap.numShuffle; ++ii ) {
		auto pkfr = shuffleTrialAve( data, tave, numTrials, numCells, shuffler, rng, ap );
		for(unsigned int cc = 0; cc < numCells; cc++){
			double rr = r2b( tave[cc], pkfr[cc] );
			if ( !isnan( rr ) ) {
				sumShuff[cc] += rr;
			}
			r2bBootstrap[cc] += ( r2bOriginal[cc] > rr );
		}
	}
	for(unsigned int cc = 0; cc < numCells; cc++) {
		if ( abs( sumShuff[cc] ) > ap.epsilon ) {
			ret[cc].meanScore = sumShuff[cc] / double(ap.numShuffle);
			ret[cc].baseScore = r2bOriginal[cc];
		} else {
			ret[cc].baseScore = 0.0;
		}
		ret[cc].percentileScore = r2bBootstrap[cc] / double(ap.numShuffle);
		ret[cc].sigMean = ( ret[cc].baseScore > (r2bThresh*ret[cc].meanScore) );
		ret[cc].sigBootstrap = ( ret[cc].percentileScore > r2bPercentile/100.0 );
	}
	return ret;
}
