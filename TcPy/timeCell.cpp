#include <vector>
#include <random>
#include <numeric>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Python.h>

PYBIND11_MAKE_OPAQUE(std::vector<double>);
using namespace std;
namespace py = pybind11;
// py::module nn = py::module::import("iteration");
#define BASELINE_START_FRAME 20
#define CS_ONSET_FRAME 75
#define US_ONSET_FRAME 190
/// CIRC_PAD is padding on either side of CS and US.
#define CIRC_PAD 20

/// CIRC_SHUFFLE_FRAMES is the size of the circular dataset to be shuffled
#define CIRC_SHUFFLE_FRAMES CIRC_PAD*2 + US_ONSET_FRAME - CS_ONSET_FRAME
#define NUM_R2B_SHUFFLE 1000
#define PERCENTILE 10
#define EPSILON 1.0e-6
#define NUM_TI_SHUFF 1000

// From Mau et al 2018
#define TI_TRANSIENT_THRESH 2.0
#define TI_PERCENTILE 99.0
#define TI_BIN_FRAMES 3
#define TI_FRAC_TRIALS_FIRED 0.25


struct CellScore {
	double baseTI;
	double percentileTI;
	double fracTrialsFired;
	vector< double > tuningCurve;
	vector< vector< unsigned int > > transients; // Temporary, for debug.
	unsigned int tuningPk;
	bool tuningIsSig;
	bool temporalInfoIsSig;
};

/// Returns the mean dfbf[ frame# ] for the specified shuffle[ trial# ]
/// using data[ trial# ][ frame# ]. Operates on data for a single cell.
vector< double > aveOfTrials( 
				const vector< vector< double > >& data, 
				const vector< unsigned int >& shuff, 
				unsigned int numFrames )
{
	unsigned int numTrials = data.size();
	vector< double > aveBin( numFrames, 0.0 );
	assert( numTrials == data.size() && numTrials > 0  );
	for ( unsigned int ii = 0; ii < numTrials; ++ii ) {
		vector< double >::iterator aptr = aveBin.begin();
		auto dd = data[ii].begin();
		for ( unsigned int bb = shuff[ii]; bb < shuff[ii] + numFrames; bb++, aptr++) {
			*aptr += *( dd + bb%numFrames );
		}
	}
	for( vector< double >::iterator aptr = aveBin.begin(); aptr != aveBin.end(); aptr++ ) {
		*aptr /= numFrames;
	}
	return aveBin;
}

/*
vector< double > aveOfTrials( const vector< unsigned int >& shuff, const vector< vector< double > >& data, unsigned int numFrames )
{
	vector< double > aveBin( numFrames, 0.0 );
	assert( shuff.size() == data.size() && shuff.size() > 0  );
	assert( data[0].size() >= numFrames * 2 ); // Space for circular shuffle
	for ( unsigned int ii = 0; ii < shuff.size(); ++ii ) {
		vector< double >::iterator aptr = aveBin.begin();
		auto bb = data[ii].begin() + shuff[ii];
		auto ee = data[ii].begin() + shuff[ii] + numFrames;
		for ( vector< double >::const_iterator dd = bb; dd != ee; ++dd ) {
			*aptr++ += *dd;
		}
	}
	for( auto aptr = aveBin.begin(); aptr != aveBin.end(); aptr++ ) {
		*aptr /= numFrames;
	}
	return aveBin;
}
*/

/// returns binned data for a single trial.
vector< double > binner( const vector< double >& orig, unsigned int binFrames )
{
	if (binFrames <= 1 || binFrames >= orig.size()/2 ) return orig;
	unsigned int numRet = orig.size() / binFrames;
	vector< double > ret( numRet, 0.0 );
	for ( unsigned int idx = 0; idx < numRet; idx++ ) {
		for ( unsigned int bf = 0; bf < binFrames; bf++ )
			ret[idx] += orig[ idx * binFrames + bf ];
	}
	return ret;
}

/// returns binned data for a single cell, all trials.
vector< vector< double > > binnedCellData( const vector< vector< double > >& data, unsigned int numFrames, unsigned int binFrames )
{
	unsigned int numTrials = data.size();
	vector< vector< double > > ret( numTrials );
	for( unsigned int trialIdx = 0; trialIdx < numTrials; trialIdx++ ) {
		ret[ trialIdx ] = binner( data[ trialIdx ], binFrames );
	}
	return ret;
}

vector< unsigned int> findPeaks( const vector< double >& data ) {
	double pk = data[0];
	unsigned int pkIdx = 0;
	unsigned int otherPkIdx = 0;
	unsigned int numBins = data.size();
	for( unsigned int ii = 1; ii < numBins; ++ii ) {
		if ( pk < data[ii] ) {
			pk = data[ii];
			pkIdx = ii;
		}
	}
	if (pkIdx == 0) {
		otherPkIdx = 1;
	} else if (pkIdx == numBins - 1) {
		otherPkIdx = numBins - 2;
	} else if ( data[pkIdx-1] > data[pkIdx+1] ) {
		otherPkIdx = pkIdx-1;
	} else {
		otherPkIdx = pkIdx+1;
	}
	vector< unsigned int > ret( {pkIdx, otherPkIdx} );
	return ret;
}

/// Return the average tuning curve of the cell, and pass back an arg with
/// the bin# of the peak iff significant. Else return ~0U.
void tuningCurve( const vector< vector< double >>& data, const vector< vector< unsigned int > >& shuff, unsigned int numFrames, unsigned int binFrames, CellScore& cs )
{
	// Get the binned data
	unsigned int numTrials = data.size();
	vector< vector< double >> binnedCellData_ = binnedCellData( data, numFrames, binFrames );
	unsigned int numBins = binnedCellData_[0].size();
	vector< unsigned int > zeros( numTrials, 0 );
	// Get the reference averaged binned trial data.
	cs.tuningCurve = aveOfTrials( binnedCellData_, zeros, numBins );

	// Find the 2 consecutive biggest. Assume peak is one of them.
	vector< unsigned int > pks = findPeaks( cs.tuningCurve );
	cs.tuningPk = pks[0];

	// Now I could be efficient and just average over the two bins
	// here, as per shuffling. Or I could use the simple aveOfTrials.
	unsigned int numOK = 0;
	for ( unsigned int ii = 0; ii < NUM_TI_SHUFF; ii++ ) {
		vector< double > temp = aveOfTrials( binnedCellData_, shuff[ii], numBins );
		numOK += ( cs.tuningCurve[pks[0]] > temp[pks[0]] ) && 
				( cs.tuningCurve[pks[1]] > temp[pks[1]] );
	}
	cs.tuningIsSig = (  (numOK * 100) > (99 * NUM_TI_SHUFF) );
}

/// args: data[ frame# ][ trial# ][ cell# ]
/// Pass-back arg: ret[ cell# ][ trial# ][ frame# ]
/// Fills in only the data in the window around the stimuli, from 
/// CS_ONSET_FRAME - CIRC_PAD to US_ONSET_FRAME + CIRC_PAD.
void reorderData( const double* data, unsigned int numCells, unsigned int numTrials, unsigned int numFrames, vector< vector< vector< double >>>& ret ) 
{
	ret.clear();
	ret.resize( numCells );
	unsigned int ncxnt = numCells * numTrials;
	unsigned int nf = CIRC_SHUFFLE_FRAMES;
	for ( unsigned int cell = 0; cell < numCells; cell++ ) {
		vector< vector< double > >& rc = ret[cell];
		rc.resize( numTrials );
		for( unsigned int tt = 0; tt < numTrials; tt++ ) {
			vector< double >& rct = rc[ tt ];
			rct.resize( nf,  0.0 );
			for( unsigned int ff = 0; ff < nf; ++ff ) {
				rct[ff] = data[ cell + tt * numCells + (ff + CS_ONSET_FRAME - CIRC_PAD) * ncxnt];
			}
		}
	}
}

// Returns threshold for sig peak for data[trial#][frame#]. The
// data has already been sub-selected for a given cell.
double findCellThresh(
		const vector< vector< double > >& data, unsigned int numFrames )
{
	double sum = 0.0;
	double sq = 0.0;
	for( auto trialD = data.begin(); trialD != data.end(); ++trialD ) {
		for( auto frameD = trialD->begin(); frameD != trialD->end(); ++frameD ) {
			sum += *frameD;
			sq += *frameD * *frameD;
		}
	}
	double numSamples = numFrames * data.size(); //numFrames * numTrials
	double mean = sum/numSamples;
	return mean + TI_TRANSIENT_THRESH * sqrt( sq / numSamples - mean*mean );
}

/// Returns time stamps of transients trial# ][ transient# ]. 
// given the reordered dfbf data[ trial# ][ frame# ]
// A given trial may have zero or more transients.
vector< vector< unsigned int > > findTransients(
		const vector< vector< double > >& data, double cellThresh, unsigned int numFrames )
{
	unsigned int numTrials = data.size();
	vector< vector< unsigned int > > ret ( numTrials );
	vector< vector< unsigned int > >::iterator trialT = ret.begin();
	for( auto trialD = data.begin(); trialD != data.end(); ++trialD, ++trialT ) {
		trialT->clear();
		double lastFrame = 0.0;
		bool refractory = 0;
		for( unsigned int ii = 0; ii < numFrames; ++ii ) {
			double frameD = (*trialD)[ii];
			// Don't permit another transient till signal goes < cellThresh
			if (frameD > cellThresh && frameD > lastFrame && !refractory) {
				trialT->push_back( ii );
				refractory = 1;
			} else {
				refractory = (frameD > cellThresh);
			}
			lastFrame = frameD;
		}
	}
	return ret;
}

// Returns transient rates for the current
double findAveTransientRate( const vector< vector< unsigned int > >& transients, double trialDuration, unsigned int numTrials ) 
{
	double ret = 0.0;
	for( auto trial: transients ) {
		for ( auto bin: trial ) {
			ret += bin;
		}
	}
	return ret/ ( trialDuration * numTrials );
}

/// Returns temporal information for the current cell.
double findTemporalInformation( const vector< vector< unsigned int > >& transients, double aveTransientRate, const vector< unsigned int >& shuff, unsigned int numFrames )
{
	double ret = 0.0;
	unsigned int numTrials = transients.size();

	vector< double > binAve( numFrames, 0.0 );
	for( unsigned int trialIdx = 0; trialIdx < numTrials; ++trialIdx ) {
		for ( unsigned int bin: transients[trialIdx] ) {
			unsigned int b = ( bin + shuff[trialIdx] ) % numFrames;
			binAve[b] += 1.0;
		}
	}
	for ( auto bb: binAve ) {
		if ( bb > 0.0 ) {
			ret += bb * log2( bb / aveTransientRate );
		}
	}
	return ret / ( aveTransientRate * numFrames );
}

// Returns the specified percentile value of shuffled TIs, for current cell
double percentileTI( const vector< vector< unsigned int > >& transients, double aveTransientRate, unsigned int numFrames, double percentile )
{
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> shuffler(0,numFrames);
	unsigned int numTrials = transients.size();
	assert( numTrials > 0 );

	// allocate tiRank[iter#], for sorting.
	vector< double > tiRank( NUM_TI_SHUFF ); 

	for ( unsigned int iter = 0; iter < NUM_TI_SHUFF; ++iter ) {
		vector< unsigned int > shuff( numTrials );
		for ( unsigned int tt = 0; tt < numTrials; ++tt ) {
			shuff[tt] = shuffler( rng );
		}
		tiRank[iter] = findTemporalInformation( transients, aveTransientRate, shuff, CIRC_SHUFFLE_FRAMES  );
	}

	unsigned int percIdx = (percentile * NUM_TI_SHUFF) / 100;
	sort( tiRank.begin(), tiRank.end() );
	return tiRank[percIdx];
}

// Fraction of trials in which there was >= 1 transient. frac[ cell# ]
double fracTrialsFired( const vector< vector< unsigned int > >& transients )
{
	assert( transients.size() > 0 );
	double ret = 0.0;
	for( auto tt : transients ) {
		ret += (tt.size() > 0 );
	}
	return ret / transients.size();
}

/// Computes TI score for a given neuron. zero means it fails a criterion.
// Should really return a lot of stuff: the mean bin vector, the
// baseTI, the PTI, the frac fired, the percentile bin height.
// Suggest return a 2-D vector or a small struct.
CellScore cellTIscore( const vector< vector< double > >& data, double frameDt )
{
	CellScore cs;
	unsigned int numTrials = data.size();
	assert( numTrials > 0 );
	unsigned int numFrames = data[0].size();
	vector< unsigned int > nonShuff( numTrials, 0 );
	vector< double > aveOfTrials_ = aveOfTrials( data, nonShuff, CIRC_SHUFFLE_FRAMES );

	double cellThresh = findCellThresh( data, CIRC_SHUFFLE_FRAMES );
	vector< vector< unsigned int > > transients =
			findTransients( data, cellThresh, CIRC_SHUFFLE_FRAMES  );

	double trialDuration = frameDt * CIRC_SHUFFLE_FRAMES;
	double aveTransientRate = findAveTransientRate( transients, trialDuration, numTrials );
	cs.baseTI = findTemporalInformation( transients, aveTransientRate, nonShuff, CIRC_SHUFFLE_FRAMES  );
	cs.percentileTI = percentileTI( transients, aveTransientRate, CIRC_SHUFFLE_FRAMES, TI_PERCENTILE );
	cs.fracTrialsFired = fracTrialsFired( transients );
	cs.transients = transients;

	// Fill up the RNGs
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> shuffler(0, numFrames / TI_BIN_FRAMES );
	vector< vector< unsigned int > > binShuff( NUM_TI_SHUFF );
	for ( unsigned int ii = 0; ii < NUM_TI_SHUFF; ii++ ) {
		binShuff[ii].resize( numTrials );
		for ( unsigned int jj = 0; jj < numTrials; jj++ ) {
			binShuff[ii][jj] = shuffler( rng );
		}
	}
	tuningCurve( data, binShuff, numFrames, TI_BIN_FRAMES, cs );

	cs.temporalInfoIsSig = ( cs.baseTI > cs.percentileTI ) && (cs.fracTrialsFired > TI_FRAC_TRIALS_FIRED );

	return cs;
}


// Returns the reliability index.
vector< double > tiScore( py::array_t<double> xs, double frameDt )
{
	py::buffer_info info = xs.request();
	auto data = static_cast< double* >( info.ptr);
	unsigned int numCells = 1;
	// unsigned int numDat = info.itemsize;
	unsigned int numTrials = 1;
	unsigned int numFrames = 1;
	assert( info.ndim == 3 );
	numCells = info.shape[2];
	// numDat = info.shape[0] * info.shape[1];
	numTrials = info.shape[1];
	numFrames = info.shape[0];

	vector< vector< vector< double >>> reorderedData;
	reorderData( data, numCells, numTrials, numFrames, reorderedData );

	vector< double > ret( numCells );
	for( unsigned int cellIdx = 0; cellIdx < numCells; cellIdx++ ) {
		CellScore cs = cellTIscore( reorderedData[cellIdx], frameDt );
		ret[cellIdx] = cs.tuningIsSig * 2 + cs.temporalInfoIsSig;
	}

	return ret;
}


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
		unsigned int numTrials, unsigned int numCells )
{
		// unsigned int offset, std::uniform_int_distribution<std::mt19937::result_type>* shuffler, std::mt19937* rng ) {
	// offset is for odd/even trials
	unsigned int nf = CIRC_SHUFFLE_FRAMES;
	ret.resize( numCells );
	for ( auto rr = ret.begin(); rr != ret.end(); ++rr )
		rr->resize( nf, 0.0 );
	vector< unsigned int > pk;
	
	for ( unsigned int cell = 0; cell < numCells; cell++ ) {
		vector< double >& rr = ret[cell];
		vector< double > pksum( nf, 0.0 );
		for( unsigned int ff = 0; ff < nf; ++ff ) {
			const double* d = data + cell + (ff + CS_ONSET_FRAME-CIRC_PAD)*numCells*numTrials;
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
	   	std::mt19937& rng ) 
{
	unsigned int nf = CIRC_SHUFFLE_FRAMES;
	ret.resize( numCells );
	for ( auto rr = ret.begin(); rr != ret.end(); ++rr )
		rr->resize( nf, 0.0 );
	vector< unsigned int > pk;
	unsigned int ncxnt = numCells*numTrials;

	vector< int > shuff;
	for( unsigned int tt = 0; tt < numTrials; ++tt ) {
		shuff.push_back( shuffler(rng) );
		// cout << shuff.back() << endl;
	}
	
	for ( unsigned int cell = 0; cell < numCells; cell++ ) {
		vector< double >& rr = ret[cell];
		vector< double > pksum( nf, 0.0 );
		for( unsigned int tt = 0; tt < numTrials; tt+=2 ) {
			for( unsigned int ff = 0; ff < nf; ++ff ) {
				unsigned int circff = CS_ONSET_FRAME - CIRC_PAD + (ff+shuff[tt]) % nf;
				unsigned int offset = cell + circff * ncxnt;
				pksum[ff] += data[ offset + tt * numCells ];
				rr[ff] += data[ offset + (tt+1) * numCells ];
			}
		}
		/*
		for( unsigned int ff = 0; ff < nf; ++ff ) {
			unsigned int offset = cell + (ff + CS_ONSET_FRAME-PAD_FRAMES)*ncxnt;
			for( unsigned int tt = 0; tt < numTrials; tt+=2 ) {
				pksum[ff] += data[ offset + tt * numCells + shuff[tt] * ncxnt ];
			}
			for( unsigned int tt = 1; tt < numTrials; tt+=2 ) {
				rr[ff] += data[ offset + tt * numCells + shuff[tt] * ncxnt];
			}
		}
		*/
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
		// cout << "-----ve Ridge= " << ridge << ",	bg= " << background << endl;
	}
	// vector< double >ret( { double(pkfr), ridge, background, ridge/background } );
	return ridge/background;
}

// Returns vector of len 2*numCells. First numCells is reliability ratio.
// Second numCells is bootstrap percentile, where anything above 0.99 is sig
vector< double > r2bScore( py::array_t<double> xs )
{
	py::buffer_info info = xs.request();
	auto data = static_cast< double* >( info.ptr);
	unsigned int numCells = 1;
	// unsigned int numDat = info.itemsize;
	unsigned int numTrials = 1;
	// unsigned int numFrames = 1;
	assert( info.ndim == 3 );
	numCells = info.shape[2];
	// numDat = info.shape[0] * info.shape[1];
	numTrials = info.shape[1];
	// numFrames = info.shape[0];
	std::random_device dev;
	std::mt19937 rng(dev());
	// Do I need to seed this?
	std::uniform_int_distribution<std::mt19937::result_type> shuffler(0, CIRC_SHUFFLE_FRAMES);

	vector< vector< double > > tave;
	auto pkfr = trialAve( data, tave, numTrials, numCells );
	vector< double > r2bOriginal;;

	// cout << "Cell #		r2b original\n";
	for(unsigned int cc = 0; cc < numCells; cc++){
		r2bOriginal.push_back( r2b( tave[cc], pkfr[cc] ) );
    }
	
	vector< double > sumShuff( numCells, 0.0 );
	vector< double > r2bBootstrap( numCells, 0.0 );
	for( unsigned int ii = 0; ii < NUM_R2B_SHUFFLE; ++ii ) {
		auto pkfr = shuffleTrialAve( data, tave, numTrials, numCells, shuffler, rng );
		for(unsigned int cc = 0; cc < numCells; cc++){
			double rr = r2b( tave[cc], pkfr[cc] );
			if ( !isnan( rr ) ) {
				sumShuff[cc] += rr;
			}
			r2bBootstrap[cc] += ( r2bOriginal[cc] > rr );
		}
	}
	vector< double > ret;
	for(unsigned int cc = 0; cc < numCells; cc++) {
		if ( abs( sumShuff[cc] ) > EPSILON ) {
			ret.push_back( double(NUM_R2B_SHUFFLE) * r2bOriginal[cc] / sumShuff[cc] );
		} else {
			ret.push_back( 0.0 );
		}
		r2bBootstrap[cc] /= NUM_R2B_SHUFFLE;
	}
	ret.insert( ret.end(), r2bBootstrap.begin(), r2bBootstrap.end() );
	return ret;
}

PYBIND11_MODULE(tc, m) {
	py::bind_vector< std::vector< double >>( m, "VectorDouble");
    m.doc() = "pybind11 for time cells. Take a 2d vector of (cells, trials*frames) and return a vector of time-cell scores"; // module docstring

    m.def("r2bScore", &r2bScore, "A function which computes r2b time-cell score for this block of data", py::arg( "data" ) );
    m.def("tiScore", &tiScore, "A function which computes temporal information time-cell score for this block of data", py::arg( "data" ), py::arg( "frameDt" ) );
}

