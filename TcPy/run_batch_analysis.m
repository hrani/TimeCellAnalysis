% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 3, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License 
% along with this program; see the file COPYING.  If not, write to
% the Free Software Foundation, Inc., 51 Franklin Street, Fifth
% Floor, Boston, MA 02110-1301, USA.
% 
%*******************************************************************
%* File:            run_batch_analysis.m
%* Description:     Demo for running timeCell python module in Matlab
%* Author:          K. Ananthamurthy, Upinder S. Bhalla
%* E-mail:          ananthamurthy@ncbs.res.in, bhalla@ncbs.res.in
%********************************************************************/

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

clear
close all

pe = pyenv;
if pe.Status == "NotLoaded"
    pyenv('Version', '/usr/bin/python3')
    pyenv("ExecutionMode","InProcess"); 
end

% Define a bunch of parameters for the analysis.
csOnset = 75;
usOnset = 190;
circPad = 20;
circShuffleFrames = 2 * circPad + usOnset - csOnset;
binFrames = 3;
numShuffle = 1000;

transientThresh = 2.0;
tiPercentile = 99.0;
fracTrialsFiredThresh = 0.25;
frameDt = 1.0 / 14.5;

r2bThresh = 3.0;
r2bPercentile = 99.5;

%saveFolder = '/media/ananth/Storage/Work/Analysis/RHO/M26/20180514';
saveFolder = '/home/bhalla/ananthamurthy/Work/Analysis/Imaging/M26/20180514';
bigfname = 'synthDATA_20220405_gRun3_batch_537.mat';
smallfname = 'synthDATA_20220405_gRun1_batch_36.mat';
sdcpStart = 1;
sdcpEnd = 1;        % Just doing the first dataset here.

load(strcat( saveFolder, '/', smallfname), 'sdo_batch');
for runi = sdcpStart:1:sdcpEnd
    fprintf('Dataset: %i\n', runi)

    DATA = sdo_batch(runi).syntheticDATA;
    nCells = size(DATA, 1);

    %Method A - Mehrab's Reliability Analysis (Bhalla Lab)
    [mAOutput_py] = py.pyBindMat.runR2Banalysis( DATA, csOnset, usOnset, circPad, circShuffleFrames, binFrames, numShuffle, r2bThresh, r2bPercentile );
    disp( "Finished r2b" );
    % mAOutput is a list of lists. Each cell has a list of the form:
    % [shuffled_mean, raw_r2b_ratio, r2b_bootstrap, is_mean_sig, is_bootstrap_sig, fracTrialsFired, frame_num_of_peak ]

    %Method B - William Mau's Analysis (Eichenbaum Lab)
    [mBOutput_py] = py.pyBindMat.runTIanalysis( DATA, csOnset, usOnset, circPad, circShuffleFrames, binFrames, numShuffle, transientThresh, tiPercentile, fracTrialsFiredThresh, frameDt );
    disp( "Finished TI" );
	fprintf( "      Modi et al Method   |   Mau et al method \n")
	fprintf( "cell#   meanSig bootSig   |  meanSig  bootSig\n")
    % mBOutput is a list of lists. Each cell has a list of the form:
    % [pk of mean trace, raw TI score, temporal info, is_mean_sig, is_bootstrap_sig, fracTrialsFired, frame_num_of_peak ]

    for celli = 1:nCells
        val = cell(mAOutput_py{1, celli});
        mAOutput.Q1(celli) = val{2};
        mAOutput.timeCells1(celli) = val{4};
        mAOutput.timeCells2(celli) = val{5};
        clear val
        val = cell(mBOutput_py{1, celli});
        mBOutput.Q1(celli) = val{2};
        mBOutput.timeCells1(celli) = val{4};
        mBOutput.timeCells2(celli) = val{5};
        clear val
		fprintf( "%4i      %1.1f     %1.1f     |   %1.1f       %1.1f\n", celli, mAOutput.timeCells1(celli), mAOutput.timeCells2(celli), mBOutput.timeCells1(celli), mBOutput.timeCells2(celli) )
    end
end

