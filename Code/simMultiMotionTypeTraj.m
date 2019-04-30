function [traj,trackMatrix,trajTransDiffClass,errFlag] = ...
    simMultiMotionTypeTraj(numParticles,volumeEdges,totalTime,timeStep,...
    diffCoefRange,confRadRange,driftVelRange,durationRange,strictSwitch,immPrecision)
%SIMMULTIMOTIONTYPETRAJ generates trajectories that exhibit a combination of confined diffusion, free diffusion and diffusion with drift
%
%SYNOPSIS [traj,trackMatrix,trajTransDiffClass,errFlag] = ...
%    simMultiMotionTypeTraj(numParticles,volumeEdges,totalTime,timeStep,...
%    diffCoefRange,confRadRange,driftVelRange,durationRange,strictSwitch,immPrecision)
%
%INPUT  numParticles : Number of particles in simulation.
%       volumeEdges  : Edges of volume in which particles reside. A row
%                      vector with number of entries = system
%                      dimensionality.
%       totalTime    : Total time of simulation [time units].
%       timeStep     : Simulation time step [time units].
%       diffCoefRange: Row vector with 2 entries indicating diffusion
%                      coefficient range [(space units)^2/(unit time)].
%       confRadRange : Row vector with 2 entries indicating confinement
%                      radius range [space units].
%       driftVelRange: Row vector with 2 entries indicating drift speed
%                      range [space units/unit time]. Direction chosen
%                      randomly by algorithm.
%       durationRange: 4-by-2 array indicating range of time spent in each
%                      motion category. 1st row: confined diffusion; 2nd
%                      row: free diffusion; 3rd row: drift; 4th row immobile.
%                      1st column: shortest duration per category; 2nd column: longest
%                      duration per category. To exclude a category, put
%                      zeros in its row. 
%       strictSwitch : 1 to force switching to a different motion type
%                      between segments, 0 to allow switching to the same
%                      motion type but possibly with different parameters.
%
%       immPrecision: Localization precision for simulating immobile tracks
%                     (default : 0.8)
%
%OUTPUT traj       : (Number of time points) - by - (dimension) - by -
%                    (number of particles) array of particle trajectories.
%       trackMatrix: Trajectories in matrix format, as described for the
%                    input in basicTransientDiffusionAnalysisv1
%       trajTransDiffClass: Structure indicating trajectory segment diffusion
%                    classification. Format as described for the output in
%                    basicTransientDiffusionAnalysisv1.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2009
% Edit: Tony Vega 2015 add immobile option
%
% Copyright (C) 2018, Jaqaman & Danuser Labs - UTSouthwestern 
%
% This file is part of DC-MSS.
% 
% DC-MSS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% DC-MSS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with DC-MSS.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%% Output

errFlag = 0;
traj = [];
trackMatrix = [];
trajTransDiffClass = [];

%% Input

%check if correct number of arguments was used when function was called
if nargin < nargin('simMultiMotionTypeTraj')
    disp('--simMultiMotionTypeTraj: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

if nargin < 10 || isempty(immPrecision)
    immPrecision = 0.8;
end
%get system dimensionality
dimension = length(volumeEdges);

%determine total number of iterations to perform
numIterations = ceil(totalTime/timeStep) + 1;

%get possible motion types
rangeSum = sum(durationRange,2);
possibleMotionTypes = find(rangeSum~=0);
numMotionTypes = length(possibleMotionTypes);

%calculate maximum possible number of motion segments a particle can have
rangeNonZero = durationRange(durationRange~=0);
maxPossibelSegments = ceil( totalTime / min(rangeNonZero) );

%complain if there is only one possible motion type, strictSwitch = 1, and
%motion duration does not cover full trajectory
if numMotionTypes == 1 && strictSwitch == 1 && maxPossibelSegments > 1
    disp('--simMultiMotionTypeTraj: Inconsistency - please make motion duration longer than totalTime, or make strictSwitch = 0!');
    errFlag = 1;
    return
end

%% Trajectory generation

%reserve memory for output
traj = zeros(numIterations,dimension,numParticles);
trackMatrix = zeros(numParticles,8*numIterations);
segmentClass = struct('asymmetry',[],'momentScalingSpectrum',[],...
    'asymmetryAfterMSS',[]);
trajTransDiffClass = repmat(struct('segmentClass',segmentClass),numParticles,1);

%go over all particles
for iParticle = 1 : numParticles

    %get required random numbers
    randNumArray = rand(maxPossibelSegments,7);
    
    %determine the sequence of motion type segments for this particle
    %1 = confined diffusion
    %2 = free diffusion
    %3 = diffusion with drift
    %4 = immobile
    motionType = possibleMotionTypes(ceil(randNumArray(:,1)*numMotionTypes));

    if strictSwitch
        for iSegment = 2 : maxPossibelSegments
            possibleTypes = setdiff(possibleMotionTypes,motionType(iSegment-1));
            motionType(iSegment) = possibleTypes(ceil(randNumArray(iSegment,1)*(numMotionTypes-1)));
        end
    end    

    %choose a duration for each segment
    minDuration = durationRange(motionType,1);
    maxMinusMin = durationRange(motionType,2) - minDuration;
    segmentDuration = randNumArray(:,2) .* maxMinusMin + minDuration;
    
    %initialize segment classification array
    segmentClass = NaN(maxPossibelSegments,22);
    segmentClass(:,3) = motionType;

    %choose a diffusion coefficient, confinement radius and drift velocity
    %for each segment
    diffCoefSegment = randNumArray(:,3) * diff(diffCoefRange) + ...
        diffCoefRange(1);
    confRadSegment  = randNumArray(:,4) * diff(confRadRange) + ...
        confRadRange(1);
    driftVelSegment = randNumArray(:,5) * diff(driftVelRange) + ...
        driftVelRange(1);
    switch dimension
        case 1
            thetaSegment = zeros(maxPossibelSegments,0);
            phiSegment = zeros(maxPossibelSegments,0);
        case 2
            thetaSegment = randNumArray(:,6) * 360;
            phiSegment = zeros(maxPossibelSegments,0);
        case 3
            thetaSegment = randNumArray(:,6) * 360;
            phiSegment = randNumArray(:,7) * 180;
    end
    driftVelSegment = [driftVelSegment thetaSegment phiSegment]; %#ok<AGROW>

    %choose a random initial position for particle
    trajTmp = rand(1,dimension) .* volumeEdges;
    
    %loop through segments while trajectory is still shorter than totalTime
    iSegment = 0;
    while size(trajTmp,1) < numIterations
        
        %update segment number
        iSegment = iSegment + 1;
        
        %get final position of last segment, to shift trajectory
        trajShift = trajTmp(end,:);

        %generate segment
        switch motionType(iSegment)

            case 1 %confined
                [trajSegment,errFlag] = brownianMotion(dimension,...
                    diffCoefSegment(iSegment),segmentDuration(iSegment),...
                    timeStep,1,confRadSegment(iSegment));
                segmentClass(iSegment,20) = confRadSegment(iSegment);
                segmentClass(iSegment,21:22) = trajShift;
            case 2 %free
                [trajSegment,errFlag] = brownianMotion(dimension,...
                    diffCoefSegment(iSegment),segmentDuration(iSegment),...
                    timeStep);
            case 3 %drift
                [trajSegment,errFlag] = brownianMotion(dimension,...
                    diffCoefSegment(iSegment),segmentDuration(iSegment),...
                    timeStep,[],[],driftVelSegment(iSegment,:));
                
            case 4 %immobile
                %Add a bit of noise, beneath limit of confinement radius
                trajSegment = zeros(round(segmentDuration(iSegment)./timeStep)+1,2);
                trajSegment(2:end,:) = randn(round(segmentDuration(iSegment)./timeStep),2).*immPrecision;%+randi(2,1,1);%+3*ones(round(segmentDuration(iSegment)),2);
        end

        %shift segment based on last particle positions
        trajSegment = trajSegment + repmat(trajShift,size(trajSegment,1),1);
        if motionType(iSegment) == 4
            segmentClass(iSegment,3) = 0;
        end
        %store start of segment
        segmentClass(iSegment,1) = size(trajTmp,1) + 1;
        
        %append segment to end of trajectory
        trajTmp = [trajTmp; trajSegment(2:end,:)]; %#ok<AGROW>
        
        %store end of segment
        segmentClass(iSegment,2) = size(trajTmp,1);

    end %(while size(trajTmp,1) < numIterations)

    %keep only filled rows in segmentClass and modify first start and last
    %end
    segmentClass = segmentClass(~isnan(segmentClass(:,1)),:);
    segmentClass(1,1) = 1;
    segmentClass(end,2) = numIterations;
    
    %store trajectory information
    traj(:,:,iParticle) = trajTmp(1:numIterations,:);
    trackMatrix(iParticle,1:8:end) = trajTmp(1:numIterations,1);
    if dimension >= 2
        trackMatrix(iParticle,2:8:end) = trajTmp(1:numIterations,2);
    end
    if dimension >= 3
        trackMatrix(iParticle,3:8:end) = trajTmp(1:numIterations,3);
    end
    trajTransDiffClass(iParticle).segmentClass.asymmetry = [1 numIterations NaN];
    trajTransDiffClass(iParticle).segmentClass.momentScalingSpectrum = segmentClass;
    segmentClass(:,3) = NaN;
    trajTransDiffClass(iParticle).segmentClass.asymmetryAfterMSS = segmentClass(:,1:3);
    
end %(for iParticle = 1 : numParticles)
