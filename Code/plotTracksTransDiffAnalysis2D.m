function plotTracksTransDiffAnalysis2D(trackedFeatureInfo,transDiffAnalysisRes,timeRange,...
    newFigure,image)
%PLOTTRACKSTRANSDIFFANALYSIS plots tracks in 2D highlighting the different diffusion segments within each track
%
%SYNOPSIS plotTracksTransDiffAnalysis2D(trackedFeatureInfo,transDiffAnalysisRes,timeRange,...
%    newFigure,image)
%
%INPUT  trackedFeatureInfo: Same input, tracks, used for function, 
%                           basicTransientDiffusionAnalysisv1
%       transDiffAnalysisRes   : Output of the function,
%                           basicTransientDiffusionAnalysisv1
%       timeRange         : 2-element row vector indicating time range to plot. 
%                           Optional. Default: whole movie.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image.
%
%OUTPUT The plot.
%       Color coding:
%        Brown: Immobile
%        Blue: Confined diffusion
%        Cyan: Free diffusion
%        Magenta: Directed diffusion
%        Black: unclassified
%
%Khuloud Jaqaman, April 2009
%Edited: Tony Vega 2016
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

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--plotTracksTransDiffAnalysis2D: Incorrect number of input arguments!');
    return
end

%get number of tracks and number of time points
if isstruct(trackedFeatureInfo) %if tracks are in structure format
    numTracks = length(trackedFeatureInfo);
    tmp = vertcat(trackedFeatureInfo.seqOfEvents);
    numTimePoints = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
end

errFlag = 0;

%check whether a time range for plotting was input
if nargin < 3 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracksTransDiffAnalysis2D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 4 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1
        disp('--plotTracksTransDiffAnalysis2D: newFigure should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether user supplied an image
if nargin < 5 || isempty(image)
    image = [];
end

%check whether to plot confinement radii
% if nargin < 6 || isempty(showConf)
    showConf = 0;
% end

%check whether to plot confinement radii
% if nargin < 7 || isempty(checkAsym)
    checkAsym = 0;
% end

%exit if there are problems in input variables
if errFlag
    disp('--plotTracksTransDiffAnalysis2D: Please fix input data!');
    return
end

%% Pre-processing

if isstruct(trackedFeatureInfo) %if tracks are input in structure format

    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end

    %if all tracks have only one segment ...
    if max(numSegments) == 1

        %indicate that there are no compound tracks with merging and splitting branches
        mergeSplit = 0;

        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step) 
        %in this case of course every compound track is simply one track
        %without branches
        trackStartRow = (1:numTracks)';

        %store tracks in a matrix
        trackedFeatureInfo = NaN*ones(numTracks,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(i,8*(startTime-1)+1:8*endTime) = inputStructure(i).tracksCoordAmpCG;
        end
        
    else %if some tracks have merging/splitting branches
        
        %indicate that in the variable mergeSplit
        mergeSplit = 1;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        trackStartRow = ones(numTracks,1);
        for iTrack = 2 : numTracks
            trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);            
        end
        
        %put all tracks together in a matrix
        trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
                numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
                inputStructure(i).tracksCoordAmpCG;
        end
        
    end    
    
else %if tracks are not input in structure format

    %indicate that there are no compound tracks with merging and splitting branches
    mergeSplit = 0;
    
    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';

end

%get the x,y-coordinates of features in all tracks
tracksX = trackedFeatureInfo(:,1:8:end)';
tracksY = trackedFeatureInfo(:,2:8:end)';

%find x-coordinate limits
minXCoord = min(floor(min(tracksX(:))),0);
maxXCoord =  ceil(max(tracksX(:)));

%find y-coordinate limits
minYCoord = min(floor(min(tracksY(:))),0);
maxYCoord =  ceil(max(tracksY(:)));

%get number of track segments to be plotted
numTrackSegments = size(tracksX,2);

% %% confinement radius information
% 
% %get track segment center, confinement radii and preferred direction of
% %motion
% trackSegmentCenter = catStruct(1,'diffAnalysisRes.confRadInfo.trackCenter');
% trackSegmentConfRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');
% trackSegmentPrefDir = catStruct(1,'diffAnalysisRes.confRadInfo.prefDir');
% 
% %determine indices of tracks with one confinement radius
% indxCircle = find( ~isnan(trackSegmentConfRad(:,1)) & isnan(trackSegmentConfRad(:,2)) );
% 
% %determine indices of tracks with 2 confinement radii
% indxRectangle = find( ~isnan(trackSegmentConfRad(:,2)) );

%% Plotting

%if the user wants to plot in a new figure window
if newFigure

    %open new figure window
    figure

    if ~isempty(image) %if user supplied an image
        imshow(image,[]); %plot the image
    else %if user did not supply an image
        imshow(ones(maxYCoord,maxXCoord),[]); %plot an empty image
    end

    %set figure axes limits
    axis([minXCoord maxXCoord minYCoord maxYCoord]);

    %show coordinates on axes
    ah = gca;
    set(ah,'visible','on');

    %label axes
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');

end

%hold on figure
hold on

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

% Extracting diffusion data
% Build 3D matrix to store indices for all diffusion types. Rows correspond to
% time frame, columns correspond to specific track and  z stack corresponds
% to diffusion type.
% Two matrices are built to store symmetric and asymmetric classification

base = NaN([size(tracksX,1) size(tracksX,2) 4]); %Create 3D matrix to store indices for all diffusion types
asymBase = NaN([size(tracksX,1) size(tracksX,2) 4]);
nanBase = NaN(size(tracksX));
diffTypesF = [];

%  if isempty(bayes) 
        %get track segment types from diffusion analysis
        trackSegmentType = vertcat(transDiffAnalysisRes.segmentClass);
        
        % For a given track, assign all time points a given
        % diffusion type. The assignment of a given diffusion type
        % corresponds to a specific z-stack. Thus the first level
        % of the z-stack will be populated with time points of all
        % tracks that are immobile, for example
        if checkAsym == 1
            for k = 1:length(trackSegmentType)
                for j = 1:size(trackSegmentType(k).momentScalingSpectrum,1)
                    if ~isnan(trackSegmentType(k).momentScalingSpectrum(j,3))
                        base(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k) = ...
                        trackSegmentType(k).momentScalingSpectrum(j,3);

                        asymBase(trackSegmentType(k).momentScalingSpectrum1D(j,1):trackSegmentType(k).momentScalingSpectrum1D(j,2),k) = ...
                        trackSegmentType(k).momentScalingSpectrum1D(j,3);
                        diffTypesF = [diffTypesF ; trackSegmentType(k).momentScalingSpectrum(j,3)];
                    else
                        %If 2D is unclassified, check 1D classification
                        if ~isnan(trackSegmentType(k).momentScalingSpectrum1D(j,3))
                            asymBase(trackSegmentType(k).momentScalingSpectrum1D(j,1):trackSegmentType(k).momentScalingSpectrum1D(j,2),k) = ...
                            trackSegmentType(k).momentScalingSpectrum1D(j,3);
                        
                            nanBase(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k) = -1;
                            diffTypesF = [diffTypesF ; -1];                        
                        else
                            %Otherwise, set postions in nanBase to -1 and store classification as -1;
                            nanBase(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k) = -1;
                            diffTypesF = [diffTypesF ; -1];
                        end
                    end

                end

            end
        else
            for k = 1:length(trackSegmentType) % For each track

                for j = 1:size(trackSegmentType(k).momentScalingSpectrum,1) % For each segment
                    %If this segment has a classification, continue with
                    %process and store the type of classification
                    if ~isnan(trackSegmentType(k).momentScalingSpectrum(j,3))
                        base(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k,trackSegmentType(k).momentScalingSpectrum(j,3)+1) = ...
                        trackSegmentType(k).momentScalingSpectrum(j,3);
                        diffTypesF = [diffTypesF ; trackSegmentType(k).momentScalingSpectrum(j,3)];
                    else
                        %Otherwise, set postions in nanBase to -1 and store classification as -1;
                        nanBase(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k) = -1;
                        diffTypesF = [diffTypesF ; -1];
                    end
                end

            end
        end
%  else
%    for m = 1:size(tracksXP,2)
%         results = diffAnalysisRes;
%         states = results.ML_states{1,m};
%         track = results.track{1,m};
%         steps = NaN(2,size(track,2));
%         steps(:,2:end) = results.steps{1,m};
% 
%         track(1,~isnan(steps(1,:))) = states;
%         track(2,~isnan(steps(1,:))) = states;
%         index = find(isnan(steps(2,:)));
%         track(1,index)=track(2,index+1);
%         track(2,index)=track(2,index+1);
%         base(~isnan(tracksXP(:,m)),m) = track(1,:);
%     end  
%      
%  end

numTimePlot = timeRange(2) - timeRange(1) + 1;


%Get all diffusion types present
diffTypes = unique(diffTypesF);

% Plot dashed line of entire track to allow gaps in track to be easily
% visible
copyR = size(tracksXP,1);
copyC = size(tracksXP,2);
removeGapParams.Delimeter = Inf;
removeGapParams.RemoveOtherGaps = true;
removeGapParams.Color = 'k';
removeGapParams.LineStyle = ':';
lineWithGaps(tracksXP,tracksYP,removeGapParams);

% Now go through all present diffusion types and plot one long line that 
% includes positions from all tracks with gaps between them
for k = 1:length(diffTypes)
 	
    switch diffTypes(k)
        case -1 %Unclassified
            ind = find(nanBase(timeRange(1):timeRange(2),:) ==-1);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
        
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 0]);
            end %Allow to change to white if image involved
           %asym linear and unclassified
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
                if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                    for i=1:numTimePlot-1
                        validData=~all(isnan(copyX(i:i+1,:)),1);
                        lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 1]);%[1 1 0][0 102/255 51/255]
                    end 
                end
           %symmetric and unclassified
            indA = find(asymBase ==0);
            indAF  = intersect(indA,ind);
                if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                    for i=1:numTimePlot-1
                        validData=~all(isnan(copyX(i:i+1,:)),1);
                        lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0.6 0 1]);
                    end 
                end
        case 0 %Immobile
            ind = find(base(timeRange(1):timeRange(2),:,1) ==0);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
 	
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0.5 0.3 0]);
            end

        case 1 %Confined
 	
            ind = find(base(timeRange(1):timeRange(2),:,2) ==1);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);

            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 1]);%Real
            end
            
            %Asym linear and confined
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
            if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                for i=1:numTimePlot-1
                    validData=~all(isnan(copyX(i:i+1,:)),1);
                    lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[1 0.7 0]);
                end 
            end
            
        case 2 %Free
 	
            ind = find(base(timeRange(1):timeRange(2),:,3) ==2);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
 	
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color','c');%Real
            end
            
           %asym linear and free
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
            if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                for i=1:numTimePlot-1
                    validData=~all(isnan(copyX(i:i+1,:)),1);
                    lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[1 0 0]);
                end 
            end
            
        case 3 %Directed
            ind = find(base(timeRange(1):timeRange(2),:,4) ==3);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
 	
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color','m');
            end
            
 	           %asym linear and super
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
            if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                for i=1:numTimePlot-1
                    validData=~all(isnan(copyX(i:i+1,:)),1);
                    lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 1 0]);
                end 
            end
    end
 	
end

%show merges and splits
if mergeSplit

    %go over all tracks
    for iTrack = 1 : numTracks

        %parse sequence of events of this compound track and find merges and
        %splits
        seqOfEvents = inputStructure(iTrack).seqOfEvents;
        indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
        indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';

        %go over all splits
        for iSplit = indxSplit

            %get time of splitting
            timeSplit = seqOfEvents(iSplit,1);

            %determine row where starting track is located
            rowS = trackStartRow(iTrack) + seqOfEvents(iSplit,3) - 1;

            %determine row where splitting track is located
            rowSp = trackStartRow(iTrack) + seqOfEvents(iSplit,4) - 1;

            %plot split as a dash-dotted line
            plot([tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp)], ...
                [tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp)],'k-.');

        end

        %go over all merges
        for iMerge = indxMerge

            %get time of merging
            timeMerge = seqOfEvents(iMerge,1);

            %determine row where ending track is located
            rowE = trackStartRow(iTrack) + seqOfEvents(iMerge,3) - 1;

            %determine row where merging track is located
            rowM = trackStartRow(iTrack) + seqOfEvents(iMerge,4) - 1;

            %plot merge as a dashed line
            plot([tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
                [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],'k--');

        end

    end %(for iTrack = 1 : numTracks)

end %(if mergeSplit)

% %show confinement areas if requested
% if showConf
% 
%     %generate circle to plot
%     theta = (0:pi/10:2*pi); %angle
%     xy = [cos(theta') sin(theta')]; %x and y-coordinates
%     
%     %go over symmetric tracks
%     for iTrack = indxCircle'
%         
%         %plot a circle of radius = confinement radius and centered at the
%         %center of this track
%         circleVal = xy .* trackSegmentConfRad(iTrack,1);
%         plot(trackSegmentCenter(iTrack,1)+circleVal(:,1),...
%             trackSegmentCenter(iTrack,2)+circleVal(:,2),'k');
%         
%     end
%     
%     %go over linear tracks
%     for iTrack = indxRectangle'
%     
%         %get the confinement axes
%         axisPara = trackSegmentPrefDir(iTrack,:);
%         axisPerp = [-axisPara(2) axisPara(1)] * trackSegmentConfRad(iTrack,1);
%         axisPara = axisPara * trackSegmentConfRad(iTrack,2);
%         
%         %find the 4 corners of the confinement rectangle
%         cornerCoord = [-axisPara - axisPerp; -axisPara + axisPerp; ...
%             axisPara + axisPerp; axisPara - axisPerp; -axisPara - axisPerp] ...
%             + repmat(trackSegmentCenter(iTrack,:),5,1);
%         
%         %plot the rectangle
%         plot(cornerCoord(:,1),cornerCoord(:,2),'k');
% 
%     end
%     
% end

%%%%% ~~ the end ~~ %%%%%

