function [transDiffAnalysisRes,errFlag] = basicTransientDiffusionAnalysisv1(tracks,...
    probDim,plotRes,peakAlpha)
%BASICTRANSIENTDIFFUSIONANALYSISV1 detects potential diffusion segments of a track and performs MSS analysis on these segments
%
%SYNOPSIS [transDiffAnalysisRes,errFlag] = basicTransientDiffusionAnalysisv1(tracks,...
%     probDim,plotRes,peakAlpha)
%
%INPUT  tracks: -- EITHER --
%                           Matrix indicating the positions and amplitudes
%                           of the tracked features to be plotted. Number
%                           of rows = number of tracks, while number of
%                           columns = 8*number of time points. Each row
%                           consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           NaN is used to indicate time points where a
%                           track does not exist.
%                           -- OR --
%                           Structure array with number of entries equal to
%                           number of compound tracks. Contains the fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of individual tracks in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where individual tracks
%                              do not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of individual track, 2 - end of individual track;
%                              3rd: Index of individual track that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, and number is the index
%                                   of individual track being merged with
%                                   or split from.
%       probDim     : Problem dimensionality.
%                     Optional. Default: 2.
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0.
%                     Results can be plotted only if problem is 2D.
%                     color-coding:
%                     *brown: immobile
%                     *blue: confined diffusion.
%                     *cyan: normal diffusion.
%                     *magenta: super diffusion.
%                     *black: unclassified.
%
%
%       peakAlpha   : confidence level for choosing peaks when initially
%                     segmenting track. Default : 95
%
%OUTPUT transDiffAnalysisRes : A structure array with number of rows = number
%                     of tracks (or compound tracks in the case of merging
%                     and/or splitting). It contains the field
%                     ".segmentClass," with number of rows = number of
%                     individual tracks in the compound track (which merge
%                     with and/or split from each other). In case of no
%                     merging or splitting, each compound track is one
%                     individual track and segmentClass will contain one
%                     row only. segmentClass itself contains the following
%                     fields:
%           .momentScalingSpectrum: (Number of classification
%                     subparts)-by-(20+probDim) matrix, where each row
%                     contains the information for one classification
%                     subpart, and the columns store the following:
%                     (1) Start frame of subpart.
%                     (2) End frame of subpart.
%                     (3) Classification of subpart: 1 = confined, 2 =
%                         free, 3 = directed.
%                     (4) MSS slope resulting in classification.
%                     (5-11) Generalized diffusion coefficients for moment
%                            orders 0-6.
%                     (12-18) Scaling power for moment orders 0-6.
%                     (19) Normal diffusion coefficient (from the MSD).
%                     (20) Confinement radius or localization precision, if subpart is classified as
%                          confined or immobile (NaN otherwise).
%                     (21/22/23) Coordinates of segment  center (1D,2D,3D.
%                     depending on probDim)
%           .momentScalingSpectrum1D: NOT IMPLEMENTED RIGHT NOW.
%           .asymmetry: NOT IMPLEMENTED RIGHT NOW.
%
%       errFlag         : 0 if executed normally, 1 otherwise.
%
% Adapted from rolling window approach, trackTransientDiffusionAnalysis1 Khuloud Jaqaman 2008
%Tony Vega, July 2016
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
transDiffAnalysisRes = [];
errFlag = 0;

%% Input
checkAsym = 0; %Needs work so not imposed currently
%check whether tracks were input
if nargin < 1
    disp('--basicTransientDiffusionAnalysisv1: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end


if nargin < 3 || isempty(plotRes)
    plotRes = 0;
elseif plotRes == 1 && probDim ~= 2
    disp('--basicTransientDiffusionAnalysisv1: Cannot plot tracks if problem is not 2D!');
    plotRes = 0;
end

if nargin < 4 || isempty(peakAlpha)
    peakAlpha = 95;
end


if errFlag
    disp('--trackTransientDiffusionAnalysis1: Please fix input variables');
    return
end
% Load information for probability of misclassification for all adjacent
%diffusion distributions, see paper.
    p = mfilename('fullpath');
    load(strcat(p(1:end-33),'positionConfidenceCI.mat'))
    load(strcat(p(1:end-33),'positionConfidenceFC.mat'))
    load(strcat(p(1:end-33),'positionConfidenceDF.mat'))
    
    if checkAsym
        disp('Sorry - not implemented yet!')
        errFlag = 1;
%         load(strcat(p(1:end-33),'positionConfidenceCI_1D.mat'))
%         load(strcat(p(1:end-33),'positionConfidenceFC_1D.mat'))
%         load(strcat(p(1:end-33),'positionConfidenceDF_1D.mat'))  
    end
%define window sizes
windowAsym = 5;
windowMSS =11;
windowMSSMin = 20;
halfWindowMSS = (windowMSS - 1) / 2;

%Set other variables
alphaAsym = 0.05;
minDuration = 5;
confRadMin = 0;

%specify MSS analysis moment orders
momentOrders = 0:6;

%% Track extraction for analysis

%store input tracks in a new variable
tracksInput = tracks;

%extract segments for analysis if tracks were input as a structure that
%might contain merges and splits
%the point is to reduce compound tracks that contain merges and splits into
%simple separate tracks
%thus this step is not necessary if the tracks were input as a matrix,
%which by definition does not contain unresolved compound tracks.
if isstruct(tracks)

    %get number of input tracks from structure
    numInputTracks = length(tracksInput);

    clear tracks

    [tracks,dummy,compTrackStartRow,numSegments] = ...
                convStruct2MatIgnoreMS(tracksInput);


else

    %get number of input tracks from matrix
    numInputTracks = size(tracksInput,1);

    %indicate rows where tracks start (trivial in this case)
    compTrackStartRow = (1 : numInputTracks)';

    %indicate number of segments in each track (1 for all tracks)
    numSegments = ones(numInputTracks,1);

end

%get number of track segments to be analyzed
numTrackSegments = size(tracks,1);

%get track segment start times, end times and life times
trackSEL = getTrackSEL(tracks);

%find track segments that are long enough for analysis
if checkAsym
    indx4analysis = find(trackSEL(:,3) >= windowAsym);
else
    indx4analysis = find(trackSEL(:,3) >= windowMSSMin);
end


indxNot4analysis = setdiff((1:numTrackSegments)',indx4analysis);

%reserve memory
trackSegmentClassRes = repmat(struct('asymmetry',NaN(1,3),...
    'momentScalingSpectrum',NaN(1,20+probDim),...
    'momentScalingSpectrum1D',NaN(1,20+probDim)),...
    numTrackSegments,1);

%% Step 1a. Initial track segmentation
% Rolling window of maximum pairwise displacement (MPD)to divide track segment into parts with potentially different diffusion behavior

%Reserve memory
gaussDeriv = cell(numTrackSegments,1);

%go over all analyzable track segments
for iTrack = indx4analysis'

    %get track segment start, end and life times
    trackSELCurrent = trackSEL(iTrack,:);

    %find the length of each part
    trackPartLength = trackSELCurrent(1,3);

    %get starting point of this part
    trackPartStart = trackSELCurrent(1,1);

    %get number of MSS analysis rolling windows in this part
    numRollWindows = trackPartLength - windowMSS + 1;

    %if number of rolling windows is larger than the minimum
    %required duration of a classification, proceed with rolling
    %window analysis
    if numRollWindows > minDuration

        %initialize max displacement vector
        maxDisplacement = NaN(numRollWindows,1); 
        %go over all windows and calculate MPD
        for iWindow = 1 : numRollWindows

            %get window start and end points
            startPoint = trackPartStart + iWindow - 1;
            endPoint = startPoint + windowMSS - 1;

            test = tracks(iTrack,8*(startPoint-1)+1:8*endPoint);
            xTest = test(1:8:end);
            yTest = test(2:8:end);
            X =[xTest',yTest'];
            D = pdist(X,'euclidean');
            maxDisplacement(iWindow) = max(D);



        end %(for iWindow = 1 : numRollWindows)
        
        %Smooth MPD vector with Gaussian
        h =1;
        y = maxDisplacement;
        [out] = filterGauss1D(y, 2, 'symmetric');
        %Calculate derivative
        der = diff(out)/h;
        %Get absolute value and normalize curve to put different mobility 
        %switches (e.g. immobile to confined vs. immobile to free) on the 
        %same scale to facilitate their detection. 
        absDer = abs(der);
        normFactor = ((maxDisplacement(2:end)+maxDisplacement(1:end-1))./2);
        %If part of normFactor is zero, set to NaN to avoid error. This
        %should be okay since we are only interested in regions of changes
        %in MPD i.e opposite of where normFactor would be zero
        normFactor(normFactor==0) = NaN;
        normDer = absDer./normFactor;
        % Store
        gaussDeriv{iTrack} = normDer;

    else
        %If track doesn't last for minimum duration, set value to empty
        gaussDeriv{iTrack} = [];
    end
end

%% Step 1b and Step 2: Segmenation Phase 2 and segment classification
% Once again go through tracks long enough for classification
for iTrack = indx4analysis'
    % Separate tracks that change and those that don't
    trackFull = gaussDeriv{iTrack};
%     trackFull = [];
    if ~isempty(trackFull)
        % If track has MPD vector, take maxima above threshold as initial
        % guesses for segments
        level = prctile(trackFull,peakAlpha);

        %% Asymmetry detection
        %Check track to see if any sections have asymmetry
        %this classification scheme is taken from Huet et al (BJ 2006)
        %it classifies tracks as asymmetric or not, based on the scatter of
        %positions along them
        % Notice: Currently not implemented! 
        
        if checkAsym       
            %Divide track into segments using asym minimum
            [segPointsA,peakThresh] = findCloseSegments(trackFull,level,halfWindowMSS,windowAsym,windowMSSMin);
            [segPointsA] = findDiffSegments(segPointsA,windowAsym,peakThresh);% using halfWindowMSS to look at same part of track analyzed above
            n = 1:length(segPointsA)-1;
            difference = segPointsA(n+1)-segPointsA(n);
            trackSELCurrent = trackSEL(iTrack,:);
            partClassAsym = NaN(length(n),3);
            
            %go over all of these segments and determine if they are
            %asymmetric
            for k = 1:length(segPointsA)-1 
            startPoint = trackSELCurrent(1) +segPointsA(k)-1;
            endPoint  = startPoint+ difference(k)-1;

                %get the particle positions along the track
                coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';%Not implemented yet
                coordXYZ = [coordX coordY coordZ];

                %determine whether the track is sufficiently asymmetric
                [~,asymFlag] = asymDeterm2D3D(coordXYZ(:,1:probDim),alphaAsym);

                %classify track as ...
                %1 = linear, if the asymmetry parameter is larger than the threshold
                %0 = not linear, if the asymmetry parameter is smaller than the
                %threshold
                %otherwise, keep track classification as undetermined
                partClassAsym(k,1) = startPoint;
                partClassAsym(k,2) = endPoint;
                partClassAsym(k,3) = asymFlag;
            end


        %find indices of all tracks classified as asymmetric
        indxAsym = find(partClassAsym(:,3) == 1);
        else
        % Otherwise, do not check for asymmetry and simply choose initial
        % segments and try to extend any short segments
                indxAsym = 0;
                %Determine initial segments and try to extend short
                %segments
                [segPointsA,peakThresh] = findCloseSegments(trackFull,level,halfWindowMSS,5,windowMSSMin);
                %
                [segPointsA,peakThresh] = findDiffSegments(segPointsA,windowMSSMin,peakThresh);% halfWindowMSS,windowMSSMin
        end
        
        % If there are multiple segments, check adjacents segments to see
        % if they should be merged or kept separate. Specifically, compare
        % all pairwise distances of segments using kolmogorov-smirnov test 
                       
        if size(segPointsA,1)>=2
            pointsTmp = segPointsA;
            list = [];
            sel = trackSEL(iTrack,1);
            for seg = 1:size(segPointsA,1)-2
                % Get beginning and end of segments
                    b1 = (pointsTmp(seg))+sel-1;
                    b2 =(pointsTmp(seg+1))+sel-1;
                    e1 = (pointsTmp(seg+1)-1)+sel-1;
                    e2 =(pointsTmp(seg+2)-1)+sel-1;
                % Get pairwise distribution of first segment
                    test1 = tracks(iTrack,8*(b1-1)+1:8*e1);
                    xTest1 = test1(1:8:end);
                    yTest1 = test1(2:8:end);
                    X1 =[xTest1',yTest1'];
                    D1 = pdist(X1,'euclidean');
                % Get pairwise distribution of second segment
                    test2 = tracks(iTrack,8*(b2-1)+1:8*e2);
                    xTest2 = test2(1:8:end);
                    yTest2 = test2(2:8:end);
                    X2 =[xTest2',yTest2'];
                    D2 = pdist(X2,'euclidean');
                % KS test
                    try
                        [~,p] = kstest2(D1,D2);
                    catch
                        %If segments are not long enough for test to work,
                        %set to don't merge
                        p = 0;
                    end
                 % If KS test does not find significance, mark for merging
                    if p > 0.05
                        list = [list;seg+1];
                    end
            end
            % If there are segments to be merged
            if ~isempty(list)
                %Remove separation between segments
                segPointsA(list) = [];
                peakThresh(list-1) = [];
                        if checkAsym
                            %go over all of these segments and determine if they are
                            %asymmetric
                            n = 1:length(segPointsA)-1;
                            difference = segPointsA(n+1)-segPointsA(n);
                            for k = 1:length(segPointsA)-1 
                            startPoint = trackSELCurrent(1) +segPointsA(k)-1;
                            endPoint  = startPoint+ difference(k)-1;

                                %get the particle positions along the track
                                coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';%Not implemented yet
                                coordXYZ = [coordX coordY coordZ];

                                %determine whether the track is sufficiently asymmetric
                                [~,asymFlag] = asymDeterm2D3D(coordXYZ(:,1:probDim),alphaAsym);

                                %classify track as ...
                                %1 = linear, if the asymmetry parameter is larger than the threshold
                                %0 = not linear, if the asymmetry parameter is smaller than the
                                %threshold
                                %otherwise, keep track classification as undetermined
                                partClassAsym(k,1) = startPoint;
                                partClassAsym(k,2) = endPoint;
                                partClassAsym(k,3) = asymFlag;
                            end


                        %find indices of all tracks classified as asymmetric
                        indxAsym = find(partClassAsym(:,3) == 1);
                        end
            end

        end
                %}
        %% Step 2. Initial Segment Classification
        %Initialize all variables that will be saved later
        count = length(segPointsA)-2;


        n = 1:length(segPointsA)-1;
        difference = segPointsA(n+1)-segPointsA(n);
        pointClassMSS = NaN(length(n),1);
        mssSlope = pointClassMSS;
        normDiffCoef = pointClassMSS;
        confRadTmp = pointClassMSS;
        centerTmp = NaN(length(n),probDim);
        genDiffCoef = NaN(length(n),length(momentOrders));
        scalingPower = NaN(length(n),length(momentOrders));
        partClassMSS = NaN(length(n),3);
        
        partClassMSS1D = NaN(length(n),3);
        mssSlope1D = pointClassMSS;
        genDiffCoef1D = NaN(length(n),length(momentOrders));
        scalingPower1D = NaN(length(n),length(momentOrders));
        normDiffCoef1D = pointClassMSS;
        confRadius1D = NaN(length(n),2);
        prefDir = NaN(length(n),2);
        trackCenter = NaN(length(n),2);
        
        %Now go through segments and get diffusion classification.
        trackSELCurrent = trackSEL(iTrack,:);
        
        for k = 1:length(segPointsA)-1 
            startPoint = trackSELCurrent(1) +segPointsA(k)-1;
            % If looking at last segment, correct to include last frame
            if k == length(segPointsA)-1
               endPoint  = startPoint+ difference(k)-1;  
            else
               endPoint  = startPoint+ difference(k);
            end
            
            if ismember(k,indxAsym)    

                    % If segment is asymmetric, run different analysis to
                    % get 1D classification
                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter(k,:),confRadius1D(k,:),prefDir(k,:)] = asymmetricDiffusion(startPoint,endPoint,probDim,tracks,iTrack);

                    %since not all track segments may be linear, put analysis results in their
                    %proper place among all track segments
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;              
                    partClassMSS1D(k,3) = partClassAsym(k,3);
                    partClassMSS(k,1)= startPoint;
                    partClassMSS(k,2)= endPoint;
                    partClassMSS(k,3)= pointClass;
                    mssSlope1D(k) = mssSlopeT;
                    genDiffCoef1D(k,:) = genDiffCoefT;
                    scalingPower1D(k,:) = scalingPowerT;
                    normDiffCoef1D(k) = normDiffCoefT;

            else
                % Run MSS Analysis
                    [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                        tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                        probDim,momentOrders,-0.05);

                % If segment is classified as confined or immobile, get
                % additional information
                    if pointClassMSS(k) == 1
                        %If confined, calculate confinement radius
                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);

                    elseif pointClassMSS(k) == 0
                        %If immobile, calculate localization precision
                        [confRadTmp(k)] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                        [~,centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                    else
                        %                         confRadTmp(k) = NaN;
                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);   
                    end 
                    % Store segment begining, end, and classification
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;
                    partClassMSS(k,1)= startPoint;
                    partClassMSS(k,2)= endPoint;
                    partClassMSS(k,3)= pointClassMSS(k);
            end
                    
        end
        %% Step 3. Final Segmentation and classification
        % Try to merge unclassified segments
             %Check if there are symmetric unclassified segments
                    partClassTmp = partClassMSS;
                    mssSlopeTmp = mssSlope;
                    unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                    peakThreshTmp = peakThresh;
                    %If there's more than one segment and there are
                    %unclassified segments
                    if sum(unclassCheck)>0 && size(partClassMSS,1)>1
                        rnd = find(unclassCheck);
                        check = find(unclassCheck);
                        
                        for un = 1:length(rnd)
                            %Attempt to merge unclassified segments with
                            %adjacent segments
                            [partClassMSS,partClassMSS1D,mssSlope,peakThresh] = mergeUnclassSegments(unclassCheck,check(1),partClassMSS,partClassMSS1D,mssSlopeTmp,peakThresh);
                            
                            %if parts have been merged, redo diffusion analysis on
                            %these parts. Otherwise continue
                            if size(partClassMSS,1) < size(partClassTmp,1)
        
                            [partClassMSST,mssSlopeT,normDiffCoefT,confRadTmpT,centerTmpT,genDiffCoefT,scalingPowerT,...
                                partClassMSS1DT,mssSlope1DT,normDiffCoef1DT,confRadius1DT,trackCenterT,genDiffCoef1DT,scalingPower1DT,prefDirT] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaAsym,confRadMin,checkAsym);
                                
                                % If reclassification results in more
                                % mobile segment, check if probability of
                                % misclassification has improved
                                changeSeg = find((partClassMSST(:,3)-partClassMSS(:,3)) > 0);
                                if changeSeg
                                        if partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==2 %If confined track became free
                                            
                                            if isnan(partClassMSS1D(:,3))
                                                oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            else
                                                oldConfidenceList = positionConfidenceFC_1D{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            end
                                            
                                            if isnan(partClassMSS1D(:,3))
                                                newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            else
                                                newConfidenceList = positionConfidenceFC_1D{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            end
                                            
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                        elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==1 %If immobile track became confined
                                            
                                            oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                        elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==2 %If immobile track became free
                                            oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                        elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==3 %If immobile track became directed
                                            oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceDF{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2); 

                                        elseif partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==3 %If confined track became directed
                                            oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceDF{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                        elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==3 %If free track became directed
                                            oldConfidenceList = positionConfidenceDF{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceDF{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);
                                        end
                                        if newConfidence < oldConfidence %If merge improved confidence, accept.
                                            unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                                            check = find(unclassCheck);
                                            %Accept and make all T's (temp values) into actual value
                                            partClassMSS = partClassMSST;
                                            partClassTmp = partClassMSS;
                                            mssSlope = mssSlopeT;
                                            mssSlopeTmp = mssSlope;
                                            normDiffCoef = normDiffCoefT;
                                            confRadTmp = confRadTmpT;
                                            centerTmp = centerTmpT;
                                            genDiffCoef = genDiffCoefT;
                                            scalingPower = scalingPowerT;
                                            partClassMSS1D =partClassMSS1DT;
                                            mssSlope1D = mssSlope1DT;
                                            normDiffCoef1D =normDiffCoef1DT;
                                            confRadius1D = confRadius1DT;
                                            trackCenter = trackCenterT;
                                            genDiffCoef1D = genDiffCoef1DT;
                                            scalingPower1D = scalingPower1DT;
                                            prefDir = prefDirT;
                                            clear oldConfidence newConfidence
                                        else %If not, get previous classification and redo classification
                                          partClassMSS = partClassTmp;
                                          peakThresh = peakThreshTmp;
                                        [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
                                    partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaAsym,confRadMin,checkAsym);
                                            clear oldConfidence newConfidence
                                        end      

                                else
                                    %If mobility did not increase, automatically accept
                                    %merge
                                    unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                                    check = find(unclassCheck);
                                    %Accept and make all T's into actual value
                                    partClassMSS = partClassMSST;
                                    partClassTmp = partClassMSS;
                                    mssSlope = mssSlopeT;
                                    mssSlopeTmp = mssSlope;
                                    normDiffCoef = normDiffCoefT;
                                    confRadTmp = confRadTmpT;
                                    centerTmp = centerTmpT;
                                    genDiffCoef = genDiffCoefT;
                                    scalingPower = scalingPowerT;
                                    partClassMSS1D =partClassMSS1DT;
                                    mssSlope1D = mssSlope1DT;
                                    normDiffCoef1D =normDiffCoef1DT;
                                    confRadius1D = confRadius1DT;
                                    trackCenter = trackCenterT;
                                    genDiffCoef1D = genDiffCoef1DT;
                                    scalingPower1D = scalingPower1DT;
                                    prefDir = prefDirT;
                                end
                                %}
                            end
                        end
                    end
        
            %Try to merge adjacent segments that have the same classification
            checkSeg = partClassMSS(1:end-1,3)-partClassMSS(2:end,3);
            indicateRev = 0;
            % If identical adjacent segments exist and merging has not already
            % been attempted twice, attempt to merge
            while ~isempty(find(checkSeg==0)) && indicateRev < 2
                %Store previous classification, and MSS slope for later
                %comparison
                            mssSlopeTmp = mssSlope;
                            clear mssSlope
                            partClassTmp = partClassMSS;
                            lifeTimeMSS = (partClassTmp(:,2)-partClassTmp(:,1))+1;
                            pointClassMSSAsymT = partClassMSS1D(:,3);
                            pointClassMSS = partClassMSS(:,3);
                    %Find adjacent segments that have the same 1D or 2D classification        
                    for iSubpart = 1 : length(pointClassMSS)-1
                        iSubpartPlus1 = iSubpart + 1;
                        while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) && ...
                                ( (pointClassMSSAsymT(iSubpart) == pointClassMSSAsymT(iSubpartPlus1)) || ... %same asym
                                (isnan(pointClassMSSAsymT(iSubpart)) && isnan(pointClassMSSAsymT(iSubpartPlus1))) ) )
                            
                            partClassMSS1D(iSubpart,2) = partClassMSS1D(iSubpartPlus1,2);
                            partClassMSS1D(iSubpartPlus1,1) = partClassMSS1D(iSubpart,1);
                            
                            partClassMSS(iSubpart,2) = partClassMSS(iSubpartPlus1,2);
                            partClassMSS(iSubpartPlus1,1) = partClassMSS(iSubpart,1);
                            iSubpartPlus1 = iSubpartPlus1 + 1;
                        end
                    end
                    %If segments will be merged, get weighted MSS slope of
                    %segments being merged. This will be used later to test
                    %if the merge should be accepted. MSS slope is weighted
                    %by segment length
                    [~,uniqueParts] = unique(partClassMSS(:,1));
                    for ms = 1:length(uniqueParts)
                        mssSlope(ms) = sum(lifeTimeMSS(partClassMSS(:,2) == partClassMSS(uniqueParts(ms),2)).*mssSlopeTmp(partClassMSS(:,2) == partClassMSS(uniqueParts(ms),2)))./sum(lifeTimeMSS(partClassMSS(:,2) == partClassMSS(uniqueParts(ms),2)));
                    end
                    partClassMSS = partClassMSS(uniqueParts,:);
                    partClassMSS1D = partClassMSS1D(uniqueParts,:);
                    
                    %if parts have been merged, redo diffusion analysis on
                    %these parts. Otherwise continue
                    if size(partClassMSS,1) < size(partClassTmp,1)
                        [partClassMSST,mssSlopeT,normDiffCoefT,confRadTmpT,centerTmpT,genDiffCoefT,scalingPowerT,...
                        partClassMSS1DT,mssSlope1DT,normDiffCoef1DT,confRadius1DT,trackCenterT,genDiffCoef1DT,scalingPower1DT,prefDirT] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaAsym,confRadMin,checkAsym);
                    
                        %If new classification results in
                        %more mobile class, compare new probability of
                        %misclassification using new MSS slope and weighted
                        %MSS slope of segments being merged
                        changeSegF = find((partClassMSST(:,3)-partClassMSS(:,3)) > 0);
                        badChange = [];
                        badC=1;
                        if changeSegF
                            for cs= 1:length(changeSegF)
                                changeSeg = changeSegF(cs);
                                if partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==2 %If confined track became free
                                    newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    oldConfidenceList = newConfidenceList;
                                    oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                    newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==1 %If immobile track became confined
                                    newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    oldConfidenceList = newConfidenceList;
                                    oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                    newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2); 

                                elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==2 %If immobile track became free
                                        oldConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                        newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==3 %If immobile track became directed
                                    oldConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    newConfidenceList = positionConfidenceDF{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                    newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                elseif partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==3 %If confined track became directed
                                    oldConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    newConfidenceList = positionConfidenceDF{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                    newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);

                                elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==3 %If free track became directed
                                    newConfidenceList = positionConfidenceDF{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                    oldConfidenceList = newConfidenceList;
                                    oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1)*1E4)./1E4==round(mssSlope(changeSeg)*1E4)./1E4,3);
                                    newConfidence = newConfidenceList(round(newConfidenceList(:,1)*1E4)./1E4==round(mssSlopeT(changeSeg)*1E4)./1E4,2);    
                                end
                                %Check if probability of misclassification
                                %increased, mark to be reverted after all
                                %merges have been tested
                                if newConfidence > oldConfidence
                                    badChange(badC)=changeSegF(cs);
                                    badC = badC +1;
                                end
                            end
                            %If there was a bad change, revert it back but
                            %keep good changes
                            if ~isempty(badChange)
                                
                                if size(partClassMSST,1) == length(changeSegF)
                                      partClassMSS = partClassTmp;
                                    [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
                                partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaAsym,confRadMin,checkAsym);
                                clear oldConfidence newConfidence
                                else
                                    %Get segments that will stay new (i.e did not result worse prob. of misclassification)
                                    keep = partClassMSS(~ismember(1:size(partClassMSS,1),badChange),:);
                                    
                                    %Get segments that will revert (i.e resulted in worse prob. of misclassification)
                                    segments = [uniqueParts(2:end);length(partClassTmp)+1]-uniqueParts;
                                    for segS = 1:size(partClassMSS,1)
                                        storage.Segments(segS) = {uniqueParts(segS):uniqueParts(segS)+(segments(segS)-1)};
                                    end
                                    revertSeg = cell2mat(storage.Segments(badChange));
                                    
                                    %Stich everything together and re-run,
                                    %technically we have all the information, but
                                    %this is easiest eay to stich data together
                                    partClassMSS = sortrows([keep;partClassTmp(revertSeg,:)]);
                                        [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
                                partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaAsym,confRadMin,checkAsym);
                                clear oldConfidence newConfidence
                                end
                            else
                                %If merge improved confidence, accept.

                                %Accept and make all T's into actual value
                                partClassMSS = partClassMSST;
                                partClassTmp = partClassMSS;
                                mssSlope = mssSlopeT;
                                normDiffCoef = normDiffCoefT;
                                confRadTmp = confRadTmpT;
                                centerTmp = centerTmpT;
                                genDiffCoef = genDiffCoefT;
                                scalingPower = scalingPowerT;
                                partClassMSS1D =partClassMSS1DT;
                                mssSlope1D = mssSlope1DT;
                                normDiffCoef1D =normDiffCoef1DT;
                                confRadius1D = confRadius1DT;
                                trackCenter = trackCenterT;
                                genDiffCoef1D = genDiffCoef1DT;
                                scalingPower1D = scalingPower1DT;
                                prefDir = prefDirT;
                                clear oldConfidence newConfidence

                            end                            
       
                            indicateRev = indicateRev+1;
                        else
                            %Accept and make all T's into actual value
                            partClassMSS = partClassMSST;
                            mssSlope = mssSlopeT;
                            normDiffCoef = normDiffCoefT;
                            confRadTmp = confRadTmpT;
                            centerTmp = centerTmpT;
                            genDiffCoef = genDiffCoefT;
                            scalingPower = scalingPowerT;
                            partClassMSS1D =partClassMSS1DT;
                            mssSlope1D = mssSlope1DT;
                            normDiffCoef1D =normDiffCoef1DT;
                            confRadius1D = confRadius1DT;
                            trackCenter = trackCenterT;
                            genDiffCoef1D = genDiffCoef1DT;
                            scalingPower1D = scalingPower1DT;
                            prefDir = prefDirT;
                        end
                        %}
                        
                    else
                        checkSeg = 1;
                        continue
                    end
                    checkSeg = partClassMSS(1:end-1,3)-partClassMSS(2:end,3);
            end
            

                  
    else
        %% Analyze whole track 
        trackSELCurrent = trackSEL(iTrack,:);
        
        %Fill in empty 1D data
        mssSlope1D = NaN(1,1);
        mssSlope = NaN(1,1);
        genDiffCoef1D = NaN(1,length(momentOrders));
        genDiffCoef = NaN(1,length(momentOrders));
        scalingPower1D = NaN(1,length(momentOrders));
        scalingPower = NaN(1,length(momentOrders));
        normDiffCoef1D = NaN(1,1);
        normDiffCoef = NaN(1,1);
        confRadius1D = NaN(1,2);
        confRadTmp = NaN;
        centerTmp = NaN(1,probDim);
        
        prefDir = NaN(1,2);
        trackCenter = NaN(1,2);
        partClassMSS1D = [trackSELCurrent(:,1:2) NaN];
        
        if checkAsym ==1

            %get the particle positions along the track
            coordX = (tracks(iTrack,1:8:end))';
            coordY = (tracks(iTrack,2:8:end))';
            coordZ = tracks(iTrack,3:8:end)';
            coordXYZ = [coordX coordY coordZ];

            %determine whether the track is sufficiently asymmetric
            [~,asymFlag] = asymDeterm2D3D(coordXYZ(:,1:probDim),0.05);

            if asymFlag ==1

                [pointClass,mssSlope1D,genDiffCoef1D,scalingPower1D,normDiffCoef1D,trackCenter,confRadius1D,prefDir] = asymmetricDiffusion(trackSELCurrent(:,1),trackSELCurrent(:,2),probDim,tracks,iTrack);

                partClassMSS1D = [trackSELCurrent(:,1:2) asymFlag];
                partClassMSS = [trackSELCurrent(:,1:2) pointClass];

            else
                
                trackSELCurrent = trackSEL(iTrack,:);
                [pointClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = trackMSSAnalysis(...
                    tracks(iTrack,:),...
                    probDim,momentOrders,-0.05);
                %Get confinement info if necessary
                if pointClassMSS == 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,:),probDim,confRadMin);
                elseif pointClassMSS == 0
                          [confRadTmp] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                          [~,centerTmp] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                else
                    confRadTmp = NaN;
                    [~,centerTmp] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                end 
                partClassMSS =[trackSELCurrent(:,1:2) pointClassMSS]; 
                partClassMSS1D = [trackSELCurrent(:,1:2) asymFlag];

            end
        else
                [pointClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = trackMSSAnalysis(...
                    tracks(iTrack,:),...
                    probDim,momentOrders,alphaValues(1));
        %Get confinement info if necessary
                if pointClassMSS == 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,:),probDim,confRadMin);
                elseif pointClassMSS == 0
                          [confRadTmp] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                          [~,centerTmp] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                else
                    confRadTmp = NaN;
                    [~,centerTmp] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                end 
                partClassMSS =[trackSELCurrent(:,1:2) pointClassMSS]; 
                    
        end

    end
    
    %% Store results
    trackClassMSS = [partClassMSS,mssSlope,genDiffCoef,scalingPower, normDiffCoef,confRadTmp,centerTmp];
    trackClassMSS1D = [partClassMSS1D,mssSlope1D,genDiffCoef1D,scalingPower1D, normDiffCoef1D,confRadius1D,trackCenter,prefDir]; 
    trackSegmentClassRes(iTrack).asymmetry = partClassMSS1D;
    trackSegmentClassRes(iTrack).momentScalingSpectrum = trackClassMSS;
    trackSegmentClassRes(iTrack).momentScalingSpectrum1D = trackClassMSS1D;
    
end
%% Store trivial nonclassification information for tracks that are not classifiable
if ~isempty(indxNot4analysis')
for iTrack = indxNot4analysis' 
    trackSELCurrent = trackSEL(iTrack,:);
    trackSegmentClassRes(iTrack).asymmetry(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).momentScalingSpectrum(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).momentScalingSpectrum1D(1:2) = trackSELCurrent(1:2);
end
end
%% save results in output structure

%reserve memory
segmentClass = struct('asymmetry',[],'momentScalingSpectrum',[],'momentScalingSpectrum1D',[]);
transDiffAnalysisRes = repmat(struct('segmentClass',segmentClass),numInputTracks,1);

%go over all input tracks
for iTrack = 1 : numInputTracks

    %go over the segments of each track
    for iSegment = 1 : numSegments(iTrack)

        %store the segment's classification results
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetry = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetry;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum1D = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum1D;
%         transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetryAfterMSS = ...
%             trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetryAfterMSS;

    end %(for iSegment = 1 : numSegments(iTrack))

end %(for iTrack = 1 : numInputTracks)

%% plotting

%plot results if requested
if plotRes
    plotTracksTransDiffAnalysis2D(tracksInput,transDiffAnalysisRes,[],1,[]);
end



end

%% Subfunctions
function [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaAsym,confRadMin,checkAsym)                 
% Description: Same analysis at Step 2 in main code (lines: 424-503), just created subfunction to save space since used frequently    
    numSeg = size(partClassMSS,1);
    pointClassMSS = NaN(numSeg,1);
    mssSlope = pointClassMSS;
    normDiffCoef = pointClassMSS;
    confRadTmp = pointClassMSS;
    centerTmp = NaN(numSeg,probDim);
    genDiffCoef = NaN(numSeg,length(momentOrders));
    scalingPower = NaN(numSeg,length(momentOrders));

    partClassMSS1D = NaN(size(partClassMSS,1),3);
    mssSlope1D = pointClassMSS;
    genDiffCoef1D = NaN(numSeg,length(momentOrders));
    scalingPower1D = NaN(numSeg,length(momentOrders));
    normDiffCoef1D = pointClassMSS;
    confRadius1D = NaN(numSeg,2);
    prefDir = NaN(numSeg,2);
    trackCenter = NaN(numSeg,2);
    asymFlag = 0;
    
    for k = 1:numSeg 
    startPoint = partClassMSS(k,1);
    endPoint  = partClassMSS(k,2);
                if checkAsym
                %get the particle positions along the track
                coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';
                coordXY = [coordX coordY coordZ];

                %determine whether the track is sufficiently asymmetric

                    [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);
                end
                %classify track as ...
                %1 = linear, if the asymmetry parameter is larger than the threshold
                %0 = not linear, if the asymmetry parameter is smaller than the
                %threshold
                %otherwise, keep track classification as undetermined

                if asymFlag ==1

                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter(k,:),confRadius1D(k,:),prefDir(k,:)] = asymmetricDiffusion(startPoint,endPoint,probDim,tracks,iTrack);

                    %since not all track segments are linear, put analysis results in their
                    %proper place among all track segment
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;              
                    partClassMSS1D(k,3) = asymFlag;
                    partClassMSS(k,1)= startPoint;
                    partClassMSS(k,2)= endPoint;
                    partClassMSS(k,3)= pointClass;

                    mssSlope1D(k) = mssSlopeT;
                    genDiffCoef1D(k,:) = genDiffCoefT;
                    scalingPower1D(k,:) = scalingPowerT;
                    normDiffCoef1D(k) = normDiffCoefT;
                else

                    [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                        tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                        probDim,momentOrders,-0.05);
                    
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;
                    
                    if pointClassMSS(k) == 1  
                    [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                    elseif pointClassMSS(k) == 0
                          [confRadTmp(k)] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                          [~,centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                    else
                        confRadTmp(k) = NaN;
                        [~,centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                    end
                    partClassMSS(k,3) = pointClassMSS(k);
                end
    end

end

function  [partClassMSS,partClassMSS1D,mssSlope,peakThresh] = mergeUnclassSegments(unclassCheck,check,partClassMSS,partClassMSS1D,mssSlopeTmp,peakThresh)
% Description: Goes through all unclassified segments and attempts to merge
% with adjacent segments
        mssSlope =mssSlopeTmp;
        if size(partClassMSS,1)>1
        %If short segment is first, connect to succeeding
            if check(1) == 1 
                partClassMSS(1,2) = partClassMSS(2,2);
                partClassMSS1D(1,2) = partClassMSS1D(2,2);
                partClassMSS(1,3) = partClassMSS(2,3);
                partClassMSS1D(1,3) = partClassMSS1D(2,3);
                partClassMSS(2,:) =[];
                partClassMSS1D(2,:) =[];
                mssSlope(1) = mssSlope(2);
                mssSlope(2) = [];
                peakThresh(1) = [];
             %If segment is last, connect to preceding
            elseif check(1) == length(unclassCheck) %changed from end
                partClassMSS(end-1,2) = partClassMSS(end,2);
                partClassMSS1D(end-1,2) = partClassMSS1D(end,2);
                partClassMSS(end,:) =[];
                partClassMSS1D(end,:) =[];
                mssSlope(end) = mssSlope(end-1);
                mssSlope(end-1) = [];
                peakThresh(end) = [];
                
            %If segment is somewhere in middle
            elseif ~isempty(check)
                %If preceeding and succeeding segments are the same
                %diffusion, merge everything
                if partClassMSS(check(1)+1,3) == partClassMSS(check(1)-1,3) && isnan(partClassMSS1D(check(1)+1,3)) == isnan(partClassMSS1D(check(1)-1,3))
                    partClassMSS(check(1)-1,2) = partClassMSS(check(1)+1,2);
                    partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1)+1,2);
                    partClassMSS(check(1)+1,:) =[];
                    partClassMSS1D(check(1)+1,:) =[];
                    partClassMSS(check(1),:) =[];
                    partClassMSS1D(check(1),:) =[];
                    mssSlope(check(1)-1) = max(mssSlope(check(1):check(1)+1));
                    mssSlope(check(1)+1) =[];
                    mssSlope(check(1)) =[];
                    peakThresh(check(1))=[];
                    peakThresh(check(1)-1)=[];                    
                elseif partClassMSS1D(check(1)+1,3) == partClassMSS1D(check(1)-1,3) && isnan(partClassMSS(check(1)+1,3)) == isnan(partClassMSS(check(1)-1,3))
                    partClassMSS(check(1)-1,2) = partClassMSS(check(1)+1,2);
                    partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1)+1,2);
                    partClassMSS(check(1)+1,:) =[];
                    partClassMSS1D(check(1)+1,:) =[];
                    partClassMSS(check(1),:) =[];
                    partClassMSS1D(check(1),:) =[];  
                    peakThresh(check(1))=[];
                    peakThresh(check(1)-1)=[];                    
                 %Otherwise score the segments and merge to weaker score
                else
                    if peakThresh(check(1)) > peakThresh(check(1)-1)
                        %if higher peak is in front, delete peak behind
                        partClassMSS(check(1)-1,2) = partClassMSS(check(1),2);
                        partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1),2);
                        partClassMSS(check(1),:) =[];
                        partClassMSS1D(check(1),:) =[];
                        mssSlope(check(1))=[];
                        peakThresh(check(1)-1)=[];
                    else
                        partClassMSS(check(1)+1,1) = partClassMSS(check(1),1);
                        partClassMSS1D(check(1)+1,1) = partClassMSS1D(check(1),1);
                        partClassMSS(check(1),:) =[];
                        partClassMSS1D(check(1),:) =[];
                        mssSlope(check(1))=[];
                        peakThresh(check(1))=[];
                    end
                end
                       
            end
        end
end
function [pointClass, mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter,confRadius1D,prefDir] = asymmetricDiffusion(startPoint,endPoint,probDim,tracks, iTrack)
% Description: Classifies segment based on 1d diffusion
    trackCoordX = tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint)';
    deltaCoordX = tracks(iTrack,8*(startPoint-1)+5:8:8*endPoint)';
    trackCoordY = tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint)';
    deltaCoordY = tracks(iTrack,8*(startPoint-1)+6:8:8*endPoint)';
    trackCoordZ = tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint)';
    deltaCoordZ = tracks(iTrack,8*(startPoint-1)+7:8:8*endPoint)';
    trackCoord = [trackCoordX trackCoordY trackCoordZ];
    deltaCoord = [deltaCoordX deltaCoordY deltaCoordZ];
    trackCoord = trackCoord(:,1:probDim);
    deltaCoord = deltaCoord(:,1:probDim);

    %project onto direction of motion
    [posAlongDir,deltaPosAlongDir] = projectCoordOntoDir(trackCoord,...
        deltaCoord,[],[]);

    %construct matrix of linear tracks with projected positions
    trackCoord2 = [posAlongDir zeros(length(posAlongDir),3) deltaPosAlongDir zeros(length(posAlongDir),3)]';
    trackCoord2 = trackCoord2(:)';
    momentOrders = 0 : 6;
    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT] = ...
    trackMSSAnalysis(trackCoord2,1,momentOrders,-0.05);

%% Confinement of asym 
    xyzCoord = [trackCoordX trackCoordY trackCoordZ];

    %find the eignevalues of the variance-covariance matrix of this track's
    %positions
    [eigenVec,eigenVal] = eig(nancov(xyzCoord(:,1:probDim)));
    eigenVal = diag(eigenVal);

    %calculate the confinement radius along the preferred direction of
    %motion
    confRadius1D(1,2) = sqrt( max(eigenVal) * (3) );

    %calculate the confinement radius perpendicular to the preferred
    %direction of motion
    confRadius1D(1,1) = sqrt( mean(eigenVal(eigenVal~=max(eigenVal))) * (probDim + 1) );

    %calculate the track's center
    trackCenter = nanmean(xyzCoord(:,1:probDim));

    %store the preferred direction of motion
    prefDir = eigenVec(:,eigenVal==max(eigenVal))';
end

function [confRadTmp,centerTmp] = estimConfRad(tracks,probDim,confRadMin)

    %get subpart's coordinates
    xCoord = tracks(1:8:end)';
    yCoord = tracks(2:8:end)';
    zCoord = tracks(3:8:end)';
    xyzCoord = [xCoord yCoord zCoord];
    
    if all(isnan(xyzCoord))
        
        confRadTmp = NaN;
        centerTmp = NaN;
        
    else

        %find the eigenvalues and eigenvectors of the variance-covariance
        %matrix of this track's positions
        eigenVal = eig(nancov(xyzCoord(:,1:probDim)));
        
        %calculate the track's confinement radius
        if confRadMin
            confRadTmp = sqrt( min(eigenVal) * (probDim + 2) );
        else
            confRadTmp = sqrt( mean(eigenVal) * (probDim + 2) );
        end
        
        %calculate the track's center
        centerTmp = nanmean(xyzCoord(:,1:probDim));
        
    end
    
end

function [locPreTmp] = estimLocPre(tracks,probDim)

    %get subpart's coordinates
    xCoord = tracks(1:8:end)';
    yCoord = tracks(2:8:end)';
    zCoord = tracks(3:8:end)';
    xyzCoord = [xCoord yCoord zCoord];

    %find the eigenvalues and eigenvectors of the variance-covariance
    %matrix of this track's positions
    eigenVal = eig(nancov(xyzCoord(:,1:probDim)));

    %calculate the track's confinement radius
        locPreTmp = sqrt( mean(eigenVal) );

end

function [segPoints,peakThresh] = findCloseSegments(trackFull,level,halfWindowMSS,absMin,windowMSSMin)
%Description: Function finds initial diffusion segments in track and attempts to
%Input:
% trackFull - normalized MPD derivative
% 
% level - percentile threshold for peaks found in trackFull
%
% halfWindowMSS - half the size of the rolling window used
%
% absMin - shortest a segment is allowed to be (can still be classified as asymmetric)
%
% windowMSSMin - size of the rolling window used
%
%Output
%
% segPoints - locations of the points that separate segments in a track
%
% peakThresh- peaks that were above the percentile threshold in trackFull
%
%extend any short segments (15< length <20) 
        %Find all peaks in normalized MPD derivative
        [peaks,locs] = findpeaks(trackFull);      
        offLim = [16:19]; %
        %Find peaks that are above threshold
        peakThresh = peaks(peaks >= level);
        % Make initial segments based on these peaks
        segPoints = [1;locs(peaks >= level);length(trackFull)];
        segPoints(2:end)=segPoints(2:end)+1+halfWindowMSS; %half of windowMin, 1 for deriv offset
        segPoints(end)=segPoints(end)+1+halfWindowMSS; %half of windowMin, 1 for edn correction
        n = 1:length(segPoints)-1;
        
        %Check that these segments are long enough
        difference = (segPoints(n+1)-segPoints(n))+1;
        checkT = find(difference < windowMSSMin);
        checkM = find(difference >= 16);
        indMin = ismember(checkT,checkM);
        check = checkT(indMin);
        %If there are segments that are too short, extend into adjacent
        %segments. Segments must be within 16:19 frames long and can't 
        %shorten adjacent segments to these lengths 
    if ~isempty(check) 
        for m = 1:length(check) 
        %If short segment is first, connect to succeeding
            if check(m) ==1 && difference(check(m)) >= 16
                comp1 = difference(check(m));
                %Verify that extenstion won't result in either 1) adjacent
                %segment being less than absMin (if already < 20) or 2)
                %adjacent segment length falling into 16:19 length (if
                %>=20)
                if difference(2) -(windowMSSMin-comp1) >= absMin && ~ismember(difference(2) -(windowMSSMin-comp1),offLim)
                    difference(1) = comp1 + (windowMSSMin-comp1);
                    difference(2) = difference(2) -(windowMSSMin-comp1);

                    % Keep peak but move segPoint
                    segPoints(2,:) = segPoints(2,:) + (windowMSSMin-comp1);
                    continue
                end
    
             %If segment is last, connect to preceding
            elseif check(m) == length(difference) && difference(check(m)) >= 16
                comp1 = difference(check(m));
                if difference(check(m)-1) -(windowMSSMin-comp1) >= absMin && ~ismember(difference(check(m)-1) -(windowMSSMin-comp1),offLim)
                    difference(check(m)) = comp1 + (windowMSSMin-comp1);
                    difference(check(m)-1) = difference(check(m)-1) -(windowMSSMin-comp1);
                    % Keep peak but move segPoint
                    segPoints(end-1,:) = segPoints(end-1,:) - (windowMSSMin-comp1);
                    continue
                end
                
            

            %If segment is somewhere in middle, connect to segment with lowest
            %peak
            elseif ~isempty(check) && difference(check(m)) >= 16
                if peakThresh(check(m)) > peakThresh(check(m)-1)
                    %if higher peak is in front, delete peak behind, same
                    %as if it were last segment
                    comp1 = difference(check(m));
                    if difference(check(m)-1) -(windowMSSMin-comp1) >= absMin && ~ismember(difference(check(m)-1) -(windowMSSMin-comp1),offLim)
                        difference(check(m)) = comp1 + (windowMSSMin-comp1);
                        difference(check(m)-1) = difference(check(m)-1) -(windowMSSMin-comp1);
                        
                        % Keep peak but move segPoint
                        segPoints(check(m),:) = segPoints(check(m),:) - (windowMSSMin-comp1);
                        continue
                    end
                else
                    %if higher peak is behind, delete peak in front, same
                    %as if it were first segment
                        comp1 = difference(check(m));
                        if difference(check(m)+1) -(windowMSSMin-comp1) >= absMin && ~ismember(difference(check(m)+1) -(windowMSSMin-comp1),offLim)
                            difference(check(m)) = comp1 + (windowMSSMin-comp1);
                            difference(check(m)+1) = difference(check(m)+1) -(windowMSSMin-comp1);

                            % Keep peak but move segPoint
                            segPoints(check(m)+1,:) = segPoints(check(m)+1,:) + (windowMSSMin-comp1);
                            continue

                        end
                end
                       
            end
        end
    end
end

function [segPoints,peakThresh] = findDiffSegments(segPoints,windowMSSMin, peakThresh)
% Description: Check to see if at least one segment is classfiable. If not,
% attempt to extend a segment
    n = 1:length(segPoints)-1;
    difference = (segPoints(n+1)-segPoints(n))+1;
    checkT = find(difference < windowMSSMin); %Are there segments below the the classifiable length, windowMSSMin
    indMin = difference(difference >= windowMSSMin); %Is at least one segment classifiable?
    ind = difference(difference < windowMSSMin); %Where are the segments below the the classifiable length, windowMSSMin
    %Select shortest segment found 
    check = checkT(ind  == min(ind));

    %While there are short segments and there is not one that is
    %classifiable, attempt to merge segments until at least one is
    %classifiable
    while ~isempty(check) && isempty(indMin)
    %If short segment is first, connect to succeeding
        if check(1) == 1 && difference(check(1)) <16 %And less then a minimum length
            segPoints(2,:) =[];
            peakThresh(1) =[];

            n = 1:length(segPoints)-1;
            difference = segPoints(n+1)-segPoints(n);
            checkT = find(difference < windowMSSMin);
            ind = difference(difference < windowMSSMin);
            indMin = difference(difference >= windowMSSMin);
            check = checkT(ind  == min(ind));
            continue
        elseif check(1) ==1 && difference(check(1)) >= 16
            comp1 = difference(check(1));
            if difference(2) -(windowMSSMin-comp1) >= windowMSSMin
                difference(1) = comp1 + (windowMSSMin-comp1);
                difference(2) = difference(2) -(windowMSSMin-comp1);
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
                % Keep peak but move segPoint
                segPoints(2,:) = segPoints(2,:) + (windowMSSMin-comp1);
                continue
            else
                segPoints(2,:) =[];
                peakThresh(1) =[];

                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
                continue
            end

        end
         %If segment is last, connect to preceding
        if check(1) == length(difference)  && difference(check(1)) <16
            segPoints(end-1,:) = [];
            peakThresh(end) = [];

            n = 1:length(segPoints)-1;
            difference = segPoints(n+1)-segPoints(n);
            checkT = find(difference < windowMSSMin);
            ind = difference(difference < windowMSSMin);
            indMin = difference(difference >= windowMSSMin);
            check = checkT(ind  == min(ind));
            continue
        elseif check(1) == length(difference) && difference(check(1)) >= 16
            comp1 = difference(check(1));
            if difference(check(1)-1) -(windowMSSMin-comp1) >= windowMSSMin
                difference(check(1)) = comp1 + (windowMSSMin-comp1);
                difference(check(1)-1) = difference(check(1)-1) -(windowMSSMin-comp1);
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
                % Keep peak but move segPoint
                segPoints(end-1,:) = segPoints(end-1,:) - (windowMSSMin-comp1);
                continue
            else
                segPoints(end-1,:) = [];
                peakThresh(end) = [];
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
                continue
            end

        end

        %If segment is somewhere in middle, connect to segment with lowest
        %peak
        if ~isempty(check)
            if peakThresh(check(1)) > peakThresh(check(1)-1)
                segPoints(check(1),:) =[]; %if higher peak is in front, delete peak behind
                peakThresh(check(1)-1)=[];
            else
                segPoints(check(1)+1,:) =[];
                peakThresh(check(1)) =[];
            end

            n = 1:length(segPoints)-1;
            difference = segPoints(n+1)-segPoints(n);
            checkT = find(difference < windowMSSMin);
            ind = difference(difference < windowMSSMin);
            indMin = difference(difference >= windowMSSMin);
            check = checkT(ind  == min(ind));
        end
    end
end
                    