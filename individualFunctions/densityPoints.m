%% Identify areas with low point density (outliers) - also for clustering

% Input:    pelvisNum: Numeric identifier used only for logging
%           shrink: Structure contain fields from prior steps
%           threshold: percentile in [0,100]
%           refData: Reference pelvis mesh (vertices,faces)

% Output:   shrink: struct with additional fields under shrink.inside

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [shrink] = densityPoints(shrink,pelvisNum,threshold,refData)

refFaces = refData.faces;
% Vertices / Points
insideVerticesPoints = [shrink.inside.inVertices; shrink.inside.inPoints];

% Density
bandwidth = 1; % 0.5 / 0.6 / 1
% Density estimation
tic
density = mvksdensity(insideVerticesPoints, insideVerticesPoints, 'Bandwidth', bandwidth);
shrink.inside.executionTimeDensity = toc;
shrink.inside.pointCloudDensity = density;

% Normalisation
minDensity = min(density);
maxDensity = max(density);
normDensity = (density - minDensity) / (maxDensity - minDensity);
shrink.inside.pointCloudNormDensity = normDensity;
% Filter out points with low point density
thresholdLowDensity = prctile(normDensity,threshold);
lowDensity = normDensity <= thresholdLowDensity;
shrink.inside.inoutDensity = ~lowDensity;
shrink.inside.filteredOutAll = insideVerticesPoints(lowDensity,:);
shrink.inside.filteredInAll = insideVerticesPoints(~lowDensity,:);
% Vertices / Points
numInsideVertices = size(shrink.inside.inVertices,1);
indicesVertices = 1:numInsideVertices;
indicesPoints = numInsideVertices+1 : numInsideVertices+size(shrink.inside.inPoints,1);
shrink.inside.filteredInVertices = insideVerticesPoints(~lowDensity(indicesVertices),:);
shrink.inside.filteredOutVertices = insideVerticesPoints(lowDensity(indicesVertices),:);
shrink.inside.filteredInPoints = shrink.inside.inPoints(~lowDensity(indicesPoints),:);
shrink.inside.filteredOutPoints = shrink.inside.inPoints(lowDensity(indicesPoints),:);
% Idx of reference pelvis (vertices)
shrink.inside.filteredInVerticesMask = false(size(shrink.inside.inVerticesMask,1),1);
shrink.inside.filteredInVerticesIdx = shrink.inside.inVerticesIdx(~lowDensity(indicesVertices),:);
shrink.inside.filteredInVerticesMask(shrink.inside.filteredInVerticesIdx) = true;
% Find faces of the vertices
isVertexInside = false(max(refFaces(:)), 1);
isVertexInside(shrink.inside.filteredInVerticesIdx) = true;
allRefFacesInside = all(isVertexInside(refFaces), 2); % Faces with all vertex indices inside
shrink.inside.filteredInFaces = refFaces(allRefFacesInside, :); % Pelvis data (reference)
% Idx of point base (gridPoints)
shrink.inside.filteredInPointsMask = false(size(shrink.inside.inPointsMask,1), 1); % Logical mask of point base (gridPoints)
shrink.inside.filteredInPointsIdx = shrink.inside.inPointsIdx(~lowDensity(indicesPoints),:);
shrink.inside.filteredInPointsMask(shrink.inside.filteredInPointsIdx) = true;

% Pelvis defect with viridis color
normDensity(normDensity < 0) = 0;   % Clamp values below the range to 0
normDensity(normDensity > 1) = 1;   % Clamp values above the range to 1
% Apply colormap (recommended: viridis)
% 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
colourNum = 32768;
colourMap = viridis(colourNum);
colourIdx = round(normDensity * (colourNum-1)) + 1;
rgbColour = colourMap(colourIdx,:);
shrink.inside.densityRGB = rgbColour;

disp(['identify inner/outer points (density): pelvis defect ',num2str(pelvisNum)])

end