%% Identify inner and outer points (reference vertices and points inside) - with vertex-to-nearest-neighbour to centreFaces of the Normals

% Input:    pelvisNum: Numeric identifier used only for logging
%           shrink: Structure that already contains fields produced by shrinkInside
%           importData: Combined mesh - mesh of pelvis defect and refined patches 
%           refData: Reference pelvis mesh (vertices,faces)

% Output:   shrink: struct with results added under shrink.inside

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function shrink = inoutPoints(shrink,pelvisNum,importData,refData)

refFaces = refData.faces;
% Combined mesh: mesh of pelvis defect and refined patches (without acetabulum patch)
centreFaces = importData.comCentreFaces;
normals = importData.comNormals;
verticesMesh = importData.comVertices; %%%
facesMesh = importData.comFaces; %%%

% Vector: point to nearest centreFaces -> centreFaces Normal
verticesRef = centreFaces;
verticesInside = shrink.inside.refVerticesInside;
pointsInside = shrink.inside.refPointsInside;
% Create a KD-tree searcher for verticesB
MdlB = KDTreeSearcher(verticesRef);
% Find the nearest neighbor in verticesB for each vertex in verticesA
K = 1;
[idxVertices, distVertices] = knnsearch(MdlB, verticesInside, 'K', K);
[idxPoints, distPoints] = knnsearch(MdlB, pointsInside, 'K', K);
% Calculate the distances based on the nearest neighbor indices
minDistancesVertices = sqrt(sum((verticesInside - verticesRef(idxVertices,:)).^2, 2));
minDistancesPoints = sqrt(sum((pointsInside - verticesRef(idxPoints,:)).^2, 2));

% Initialize
insideOutsideVertices = false(size(verticesInside, 1), 1);
% Loop through each vertex
for i = 1:size(verticesInside, 1)
    % Get the nearest reference point
    nearestRefVertice = verticesRef(idxVertices(i),:);
    % Calculate the vector from the nearest normal origin to the point
    vectorNormalOriginToVertice = verticesInside(i,:) - nearestRefVertice;
    % Get the corresponding normal vector
    normalVectorVertices = normals(idxVertices(i),:);
    % Calculate the dot product
    dotProductVertice = dot(normalVectorVertices, vectorNormalOriginToVertice);
    % Determine if the point is inside or outside
    if dotProductVertice > 0
        insideOutsideVertices(i) = false;  % Outside the pelvis defect
    else
        insideOutsideVertices(i) = true; % Inside the pelvis defect
    end
end
shrink.inside.inoutVertices = insideOutsideVertices; % Logical mask of reference vertices
outsideVertices = verticesInside(~insideOutsideVertices,:);
shrink.inside.outVertices = outsideVertices;
insideVertices = verticesInside(insideOutsideVertices,:);
shrink.inside.inVertices = insideVertices;
% Logical mask of reference pelvis
shrink.inside.inVerticesMask = false(size(shrink.inside.refVerticesInsideMask,1),1);
shrink.inside.inVerticesIdx = shrink.inside.refVerticesInsideIdx(insideOutsideVertices,:);
shrink.inside.inVerticesMask(shrink.inside.inVerticesIdx) = true;
% Find faces of the vertices
isVertexInside = false(max(refFaces(:)), 1);
isVertexInside(shrink.inside.inVerticesIdx) = true;
allRefFacesInside = all(isVertexInside(refFaces), 2); % Faces with all vertex indices inside
shrink.inside.inFaces = refFaces(allRefFacesInside, :); % Pelvis data (reference)

% Initialize
insideOutsidePoints = false(size(pointsInside, 1), 1);
% Loop through each vertex
for i = 1:size(pointsInside, 1)
    % Get the nearest reference point
    nearestRefPoint = verticesRef(idxPoints(i), :);
    % Calculate the vector from the nearest normal origin to the point
    vectorNormalOriginToPoint = pointsInside(i, :) - nearestRefPoint;
    % Get the corresponding normal vector
    normalVectorPoints = normals(idxPoints(i), :);
    % Calculate the dot product
    dotProductPoint = dot(normalVectorPoints, vectorNormalOriginToPoint);
    % Determine if the point is inside or outside
    if dotProductPoint > 0
        insideOutsidePoints(i) = false;  % Outside the pelvis defect
    else
        insideOutsidePoints(i) = true; % Inside the pelvis defect
    end
end
shrink.inside.inoutPoints = insideOutsidePoints; % Logical mask of reference points
outsidePoints = pointsInside(~insideOutsidePoints,:);
shrink.inside.outPoints = outsidePoints;
insidePoints = pointsInside(insideOutsidePoints,:);
shrink.inside.inPoints = insidePoints;
% Logical mask of point base (gridPoints)
shrink.inside.inPointsMask = false(size(shrink.inside.refPointsInsideMask,1), 1); % Logical mask of point base (gridPoints)
shrink.inside.inPointsIdx = shrink.inside.refPointsInsideIdx(insideOutsidePoints,:);
shrink.inside.inPointsMask(shrink.inside.inPointsIdx) = true; % Idx of points base

disp(['identify ínner/outer points (dot product): pelvis defect ',num2str(pelvisNum)])

end