%% alphaShape and find points of reference pelvis inside defect

% Input:    pelvisNum: Numeric identifier used only for logging
%           type: 'all','acentre','acetabulum'
%           shrink: Structure that already contains an alpha-shape result at shrink.(shrinkType) 
%           shrinkType: Field name selecting the hull to use, e.g. 'all', 'acentre', 'acetabulum'.
%           refData: Reference pelvis geometry (vertices,faces)
%           aRadius: Initial alpha radius %
%           addPoints: Extra points appended to the point set
%           refGeometryIdx: Indices in refData.vertices %
%           pointsBoxNum: Total number of points represented by addPointsMaskIdx
%           addPointsMaskIdx: For each row in addPoints, the corresponding index in a global/boxed point list of length pointsBoxNum.

% Output:   shrink: struct with results added under shrink.inside

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = shrinkInside(shrink,pelvisNum,shrinkType,refData,addPoints,refGeometryIdx, ...
    pointsBoxNum,addPointsMaskIdx) % point base (gridPoints)

refVertices = refData.vertices;
refFaces = refData.faces;

% ShrinkWrap (alphaShape) for rough hull
% Use of previous function
%obj = obj.shrinkAlpha(...);

% Find vertices of the reference pelvis inside hull (alphaShape)
FV.faces = shrink.(shrinkType).faces;
FV.vertices = shrink.(shrinkType).allVerticesPoints;
refVerticesInsideLog = inpolyhedron(FV,refVertices);
% Correction:
% Part of the reference geometry part of the alpha shape
% Not all points on the surface recognised as interior points with the inpolyhedron function
refVerticesInsideLog(refGeometryIdx) = true;
refVerticesInside = refVertices(refVerticesInsideLog,:); % Pelvis data (reference)

% Find faces of the vertices
allRefFacesInside = all(refVerticesInsideLog(refFaces), 2); % Faces with all vertex indices inside
refFacesInside = refFaces(allRefFacesInside, :); % Pelvis data (reference)

% Store results in the object
shrink.inside.refVerticesInside = refVerticesInside;
shrink.inside.refVerticesInsideMask = refVerticesInsideLog;
shrink.inside.refVerticesInsideIdx = find(refVerticesInsideLog);
shrink.inside.refFacesInside = refFacesInside;

% Check which of the points inside reference pelvis are inside the alphaShape
addInsideLog = inpolyhedron(FV, addPoints);
addInside = addPoints(addInsideLog, :);
shrink.inside.refPointsInside = addInside;
shrink.inside.refPointsInsideMask = false(pointsBoxNum, 1);
shrink.inside.refPointsInsideIdx = addPointsMaskIdx(addInsideLog,:);
shrink.inside.refPointsInsideMask(shrink.inside.refPointsInsideIdx) = true;

disp(['shrinked mesh and found points inside: pelvis defect ',num2str(pelvisNum)])

end