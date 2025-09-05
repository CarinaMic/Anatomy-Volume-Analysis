%% Thicken mesh (point cloud)

% Input:    pelvisNum: Numeric identifier used only for logging
%           importData:  Mesh sampling providing normals (comCentreFaces, comNormals)
%           thickness: Total inward offset distance to generate
%           step_size: Step size for each offset layer (positive)

% Output:   shrink: struct with additional fields under shrink.inside

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = thickenMesh(pelvisNum, importData, thickness, step_size)

centreFaces = importData.comCentreFaces;
normals = importData.comNormals;
num_steps = thickness / step_size;

% Thicken mesh
verticesNew = [];
for i = 1:num_steps
    % Move vertices inward along the normal direction
    verticesSteps = centreFaces - (i * step_size * normals); % Move points by step_size along normals (inside)
    verticesNew = [verticesNew; verticesSteps];
end

% Filter out points outside mesh
% Vector: point to nearest centreFaces -> centreFaces Normal
verticesRef = centreFaces;
pointsInside = verticesNew;
% Create a KD-tree searcher for verticesB
MdlB = KDTreeSearcher(verticesRef);
% Find the nearest neighbor in verticesB for each vertex in verticesA
K = 1;
[idxPoints, distPoints] = knnsearch(MdlB, pointsInside, 'K', K);
% Calculate the distances based on the nearest neighbor indices
minDistancesPoints = sqrt(sum((pointsInside - verticesRef(idxPoints,:)).^2, 2));

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
shrink.inside.inoutVerticesThicken = insideOutsidePoints; % Logical mask
insidePoints = pointsInside(insideOutsidePoints,:);
shrink.inside.verticesThicken = insidePoints;

disp(['thicken mesh: pelvis defect ',num2str(pelvisNum)])
end