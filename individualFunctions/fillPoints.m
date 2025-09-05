%% Fill mesh with points (for intersection)

% Input:    pelvisNum: Numeric identifier used only for logging
%           refData: Struct with triangulated reference surface (vertices, faces)
%           minDist: Spacing for the regular point grid (same units as refData)
%           box: Oriented bounding-box struct with fields
%                   - edgeVector: box edge vectors (local axes, RH rule)
%                   - cornerpoints: ordered corner points in world coords
%                   - tri: (optional) faces to render the box in a check-plot

% Output:   gridPoints: Struct with generated grid points:
%                   - pointsBox: all grid points inside box
%                   - pointsBoxNum: number of points in pointsBox
%                   - inside: subset of pointsBox that lie inside refData mesh
%                   - insideMask: logical mask over pointsBox for "inside"
%                   - insideMaskIdx: indices into pointsBox of inside points

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [gridPoints] = fillPoints(pelvisNum,refData,minDist,box)

% Bounding Box (not aligned to axis)
% Assuming: cornerPoints is a 8x3 matrix where each row is a corner point
% Bounding box's axis directions
dir1 = box.edgeVector(1,:)/norm(box.edgeVector(1,:)); % x % One box edge; right-hand-rule
dir2 = box.edgeVector(2,:)/norm(box.edgeVector(2,:)); % y
dir3 = box.edgeVector(3,:)/norm(box.edgeVector(3,:)); % z

% Create the rotation matrix local -> global (box)
rotMatrix = [dir1; dir2; dir3];

% Bounding box's dimensions
boxSize = [norm(box.edgeVector(1,:)), norm(box.edgeVector(2,:)), norm(box.edgeVector(3,:))];
% Determine how many points fit in each dimension based on minDist
% Desired minimum distance between points: minDist
numCells = ceil(boxSize ./ minDist);
% Generate points in bounding box's local coordinate system
numPoints = prod(numCells); % Total number of points
pointsLocal = zeros(numPoints, 3); % Preallocate the matrix for efficiency
for idx = 1:numPoints
    [i, j, k] = ind2sub(numCells, idx); % Convert linear index to 3D index
    x = (i-1)*minDist;
    y = (j-1)*minDist;
    z = (k-1)*minDist;
    pointsLocal(idx,:) = [x, y, z];
end

% Transform Points to global cosy
origin = box.cornerpoints(1,:);  % assuming this is the origin of the bounding box in world coordinates
pointsAll = pointsLocal * rotMatrix;
pointsAll = pointsAll + origin;
gridPoints.pointsBox = pointsAll;
gridPoints.pointsBoxNum = numPoints;

% Check

figure
plot3(pointsAll(:,1), pointsAll(:,2), ...
    pointsAll(:,3), '.','Color','b', 'MarkerSize', 5);
hold on
trisurf(box.tri,...
    box.cornerpoints(:,1),...
    box.cornerpoints(:,2),...
    box.cornerpoints(:,3),...
    'FaceColor','r','EdgeColor','r','FaceAlpha',0.25);
% Vertices (V) and faces (F) define the mesh
V = refData.vertices;
F = refData.faces;
% Check which points are inside the mesh
inside = inpolyhedron(F, V, pointsAll);
pointsInside = pointsAll(inside, :);
gridPoints.inside = pointsInside;
gridPoints.insideMask = inside; % base: points in box (pointsBox)
gridPoints.insideMaskIdx = find(inside);

disp(['pelvis filled with points: pelvis ', num2str(pelvisNum)]);

end