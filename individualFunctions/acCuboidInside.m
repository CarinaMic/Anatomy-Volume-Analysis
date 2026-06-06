%% Acetabulum region: points/vertices inside sub-cuboids

% Input:    acRegion: Struct with acetabulum cuboid region
%                   - cuboid.c1.limits: limits of sub-cuboid y+ / z+
%                   - cuboid.c2.limits: limits of sub-cuboid y+ / z-
%                   - cuboid.c3.limits: limits of sub-cuboid y- / z-
%                   - cuboid.c4.limits: limits of sub-cuboid y- / z+
%           insidePoints: Points inside acetabulum region / mesh
%           insidePointsMask: Indices of insidePoints relative to all grid points / pointsBox
%           vertices: Mesh vertices to be assigned to sub-cuboids

% Output:   acRegion: Struct with updated acetabulum cuboid region:
%                   - cuboid.c1.gridIdx: grid point indices in sub-cuboid c1
%                   - cuboid.c2.gridIdx: grid point indices in sub-cuboid c2
%                   - cuboid.c3.gridIdx: grid point indices in sub-cuboid c3
%                   - cuboid.c4.gridIdx: grid point indices in sub-cuboid c4
%                   - cuboid.c1.vertexIdx: vertex indices in sub-cuboid c1
%                   - cuboid.c2.vertexIdx: vertex indices in sub-cuboid c2
%                   - cuboid.c3.vertexIdx: vertex indices in sub-cuboid c3
%                   - cuboid.c4.vertexIdx: vertex indices in sub-cuboid c4

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function acRegion = acCuboidInside(acRegion, insidePoints, insidePointsMask, vertices)

% Cuboid limits
c1 = acRegion.cuboid.c1.limits;
c2 = acRegion.cuboid.c2.limits;
c3 = acRegion.cuboid.c3.limits;
c4 = acRegion.cuboid.c4.limits;

% Grid points inside sub-cuboids (clear assignment of limit values)
% c1: y+ / z+
isInC1_points = insidePoints(:,1) >= c1.xmin & insidePoints(:,1) <= c1.xmax & ...
    insidePoints(:,2) >= c1.ymin & insidePoints(:,2) <= c1.ymax & ...
    insidePoints(:,3) >= c1.zmin & insidePoints(:,3) <= c1.zmax;
% c2: y+ / z-
isInC2_points = insidePoints(:,1) >= c2.xmin & insidePoints(:,1) <= c2.xmax & ...
    insidePoints(:,2) >= c2.ymin & insidePoints(:,2) <= c2.ymax & ...
    insidePoints(:,3) >= c2.zmin & insidePoints(:,3) <  c2.zmax;
% c3: y- / z-
isInC3_points = insidePoints(:,1) >= c3.xmin & insidePoints(:,1) <= c3.xmax & ...
    insidePoints(:,2) >= c3.ymin & insidePoints(:,2) <  c3.ymax & ...
    insidePoints(:,3) >= c3.zmin & insidePoints(:,3) <  c3.zmax;
% c4: y- / z+
isInC4_points = insidePoints(:,1) >= c4.xmin & insidePoints(:,1) <= c4.xmax & ...
    insidePoints(:,2) >= c4.ymin & insidePoints(:,2) <  c4.ymax & ...
    insidePoints(:,3) >= c4.zmin & insidePoints(:,3) <= c4.zmax;

% Global indices relative to all grid points / pointsBox
idxC1_grid = insidePointsMask(isInC1_points);
idxC2_grid = insidePointsMask(isInC2_points);
idxC3_grid = insidePointsMask(isInC3_points);
idxC4_grid = insidePointsMask(isInC4_points);

% Vertices inside sub-cuboids
% c1: y+ / z+
isInC1_vertices = vertices(:,1) >= c1.xmin & vertices(:,1) <= c1.xmax & ...
    vertices(:,2) >= c1.ymin & vertices(:,2) <= c1.ymax & ...
    vertices(:,3) >= c1.zmin & vertices(:,3) <= c1.zmax;
% c2: y+ / z-
isInC2_vertices = vertices(:,1) >= c2.xmin & vertices(:,1) <= c2.xmax & ...
    vertices(:,2) >= c2.ymin & vertices(:,2) <= c2.ymax & ...
    vertices(:,3) >= c2.zmin & vertices(:,3) <  c2.zmax;
% c3: y- / z-
isInC3_vertices = vertices(:,1) >= c3.xmin & vertices(:,1) <= c3.xmax & ...
    vertices(:,2) >= c3.ymin & vertices(:,2) <  c3.ymax & ...
    vertices(:,3) >= c3.zmin & vertices(:,3) <  c3.zmax;
% c4: y- / z+
isInC4_vertices = vertices(:,1) >= c4.xmin & vertices(:,1) <= c4.xmax & ...
    vertices(:,2) >= c4.ymin & vertices(:,2) <  c4.ymax & ...
    vertices(:,3) >= c4.zmin & vertices(:,3) <= c4.zmax;

% Vertex indices relative to vertices array
idxC1_vertex = find(isInC1_vertices);
idxC2_vertex = find(isInC2_vertices);
idxC3_vertex = find(isInC3_vertices);
idxC4_vertex = find(isInC4_vertices);

% Save results
% Grid points / pointsBox
acRegion.cuboid.c1.gridIdx = idxC1_grid;
acRegion.cuboid.c2.gridIdx = idxC2_grid;
acRegion.cuboid.c3.gridIdx = idxC3_grid;
acRegion.cuboid.c4.gridIdx = idxC4_grid;
% Vertices
acRegion.cuboid.c1.vertexIdx = idxC1_vertex;
acRegion.cuboid.c2.vertexIdx = idxC2_vertex;
acRegion.cuboid.c3.vertexIdx = idxC3_vertex;
acRegion.cuboid.c4.vertexIdx = idxC4_vertex;

disp('grid points / vertices assigned to acetabulum sub-cuboids')

end