%% Acetabulum region: hybrid = cuboids with central sector

% Input:    acRegion: Struct with acetabulum region to be updated
%           cuboidRef: Struct with reference acetabulum cuboids
%                   - all.limits: limits of full cuboid region
%                   - c1, c2, c3, c4: sub-cuboids with vertices, faces, limits,
%                     gridIdx and vertexIdx
%           sectorRef: Struct with reference acetabulum sectors
%                   - limits.centre: centre point of sector system [x y z]
%                   - s4: central sector with rmin, rmax, xmin, xmax,
%                     gridIdx and vertexIdx
%           insidePoints: Points inside acetabulum region / mesh
%           insidePointsIdx: Indices of insidePoints relative to all grid points / pointsBox
%           vertices: Mesh vertices to be assigned to hybrid regions
%
% Output:   acRegion: Struct with updated acetabulum hybrid regions:
%                   - hybrid.limits: limits of full hybrid region
%                   - hybrid.mc1: modified cuboid c1 without central sector mc5
%                   - hybrid.mc2: modified cuboid c2 without central sector mc5
%                   - hybrid.mc3: modified cuboid c3 without central sector mc5
%                   - hybrid.mc4: modified cuboid c4 without central sector mc5
%                   - hybrid.mc5: central sector region
%               Each hybrid region contains:
%                   - gridIdx: grid point indices in hybrid sub-region
%                   - vertexIdx: vertex indices in hybrid sub-region
%               Hybrid regions mc1-mc4 additionally contain:
%                   - vertices: cuboid corner points for display
%                   - faces: cuboid faces for plotting/rendering
%                   - limits: coordinate limits of the cuboid
%               Hybrid region mc5 additionally contains:
%                   - type: sector type
%                   - centre: centre point of sector system
%                   - xmin, xmax: adapted x-depth
%                   - rmin, rmax: minimum and maximum radius
%                   - thetaMinDeg, thetaMaxDeg: angular limits in degrees

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich



function acRegion = acHybridInside(acRegion, cuboidRef, sectorRef, insidePoints, insidePointsIdx, vertices)

% Copy geometry/limits from cuboid
acRegion.hybrid.limits = cuboidRef.all.limits;

% Keep cuboid geometry for display
acRegion.hybrid.mc1.vertices = cuboidRef.c1.vertices;
acRegion.hybrid.mc1.faces    = cuboidRef.c1.faces;
acRegion.hybrid.mc1.limits   = cuboidRef.c1.limits;

acRegion.hybrid.mc2.vertices = cuboidRef.c2.vertices;
acRegion.hybrid.mc2.faces    = cuboidRef.c2.faces;
acRegion.hybrid.mc2.limits   = cuboidRef.c2.limits;

acRegion.hybrid.mc3.vertices = cuboidRef.c3.vertices;
acRegion.hybrid.mc3.faces    = cuboidRef.c3.faces;
acRegion.hybrid.mc3.limits   = cuboidRef.c3.limits;

acRegion.hybrid.mc4.vertices = cuboidRef.c4.vertices;
acRegion.hybrid.mc4.faces    = cuboidRef.c4.faces;
acRegion.hybrid.mc4.limits   = cuboidRef.c4.limits;

% Central sector mc5 geometry
acRegion.hybrid.mc5 = sectorRef.s4;
% Adapt central sector depth to cuboid depth, e.g. 2.2 * acentreR
acRegion.hybrid.mc5.xmin = cuboidRef.c1.limits.xmin;
acRegion.hybrid.mc5.xmax = cuboidRef.c1.limits.xmax;

% s4 indices to remove from cuboids
%mc5GridIdx   = sectorRef.s4.gridIdx;
%mc5VertexIdx = sectorRef.s4.vertexIdx;
% Recalculate mc5 indices with adapted x-depth
centre = sectorRef.limits.centre;
x0 = centre(1);
y0 = centre(2);
z0 = centre(3);
mc5 = acRegion.hybrid.mc5;
% Grid points
x = insidePoints(:,1);
y = insidePoints(:,2);
z = insidePoints(:,3);
rYZ = sqrt((y - y0).^2 + (z - z0).^2);
isInMc5 = x >= mc5.xmin & x <= mc5.xmax & ...
    rYZ >= mc5.rmin & rYZ <= mc5.rmax;
mc5GridIdx = insidePointsIdx(isInMc5);
% Vertices
xv = vertices(:,1);
yv = vertices(:,2);
zv = vertices(:,3);
rYZv = sqrt((yv - y0).^2 + (zv - z0).^2);
isInMc5v = xv >= mc5.xmin & xv <= mc5.xmax & ...
    rYZv >= mc5.rmin & rYZv <= mc5.rmax;
mc5VertexIdx = find(isInMc5v);

% Grid indices: cuboid minus mc5/s4
acRegion.hybrid.mc1.gridIdx = setdiff(cuboidRef.c1.gridIdx, mc5GridIdx);
acRegion.hybrid.mc2.gridIdx = setdiff(cuboidRef.c2.gridIdx, mc5GridIdx);
acRegion.hybrid.mc3.gridIdx = setdiff(cuboidRef.c3.gridIdx, mc5GridIdx);
acRegion.hybrid.mc4.gridIdx = setdiff(cuboidRef.c4.gridIdx, mc5GridIdx);

% Vertex indices: cuboid minus mc5/s4
acRegion.hybrid.mc1.vertexIdx = setdiff(cuboidRef.c1.vertexIdx, mc5VertexIdx);
acRegion.hybrid.mc2.vertexIdx = setdiff(cuboidRef.c2.vertexIdx, mc5VertexIdx);
acRegion.hybrid.mc3.vertexIdx = setdiff(cuboidRef.c3.vertexIdx, mc5VertexIdx);
acRegion.hybrid.mc4.vertexIdx = setdiff(cuboidRef.c4.vertexIdx, mc5VertexIdx);

% Central sector remains as own region
acRegion.hybrid.mc5.gridIdx   = mc5GridIdx;
acRegion.hybrid.mc5.vertexIdx = mc5VertexIdx;

disp('grid points / vertices assigned to acetabulum hybrid regions')
end