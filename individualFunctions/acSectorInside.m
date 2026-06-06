%% Acetabulum region: points / vertices inside sectors

% Input:    acRegion: Struct with acetabulum sector region
%                   - sector.limits.centre: centre point of sector system [x y z]
%                   - sector.s1: superior annular sector
%                   - sector.s2: annular sector on negative y side
%                   - sector.s3: annular sector on positive y side
%                   - sector.s4: inner full cylinder
%           insidePoints: Points inside acetabulum region / mesh
%           insidePointsIdx: Indices of insidePoints relative to all grid points / pointsBox
%           vertices: Mesh vertices to be assigned to sectors
%
% Output:   acRegion: Struct with updated acetabulum sector region:
%                   - sector.s1.gridIdx: grid point indices in sector s1
%                   - sector.s2.gridIdx: grid point indices in sector s2
%                   - sector.s3.gridIdx: grid point indices in sector s3
%                   - sector.s4.gridIdx: grid point indices in sector s4
%                   - sector.s1.vertexIdx: vertex indices in sector s1
%                   - sector.s2.vertexIdx: vertex indices in sector s2
%                   - sector.s3.vertexIdx: vertex indices in sector s3
%                   - sector.s4.vertexIdx: vertex indices in sector s4

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function acRegion = acSectorInside(acRegion, insidePoints, insidePointsIdx, vertices)

% helper
centre = acRegion.sector.limits.centre;
x0 = centre(1);
y0 = centre(2);
z0 = centre(3);
% grid points
x = insidePoints(:,1);
y = insidePoints(:,2);
z = insidePoints(:,3);

rYZ = sqrt((y - y0).^2 + (z - z0).^2);
thetaDeg = atan2d(y - y0, z - z0);   % 0° at +z

% Sector 1
s1 = acRegion.sector.s1;
isInS1 = x >= s1.xmin & x <= s1.xmax & ...
    rYZ >= s1.rmin & rYZ <= s1.rmax & ...
    thetaDeg >= s1.thetaMinDeg & thetaDeg <= s1.thetaMaxDeg;
% Sector 2
s2 = acRegion.sector.s2;
isInS2 = x >= s2.xmin & x <= s2.xmax & ...
    rYZ >= s2.rmin & rYZ <= s2.rmax & ...
    thetaDeg >= s2.thetaMinDeg & thetaDeg <= s2.thetaMaxDeg;
% Sector 3
s3 = acRegion.sector.s3;
isInS3 = x >= s3.xmin & x <= s3.xmax & ...
    rYZ >= s3.rmin & rYZ <= s3.rmax & ...
    thetaDeg >= s3.thetaMinDeg & thetaDeg <= s3.thetaMaxDeg;
% Sector 4
s4 = acRegion.sector.s4;
isInS4 = x >= s4.xmin & x <= s4.xmax & ...
    rYZ >= s4.rmin & rYZ <= s4.rmax;

% global grid indices
acRegion.sector.s1.gridIdx = insidePointsIdx(isInS1);
acRegion.sector.s2.gridIdx = insidePointsIdx(isInS2);
acRegion.sector.s3.gridIdx = insidePointsIdx(isInS3);
acRegion.sector.s4.gridIdx = insidePointsIdx(isInS4);

% vertices
xv = vertices(:,1);
yv = vertices(:,2);
zv = vertices(:,3);
rYZv = sqrt((yv - y0).^2 + (zv - z0).^2);
thetaDegV = atan2d(yv - y0, zv - z0);

% Sector 1
isInS1v = xv >= s1.xmin & xv <= s1.xmax & ...
    rYZv >= s1.rmin & rYZv <= s1.rmax & ...
    thetaDegV >= s1.thetaMinDeg & thetaDegV <= s1.thetaMaxDeg;
% Sector 2
isInS2v = xv >= s2.xmin & xv <= s2.xmax & ...
    rYZv >= s2.rmin & rYZv <= s2.rmax & ...
    thetaDegV >= s2.thetaMinDeg & thetaDegV <= s2.thetaMaxDeg;
% Sector 3
isInS3v = xv >= s3.xmin & xv <= s3.xmax & ...
    rYZv >= s3.rmin & rYZv <= s3.rmax & ...
    thetaDegV >= s3.thetaMinDeg & thetaDegV <= s3.thetaMaxDeg;
% Sector 4
isInS4v = xv >= s4.xmin & xv <= s4.xmax & ...
    rYZv >= s4.rmin & rYZv <= s4.rmax;

% global vertices indices
acRegion.sector.s1.vertexIdx = find(isInS1v);
acRegion.sector.s2.vertexIdx = find(isInS2v);
acRegion.sector.s3.vertexIdx = find(isInS3v);
acRegion.sector.s4.vertexIdx = find(isInS4v);

disp('grid points / vertices assigned to acetabulum sectors')

end