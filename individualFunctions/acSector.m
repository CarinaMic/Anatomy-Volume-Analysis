%% Acetabulum region: sectors for acetabulum region

% Input:    limits: Struct with acetabulum sector limits and settings
%                   - centre: centre point of sector system [x y z]
%                   - xmin, xmax: minimum and maximum x-coordinate / x-depth
%                   - rInner: inner radius for annular sectors
%                   - rOuter23: outer radius for sectors s2 and s3
%                   - rOuter1: outer radius for sector s1
%                   - s1ThetaMinDeg, s1ThetaMaxDeg: angular limits of sector s1
%                   - s2ThetaMinDeg, s2ThetaMaxDeg: angular limits of sector s2
%                   - s3ThetaMinDeg, s3ThetaMaxDeg: angular limits of sector s3
%                   - s4ThetaMinDeg, s4ThetaMaxDeg: angular limits of sector s4
%
% Output:   acRegion: Struct with acetabulum sector region:
%                   - sector.s1: superior annular sector
%                   - sector.s2: annular sector on negative y side
%                   - sector.s3: annular sector on positive y side
%                   - sector.s4: inner full cylinder
%               Each sector contains:
%                   - type: sector type, either 'cylinder' or 'annularSector'
%                   - centre: centre point of sector system
%                   - xmin, xmax: minimum and maximum x-coordinate / x-depth
%                   - rmin, rmax: minimum and maximum radius
%                   - thetaMinDeg, thetaMaxDeg: angular limits in degrees

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function acRegion = acSector(limits)

% Centre
centre = limits.centre;
x0 = centre(1);
y0 = centre(2);
z0 = centre(3);
% x-depth
xmin = limits.xmin;
xmax = limits.xmax;
% radii
rInner   = limits.rInner;
rOuter23 = limits.rOuter23;
rOuter1  = limits.rOuter1;

% Save common sector settings
%obj.acRegion.sector.limits = limits;
%obj.acRegion.sector.centre = centre;

% Sector 4: inner full cylinder
acRegion.sector.s4.type = 'cylinder';
acRegion.sector.s4.centre = centre;
acRegion.sector.s4.xmin = xmin;
acRegion.sector.s4.xmax = xmax;
acRegion.sector.s4.rmin = 0;
acRegion.sector.s4.rmax = rInner;
acRegion.sector.s4.thetaMinDeg = limits.s4ThetaMinDeg;
acRegion.sector.s4.thetaMaxDeg = limits.s4ThetaMaxDeg;

% Sector 2: annular sector on negative y side
acRegion.sector.s2.type = 'annularSector';
acRegion.sector.s2.centre = centre;
acRegion.sector.s2.xmin = xmin;
acRegion.sector.s2.xmax = xmax;
acRegion.sector.s2.rmin = rInner;
acRegion.sector.s2.rmax = rOuter23;
acRegion.sector.s2.thetaMinDeg = limits.s2ThetaMinDeg;
acRegion.sector.s2.thetaMaxDeg = limits.s2ThetaMaxDeg;

% Sector 3: annular sector on positive y side
acRegion.sector.s3.type = 'annularSector';
acRegion.sector.s3.centre = centre;
acRegion.sector.s3.xmin = xmin;
acRegion.sector.s3.xmax = xmax;
acRegion.sector.s3.rmin = rInner;
acRegion.sector.s3.rmax = rOuter23;
acRegion.sector.s3.thetaMinDeg = limits.s3ThetaMinDeg;
acRegion.sector.s3.thetaMaxDeg = limits.s3ThetaMaxDeg;

% Sector 1: superior annular sector
acRegion.sector.s1.type = 'annularSector';
acRegion.sector.s1.centre = centre;
acRegion.sector.s1.xmin = xmin;
acRegion.sector.s1.xmax = xmax;
acRegion.sector.s1.rmin = rInner;
acRegion.sector.s1.rmax = rOuter1;
acRegion.sector.s1.thetaMinDeg = limits.s1ThetaMinDeg;
acRegion.sector.s1.thetaMaxDeg = limits.s1ThetaMaxDeg;

disp('acetabulum sector region created')

end