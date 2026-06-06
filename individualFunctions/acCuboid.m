%% Acetabulum region: cuboids for acetabulum region

% Input:    limits: Struct with limits and midpoints of acetabulum region
%                   - xmin, xmax: minimum and maximum x-coordinate
%                   - ymin, ymax: minimum and maximum y-coordinate
%                   - zmin, zmax: minimum and maximum z-coordinate
%                   - xmid: midpoint x-coordinate
%                   - ymid: midpoint y-coordinate
%                   - zmid: midpoint z-coordinate

% Output:   acRegion: Struct with acetabulum cuboid region:
%                   - cuboid.all: large cuboid enclosing the full region
%                   - cuboid.c1: sub-cuboid y+ / z+
%                   - cuboid.c2: sub-cuboid y+ / z-
%                   - cuboid.c3: sub-cuboid y- / z-
%                   - cuboid.c4: sub-cuboid y- / z+
%               Each cuboid contains:
%                   - limits: coordinate limits of the cuboid
%                   - vertices: cuboid corner points
%                   - faces: cuboid faces for plotting/rendering

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function acRegion = acCuboid(limits)

% Large cuboid limits
xmin = limits.xmin;
xmax = limits.xmax;
ymin = limits.ymin;
ymax = limits.ymax;
zmin = limits.zmin;
zmax = limits.zmax;
xmid = limits.xmid;
ymid = limits.ymid;
zmid = limits.zmid;

% Large cuboid
acRegion.cuboid.all.limits = struct( ...
    'xmin', xmin, 'xmax', xmax, ...
    'ymin', ymin, 'ymax', ymax, ...
    'zmin', zmin, 'zmax', zmax);
acRegion.cuboid.all.vertices = [
    xmin ymin zmin;
    xmax ymin zmin;
    xmax ymax zmin;
    xmin ymax zmin;
    xmin ymin zmax;
    xmax ymin zmax;
    xmax ymax zmax;
    xmin ymax zmax];
acRegion.cuboid.all.faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

% q1: y+ / z+
acRegion.cuboid.c1.limits = struct( ...
    'xmin', xmin, 'xmax', xmax, ...
    'ymin', ymid, 'ymax', ymax, ...
    'zmin', zmid, 'zmax', zmax);
acRegion.cuboid.c1.vertices = [
    xmin ymid zmid;
    xmax ymid zmid;
    xmax ymax zmid;
    xmin ymax zmid;
    xmin ymid zmax;
    xmax ymid zmax;
    xmax ymax zmax;
    xmin ymax zmax];
acRegion.cuboid.c1.faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

% q2: y+ / z-
acRegion.cuboid.c2.limits = struct( ...
    'xmin', xmin, 'xmax', xmax, ...
    'ymin', ymid, 'ymax', ymax, ...
    'zmin', zmin, 'zmax', zmid);
acRegion.cuboid.c2.vertices = [
    xmin ymid zmin;
    xmax ymid zmin;
    xmax ymax zmin;
    xmin ymax zmin;
    xmin ymid zmid;
    xmax ymid zmid;
    xmax ymax zmid;
    xmin ymax zmid];
acRegion.cuboid.c2.faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

% q3: y- / z-
acRegion.cuboid.c3.limits = struct( ...
    'xmin', xmin, 'xmax', xmax, ...
    'ymin', ymin, 'ymax', ymid, ...
    'zmin', zmin, 'zmax', zmid);
acRegion.cuboid.c3.vertices = [
    xmin ymin zmin;
    xmax ymin zmin;
    xmax ymid zmin;
    xmin ymid zmin;
    xmin ymin zmid;
    xmax ymin zmid;
    xmax ymid zmid;
    xmin ymid zmid];
acRegion.cuboid.c3.faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

% q4: y- / z+
acRegion.cuboid.c4.limits = struct( ...
    'xmin', xmin, 'xmax', xmax, ...
    'ymin', ymin, 'ymax', ymid, ...
    'zmin', zmid, 'zmax', zmax);
acRegion.cuboid.c4.vertices = [
    xmin ymin zmid;
    xmax ymin zmid;
    xmax ymid zmid;
    xmin ymid zmid;
    xmin ymin zmax;
    xmax ymin zmax;
    xmax ymid zmax;
    xmin ymid zmax];
acRegion.cuboid.c4.faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

disp('acetabulum cuboid region created')

end