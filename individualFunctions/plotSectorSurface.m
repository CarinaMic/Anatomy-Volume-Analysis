%% Function to plot sectors

% Input:    sector: Struct with sector geometry
%                   - centre: centre point of sector system [x y z]
%                   - xmin, xmax: minimum and maximum x-coordinate / x-depth
%                   - rmin, rmax: minimum and maximum radius
%                   - thetaMinDeg, thetaMaxDeg: angular limits in degrees
%           faceColor: Color used for sector surface faces and edges
%           faceAlpha: Transparency of sector surface faces
%           nTheta: Number of angular sampling points for sector surface
%                   Optional, default = 80
%
% Output:   h: Graphics handle of outer cylindrical sector surface
%                   - can be used for legend entries

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function h = plotSectorSurface(sector, faceColor, faceAlpha, nTheta)

if nargin < 5
    nTheta = 80;
end
c = sector.centre;
y0 = c(2);
z0 = c(3);
xmin = sector.xmin;
xmax = sector.xmax;
rmin = sector.rmin;
rmax = sector.rmax;
th1 = deg2rad(sector.thetaMinDeg);
th2 = deg2rad(sector.thetaMaxDeg);
theta = linspace(th1, th2, nTheta);
% Outer arc
Yout = y0 + rmax * sin(theta);
Zout = z0 + rmax * cos(theta);
% Inner arc
Yin = y0 + rmin * sin(theta);
Zin = z0 + rmin * cos(theta);
% Outer cylindrical surface -> use this handle for legend
Xo = [xmin * ones(1,nTheta); xmax * ones(1,nTheta)];
Yo = [Yout; Yout];
Zo = [Zout; Zout];
h = surf(Xo, Yo, Zo, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha, ...
    'EdgeColor', faceColor, 'LineWidth', 1.2); %
% Inner cylindrical surface
if rmin > 0
    Xi = [xmin * ones(1,nTheta); xmax * ones(1,nTheta)];
    Yi = [Yin; Yin];
    Zi = [Zin; Zin];
    surf(Xi, Yi, Zi, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha, ...
        'EdgeColor', faceColor, 'LineWidth', 1.2);
end
% Radial side surfaces
for k = [1 nTheta]
    Xr = [xmin xmax; xmin xmax];
    Yr = [Yin(k) Yin(k); Yout(k) Yout(k)];
    Zr = [Zin(k) Zin(k); Zout(k) Zout(k)];
    surf(Xr, Yr, Zr, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha, ...
        'EdgeColor', faceColor, 'LineWidth', 1.2);
end
% Front / back faces
if rmin > 0
    Yf = [Yout fliplr(Yin)];
    Zf = [Zout fliplr(Zin)];
else
    Yf = [Yout y0];
    Zf = [Zout z0];
end
Xf = xmin * ones(size(Yf));
patch(Xf, Yf, Zf, faceColor, 'FaceAlpha', faceAlpha, ...
    'EdgeColor', faceColor, 'LineWidth', 1.2);
Xb = xmax * ones(size(Yf));
patch(Xb, Yf, Zf, faceColor, 'FaceAlpha', faceAlpha, ...
    'EdgeColor', faceColor, 'LineWidth', 1.2);
end
