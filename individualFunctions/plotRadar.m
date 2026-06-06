%% Display radar plot per defect

% Input:    values: Numeric values to be displayed in the radar plot
%                   - one value per radar axis
%           plotTitle: Title of the radar plot
%           labels: Cell array/string array with axis labels
%                   - number of labels must match number of values
%
% Output:   None
%                   - creates radar plot in current axes

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function plotRadar(values, plotTitle, labels)
plotColor = [0, 101/255, 189/255]; % TUMcolors.blue300

values = values(:)';   % row vector
labels = labels(:)';   % row vector
nAxes = numel(values);

if numel(labels) ~= nAxes
    error('plotRadar:SizeMismatch', ...
        'Number of labels (%d) must match number of values (%d).', ...
        numel(labels), nAxes);
end

rMax = 100;
rTicks = [0 25 50 75 100];

% Dynamic axis positions:
% first axis at top, then clockwise
thetaDeg = 90 - (0:nAxes-1) * 360/nAxes;
theta = deg2rad(thetaDeg);

% Close polygon
thetaClosed = [theta theta(1)];
valuesClosed = [values values(1)];

% Cartesian coordinates
xPoly = valuesClosed .* cos(thetaClosed);
yPoly = valuesClosed .* sin(thetaClosed);

hold on
axis equal
axis off

% Grid polygons
for r = rTicks(2:end)
    xg = r * cos(thetaClosed);
    yg = r * sin(thetaClosed);
    plot(xg, yg, '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 1);
end

% Axes
for k = 1:nAxes
    plot([0 rMax*cos(theta(k))], ...
        [0 rMax*sin(theta(k))], '-', ...
        'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1);
end

% Filled polygon
patch(xPoly, yPoly, plotColor, ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', plotColor, ...
    'LineWidth', 2);

% Markers
plot(xPoly, yPoly, 'o', ...
    'Color', plotColor, ...
    'MarkerFaceColor', plotColor, ...
    'MarkerSize', 5);

% Labels
labelRadius = rMax * 1.12;
for k = 1:nAxes
    xl = labelRadius * cos(theta(k));
    yl = labelRadius * sin(theta(k));

    text(xl, yl, labels{k}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
end

% Tick labels near top
for r = rTicks
    text(-8, r, num2str(r), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 10, ...
        'Color', [0.1 0.1 0.1]);
end

xlim([-1.25*rMax, 1.25*rMax]);
ylim([-1.15*rMax, 1.15*rMax]);

title(plotTitle, 'FontWeight', 'bold');
hold off
end