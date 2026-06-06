%% Acetabulum region: defect points inside sub-cuboids (grid points)

% Input:    acRegion: Struct with acetabulum cuboid region to be updated
%           isPointInsideDefect: Logical mask over all grid points / pointsBox
%                   - true for grid points inside defect
%                   - false for grid points outside defect
%           refCuboid: Struct with reference acetabulum sub-cuboids
%                   - c1.gridIdx: grid point indices in sub-cuboid c1
%                   - c2.gridIdx: grid point indices in sub-cuboid c2
%                   - c3.gridIdx: grid point indices in sub-cuboid c3
%                   - c4.gridIdx: grid point indices in sub-cuboid c4
%           voxelSize: Voxel spacing in x-, y-, and z-direction
%           pelvisNum: Numeric identifier used only for logging
%
% Output:   acRegion: Struct with updated acetabulum cuboid defect results:
%                   - cuboid.totalDefect.defectGridIdx: all defect grid point indices
%                   - cuboid.totalDefect.nGridPoints: total number of defect grid points
%                   - cuboid.totalDefect.gridVolumeLoss: total defect volume loss
%                   - cuboid.c1.defectGridIdx: defect grid point indices in c1
%                   - cuboid.c2.defectGridIdx: defect grid point indices in c2
%                   - cuboid.c3.defectGridIdx: defect grid point indices in c3
%                   - cuboid.c4.defectGridIdx: defect grid point indices in c4
%                   - cuboid.c1.nGridPoints: number of defect grid points in c1
%                   - cuboid.c2.nGridPoints: number of defect grid points in c2
%                   - cuboid.c3.nGridPoints: number of defect grid points in c3
%                   - cuboid.c4.nGridPoints: number of defect grid points in c4
%                   - cuboid.c1.gridLossPercentDefect: percent of total defect in c1
%                   - cuboid.c2.gridLossPercentDefect: percent of total defect in c2
%                   - cuboid.c3.gridLossPercentDefect: percent of total defect in c3
%                   - cuboid.c4.gridLossPercentDefect: percent of total defect in c4
%                   - cuboid.c1.gridLossPercentCuboid: percent loss relative to reference c1
%                   - cuboid.c2.gridLossPercentCuboid: percent loss relative to reference c2
%                   - cuboid.c3.gridLossPercentCuboid: percent loss relative to reference c3
%                   - cuboid.c4.gridLossPercentCuboid: percent loss relative to reference c4
%                   - cuboid.c1.gridVolumeLoss: volume loss in c1
%                   - cuboid.c2.gridVolumeLoss: volume loss in c2
%                   - cuboid.c3.gridVolumeLoss: volume loss in c3
%                   - cuboid.c4.gridVolumeLoss: volume loss in c4
%                   - cuboid.outside.defectGridIdx: defect grid points outside all sub-cuboids
%                   - cuboid.outside.nGridPoints: number of outside defect grid points
%                   - cuboid.outside.gridLossPercentDefect: percent of total defect outside sub-cuboids
%                   - cuboid.outside.gridVolumeLoss: outside defect volume loss

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function acRegion = acCuboidDefect(acRegion, isPointInsideDefect, refCuboid, voxelSize, pelvisNum)

% Sub-cuboids of reference pelvis
c1.gridIdx = refCuboid.c1.gridIdx;
c2.gridIdx = refCuboid.c2.gridIdx;
c3.gridIdx = refCuboid.c3.gridIdx;
c4.gridIdx = refCuboid.c4.gridIdx;
% Reference grid points per cuboid (total number)
nRefGridC1 = numel(c1.gridIdx);
nRefGridC2 = numel(c2.gridIdx);
nRefGridC3 = numel(c3.gridIdx);
nRefGridC4 = numel(c4.gridIdx);

% Voxel size [mm]
%voxelSize = [0.5 0.5 0.5];
voxelVolume = prod(voxelSize);   % = 0.125 mm^3

% global indices of defect grid points
idxDefectGrid = find(isPointInsideDefect);
% total number of defect grid points
nDefectGrid_total = numel(idxDefectGrid);

% defect grid points assigned to each cuboid
idxDefectGridC1 = intersect(idxDefectGrid, c1.gridIdx);
idxDefectGridC2 = intersect(idxDefectGrid, c2.gridIdx);
idxDefectGridC3 = intersect(idxDefectGrid, c3.gridIdx);
idxDefectGridC4 = intersect(idxDefectGrid, c4.gridIdx);
% defect grid points outside all cuboids
idxDefectGridInCuboids = [idxDefectGridC1; idxDefectGridC2; idxDefectGridC3; idxDefectGridC4];
idxDefectGridOutside = setdiff(idxDefectGrid, idxDefectGridInCuboids);
% number of defect grid points per cuboid
nDefectGridC1 = numel(idxDefectGridC1);
nDefectGridC2 = numel(idxDefectGridC2);
nDefectGridC3 = numel(idxDefectGridC3);
nDefectGridC4 = numel(idxDefectGridC4);
nDefectGridOutside = numel(idxDefectGridOutside);

% Proportion of total defect grid points
percentDefectGridC1 = 100 * nDefectGridC1 / nDefectGrid_total;
percentDefectGridC2 = 100 * nDefectGridC2 / nDefectGrid_total;
percentDefectGridC3 = 100 * nDefectGridC3 / nDefectGrid_total;
percentDefectGridC4 = 100 * nDefectGridC4 / nDefectGrid_total;
percentDefectGridOutside = 100 * nDefectGridOutside / nDefectGrid_total;
% Loss percent relative to reference pelvis grid points inside each cuboid
percentCuboidGridC1 = 100 * nDefectGridC1 / nRefGridC1;
percentCuboidGridC2 = 100 * nDefectGridC2 / nRefGridC2;
percentCuboidGridC3 = 100 * nDefectGridC3 / nRefGridC3;
percentCuboidGridC4 = 100 * nDefectGridC4 / nRefGridC4;

% Save (points)
acRegion.cuboid.totalDefect.defectGridIdx = idxDefectGrid; % total defect, not all sub-cuboids (defect can also be outside sub-cuboids)
acRegion.cuboid.totalDefect.nGridPoints = nDefectGrid_total; 
acRegion.cuboid.totalDefect.gridVolumeLoss = nDefectGrid_total * voxelVolume; 

acRegion.cuboid.c1.defectGridIdx = idxDefectGridC1; % grid points (idx) in sub-cuboid
acRegion.cuboid.c1.nGridPoints = nDefectGridC1; % number of grid points in sub-cuboid
acRegion.cuboid.c1.gridLossPercentDefect = percentDefectGridC1; % volume loss (percent) per defect
acRegion.cuboid.c1.gridLossPercentCuboid = percentCuboidGridC1; % volume loss (percent) per sub-cuboid
acRegion.cuboid.c1.gridVolumeLoss = nDefectGridC1 * voxelVolume; % volume loss (mm3)

acRegion.cuboid.c2.defectGridIdx = idxDefectGridC2;
acRegion.cuboid.c2.nGridPoints = nDefectGridC2;
acRegion.cuboid.c2.gridLossPercentDefect = percentDefectGridC2;
acRegion.cuboid.c2.gridLossPercentCuboid = percentCuboidGridC2;
acRegion.cuboid.c2.gridVolumeLoss = nDefectGridC2 * voxelVolume;

acRegion.cuboid.c3.defectGridIdx = idxDefectGridC3;
acRegion.cuboid.c3.nGridPoints = nDefectGridC3;
acRegion.cuboid.c3.gridLossPercentDefect = percentDefectGridC3;
acRegion.cuboid.c3.gridLossPercentCuboid = percentCuboidGridC3;
acRegion.cuboid.c3.gridVolumeLoss = nDefectGridC3 * voxelVolume;

acRegion.cuboid.c4.defectGridIdx = idxDefectGridC4;
acRegion.cuboid.c4.nGridPoints = nDefectGridC4;
acRegion.cuboid.c4.gridLossPercentDefect = percentDefectGridC4;
acRegion.cuboid.c4.gridLossPercentCuboid = percentCuboidGridC4;
acRegion.cuboid.c4.gridVolumeLoss = nDefectGridC4 * voxelVolume;

acRegion.cuboid.outside.defectGridIdx = idxDefectGridOutside;
acRegion.cuboid.outside.nGridPoints = nDefectGridOutside;
acRegion.cuboid.outside.gridLossPercentDefect = percentDefectGridOutside;
acRegion.cuboid.outside.gridVolumeLoss = nDefectGridOutside * voxelVolume;

disp(['defect grid points/vertices assigned to acetabulum sub-cuboids: pelvis defect ', num2str(pelvisNum)]);

end