%% Acetabulum region: defect points inside sectors

% Input:    acRegion: Struct with acetabulum sector region to be updated
%           isPointInsideDefect: Logical mask over all grid points / pointsBox
%                   - true for grid points inside defect
%                   - false for grid points outside defect
%           refSector: Struct with reference acetabulum sectors
%                   - s1.gridIdx: grid point indices in sector s1
%                   - s2.gridIdx: grid point indices in sector s2
%                   - s3.gridIdx: grid point indices in sector s3
%                   - s4.gridIdx: grid point indices in sector s4
%           voxelSize: Voxel spacing in x-, y-, and z-direction
%           pelvisNum: Numeric identifier used only for logging
%
% Output:   acRegion: Struct with updated acetabulum sector defect results:
%                   - sector.totalDefect.defectGridIdx: all defect grid point indices
%                   - sector.totalDefect.nGridPoints: total number of defect grid points
%                   - sector.totalDefect.gridVolumeLoss: total defect volume loss
%                   - sector.s1.defectGridIdx: defect grid point indices in sector s1
%                   - sector.s2.defectGridIdx: defect grid point indices in sector s2
%                   - sector.s3.defectGridIdx: defect grid point indices in sector s3
%                   - sector.s4.defectGridIdx: defect grid point indices in sector s4
%                   - sector.s1.nGridPoints: number of defect grid points in sector s1
%                   - sector.s2.nGridPoints: number of defect grid points in sector s2
%                   - sector.s3.nGridPoints: number of defect grid points in sector s3
%                   - sector.s4.nGridPoints: number of defect grid points in sector s4
%                   - sector.s1.gridLossPercentDefect: percent of total defect in sector s1
%                   - sector.s2.gridLossPercentDefect: percent of total defect in sector s2
%                   - sector.s3.gridLossPercentDefect: percent of total defect in sector s3
%                   - sector.s4.gridLossPercentDefect: percent of total defect in sector s4
%                   - sector.s1.gridLossPercentSector: percent loss relative to reference sector s1
%                   - sector.s2.gridLossPercentSector: percent loss relative to reference sector s2
%                   - sector.s3.gridLossPercentSector: percent loss relative to reference sector s3
%                   - sector.s4.gridLossPercentSector: percent loss relative to reference sector s4
%                   - sector.s1.gridVolumeLoss: volume loss in sector s1
%                   - sector.s2.gridVolumeLoss: volume loss in sector s2
%                   - sector.s3.gridVolumeLoss: volume loss in sector s3
%                   - sector.s4.gridVolumeLoss: volume loss in sector s4
%                   - sector.outside.defectGridIdx: defect grid points outside all sectors
%                   - sector.outside.nGridPoints: number of outside defect grid points
%                   - sector.outside.gridLossPercentDefect: percent of total defect outside sectors
%                   - sector.outside.gridVolumeLoss: outside defect volume loss

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function acRegion = acSectorDefect(acRegion, isPointInsideDefect, refSector, voxelSize, pelvisNum)

% Reference sectors of pelvis
s1.gridIdx = refSector.s1.gridIdx;
s2.gridIdx = refSector.s2.gridIdx;
s3.gridIdx = refSector.s3.gridIdx;
s4.gridIdx = refSector.s4.gridIdx;

% Number of reference grid points per sector
nRefGridS1 = numel(s1.gridIdx);
nRefGridS2 = numel(s2.gridIdx);
nRefGridS3 = numel(s3.gridIdx);
nRefGridS4 = numel(s4.gridIdx);

% Voxel size [mm]
voxelVolume = prod(voxelSize);

% Global indices of defect grid points
idxDefectGrid = find(isPointInsideDefect);

% Total number of defect grid points
nDefectGrid_total = numel(idxDefectGrid);

% Defect grid points assigned to each sector
idxDefectGridS1 = intersect(idxDefectGrid, s1.gridIdx);
idxDefectGridS2 = intersect(idxDefectGrid, s2.gridIdx);
idxDefectGridS3 = intersect(idxDefectGrid, s3.gridIdx);
idxDefectGridS4 = intersect(idxDefectGrid, s4.gridIdx);

% Defect grid points outside all sectors
idxDefectGridInSectors = [idxDefectGridS1; idxDefectGridS2; idxDefectGridS3; idxDefectGridS4];
idxDefectGridOutside = setdiff(idxDefectGrid, idxDefectGridInSectors);

% Number of defect grid points per sector
nDefectGridS1 = numel(idxDefectGridS1);
nDefectGridS2 = numel(idxDefectGridS2);
nDefectGridS3 = numel(idxDefectGridS3);
nDefectGridS4 = numel(idxDefectGridS4);
nDefectGridOutside = numel(idxDefectGridOutside);

% Proportion of total defect grid points
percentDefectGridS1 = 100 * nDefectGridS1 / nDefectGrid_total;
percentDefectGridS2 = 100 * nDefectGridS2 / nDefectGrid_total;
percentDefectGridS3 = 100 * nDefectGridS3 / nDefectGrid_total;
percentDefectGridS4 = 100 * nDefectGridS4 / nDefectGrid_total;
percentDefectGridOutside = 100 * nDefectGridOutside / nDefectGrid_total;

% Loss percent relative to reference pelvis grid points inside each sector
percentSectorGridS1 = 100 * nDefectGridS1 / nRefGridS1;
percentSectorGridS2 = 100 * nDefectGridS2 / nRefGridS2;
percentSectorGridS3 = 100 * nDefectGridS3 / nRefGridS3;
percentSectorGridS4 = 100 * nDefectGridS4 / nRefGridS4;

% Save
acRegion.sector.totalDefect.defectGridIdx = idxDefectGrid; % total defect, not all sub-cuboids (defect can also be outside sectors)
acRegion.sector.totalDefect.nGridPoints = nDefectGrid_total;
acRegion.sector.totalDefect.gridVolumeLoss = nDefectGrid_total * voxelVolume;

acRegion.sector.s1.defectGridIdx = idxDefectGridS1; % grid points (idx) in sector
acRegion.sector.s1.nGridPoints = nDefectGridS1; % number of grid points in sector
acRegion.sector.s1.gridLossPercentDefect = percentDefectGridS1; % volume loss (percent) per defect
acRegion.sector.s1.gridLossPercentSector = percentSectorGridS1; % volume loss (percent) per sector
acRegion.sector.s1.gridVolumeLoss = nDefectGridS1 * voxelVolume; % volume loss (mm3)

acRegion.sector.s2.defectGridIdx = idxDefectGridS2;
acRegion.sector.s2.nGridPoints = nDefectGridS2;
acRegion.sector.s2.gridLossPercentDefect = percentDefectGridS2;
acRegion.sector.s2.gridLossPercentSector = percentSectorGridS2;
acRegion.sector.s2.gridVolumeLoss = nDefectGridS2 * voxelVolume;

acRegion.sector.s3.defectGridIdx = idxDefectGridS3;
acRegion.sector.s3.nGridPoints = nDefectGridS3;
acRegion.sector.s3.gridLossPercentDefect = percentDefectGridS3;
acRegion.sector.s3.gridLossPercentSector = percentSectorGridS3;
acRegion.sector.s3.gridVolumeLoss = nDefectGridS3 * voxelVolume;

acRegion.sector.s4.defectGridIdx = idxDefectGridS4;
acRegion.sector.s4.nGridPoints = nDefectGridS4;
acRegion.sector.s4.gridLossPercentDefect = percentDefectGridS4;
acRegion.sector.s4.gridLossPercentSector = percentSectorGridS4;
acRegion.sector.s4.gridVolumeLoss = nDefectGridS4 * voxelVolume;

acRegion.sector.outside.defectGridIdx = idxDefectGridOutside;
acRegion.sector.outside.nGridPoints = nDefectGridOutside;
acRegion.sector.outside.gridLossPercentDefect = percentDefectGridOutside;
acRegion.sector.outside.gridVolumeLoss = nDefectGridOutside * voxelVolume;

disp(['defect grid points assigned to acetabulum sectors: pelvis defect ', num2str(pelvisNum)]);

end