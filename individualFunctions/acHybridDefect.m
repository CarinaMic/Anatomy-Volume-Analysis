%% Acetabulum region: defect points inside hybrid regions

% Input:    acRegion: Struct with acetabulum hybrid region to be updated
%           isPointInsideDefect: Logical mask over all grid points / pointsBox
%                   - true for grid points inside defect
%                   - false for grid points outside defect
%           refHybrid: Struct with reference acetabulum hybrid regions
%                   - mc1.gridIdx: grid point indices in hybrid region mc1
%                   - mc2.gridIdx: grid point indices in hybrid region mc2
%                   - mc3.gridIdx: grid point indices in hybrid region mc3
%                   - mc4.gridIdx: grid point indices in hybrid region mc4
%                   - mc5.gridIdx: grid point indices in hybrid region mc5
%           voxelSize: Voxel spacing in x-, y-, and z-direction
%           pelvisNum: Numeric identifier used only for logging
%
% Output:   acRegion: Struct with updated acetabulum hybrid defect results:
%                   - hybrid.totalDefect.defectGridIdx: all defect grid point indices
%                   - hybrid.totalDefect.nGridPoints: total number of defect grid points
%                   - hybrid.totalDefect.gridVolumeLoss: total defect volume loss
%                   - hybrid.mc1.defectGridIdx: defect grid point indices in mc1
%                   - hybrid.mc2.defectGridIdx: defect grid point indices in mc2
%                   - hybrid.mc3.defectGridIdx: defect grid point indices in mc3
%                   - hybrid.mc4.defectGridIdx: defect grid point indices in mc4
%                   - hybrid.mc5.defectGridIdx: defect grid point indices in mc5
%                   - hybrid.mc1.nGridPoints: number of defect grid points in mc1
%                   - hybrid.mc2.nGridPoints: number of defect grid points in mc2
%                   - hybrid.mc3.nGridPoints: number of defect grid points in mc3
%                   - hybrid.mc4.nGridPoints: number of defect grid points in mc4
%                   - hybrid.mc5.nGridPoints: number of defect grid points in mc5
%                   - hybrid.mc1.gridLossPercentDefect: percent of total defect in mc1
%                   - hybrid.mc2.gridLossPercentDefect: percent of total defect in mc2
%                   - hybrid.mc3.gridLossPercentDefect: percent of total defect in mc3
%                   - hybrid.mc4.gridLossPercentDefect: percent of total defect in mc4
%                   - hybrid.mc5.gridLossPercentDefect: percent of total defect in mc5
%                   - hybrid.mc1.gridLossPercentRegion: percent loss relative to reference mc1
%                   - hybrid.mc2.gridLossPercentRegion: percent loss relative to reference mc2
%                   - hybrid.mc3.gridLossPercentRegion: percent loss relative to reference mc3
%                   - hybrid.mc4.gridLossPercentRegion: percent loss relative to reference mc4
%                   - hybrid.mc5.gridLossPercentRegion: percent loss relative to reference mc5
%                   - hybrid.mc1.gridVolumeLoss: volume loss in mc1
%                   - hybrid.mc2.gridVolumeLoss: volume loss in mc2
%                   - hybrid.mc3.gridVolumeLoss: volume loss in mc3
%                   - hybrid.mc4.gridVolumeLoss: volume loss in mc4
%                   - hybrid.mc5.gridVolumeLoss: volume loss in mc5
%                   - hybrid.outside.defectGridIdx: defect grid points outside all hybrid regions
%                   - hybrid.outside.nGridPoints: number of outside defect grid points
%                   - hybrid.outside.gridLossPercentDefect: percent of total defect outside hybrid regions
%                   - hybrid.outside.gridVolumeLoss: outside defect volume loss

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function acRegion = acHybridDefect(acRegion, isPointInsideDefect, refHybrid, voxelSize, pelvisNum)

% Hybrid regions of reference pelvis
mc1.gridIdx = refHybrid.mc1.gridIdx;
mc2.gridIdx = refHybrid.mc2.gridIdx;
mc3.gridIdx = refHybrid.mc3.gridIdx;
mc4.gridIdx = refHybrid.mc4.gridIdx;
mc5.gridIdx = refHybrid.mc5.gridIdx;

% Number of reference grid points per region
nRefMC1 = numel(mc1.gridIdx);
nRefMC2 = numel(mc2.gridIdx);
nRefMC3 = numel(mc3.gridIdx);
nRefMC4 = numel(mc4.gridIdx);
nRefMC5 = numel(mc5.gridIdx); %%%

% Voxel size [mm]
voxelVolume = prod(voxelSize);

% Global indices of defect grid points
idxDefectGrid = find(isPointInsideDefect);

% Total number of defect grid points
nDefectGrid_total = numel(idxDefectGrid);

% Defect grid points assigned to each hybrid region
idxDefectMC1 = intersect(idxDefectGrid, mc1.gridIdx);
idxDefectMC2 = intersect(idxDefectGrid, mc2.gridIdx);
idxDefectMC3 = intersect(idxDefectGrid, mc3.gridIdx);
idxDefectMC4 = intersect(idxDefectGrid, mc4.gridIdx);
idxDefectMC5 = intersect(idxDefectGrid, mc5.gridIdx);

% Defect grid points outside all hybrid regions
idxDefectAssigned = [idxDefectMC1; idxDefectMC2; idxDefectMC3; idxDefectMC4; idxDefectMC5];
idxDefectOutside = setdiff(idxDefectGrid, idxDefectAssigned);

% Number of defect grid points per region
nDefectMC1 = numel(idxDefectMC1);
nDefectMC2 = numel(idxDefectMC2);
nDefectMC3 = numel(idxDefectMC3);
nDefectMC4 = numel(idxDefectMC4);
nDefectMC5 = numel(idxDefectMC5);
nDefectOutside = numel(idxDefectOutside);

% Proportion of total defect grid points
percentDefectMC1 = 100 * nDefectMC1 / nDefectGrid_total;
percentDefectMC2 = 100 * nDefectMC2 / nDefectGrid_total;
percentDefectMC3 = 100 * nDefectMC3 / nDefectGrid_total;
percentDefectMC4 = 100 * nDefectMC4 / nDefectGrid_total;
percentDefectMC5 = 100 * nDefectMC5 / nDefectGrid_total;
percentDefectOutside = 100 * nDefectOutside / nDefectGrid_total;

% Loss percent relative to reference pelvis grid points inside each region
percentRegionMC1 = 100 * nDefectMC1 / nRefMC1;
percentRegionMC2 = 100 * nDefectMC2 / nRefMC2;
percentRegionMC3 = 100 * nDefectMC3 / nRefMC3;
percentRegionMC4 = 100 * nDefectMC4 / nRefMC4;
percentRegionMC5 = 100 * nDefectMC5 / nRefMC5;

% Save
acRegion.hybrid.totalDefect.defectGridIdx = idxDefectGrid;
acRegion.hybrid.totalDefect.nGridPoints = nDefectGrid_total;
acRegion.hybrid.totalDefect.gridVolumeLoss = nDefectGrid_total * voxelVolume;

acRegion.hybrid.mc1.defectGridIdx = idxDefectMC1;
acRegion.hybrid.mc1.nGridPoints = nDefectMC1;
acRegion.hybrid.mc1.gridLossPercentDefect = percentDefectMC1;
acRegion.hybrid.mc1.gridLossPercentRegion = percentRegionMC1;
acRegion.hybrid.mc1.gridVolumeLoss = nDefectMC1 * voxelVolume;

acRegion.hybrid.mc2.defectGridIdx = idxDefectMC2;
acRegion.hybrid.mc2.nGridPoints = nDefectMC2;
acRegion.hybrid.mc2.gridLossPercentDefect = percentDefectMC2;
acRegion.hybrid.mc2.gridLossPercentRegion = percentRegionMC2;
acRegion.hybrid.mc2.gridVolumeLoss = nDefectMC2 * voxelVolume;

acRegion.hybrid.mc3.defectGridIdx = idxDefectMC3;
acRegion.hybrid.mc3.nGridPoints = nDefectMC3;
acRegion.hybrid.mc3.gridLossPercentDefect = percentDefectMC3;
acRegion.hybrid.mc3.gridLossPercentRegion = percentRegionMC3;
acRegion.hybrid.mc3.gridVolumeLoss = nDefectMC3 * voxelVolume;

acRegion.hybrid.mc4.defectGridIdx = idxDefectMC4;
acRegion.hybrid.mc4.nGridPoints = nDefectMC4;
acRegion.hybrid.mc4.gridLossPercentDefect = percentDefectMC4;
acRegion.hybrid.mc4.gridLossPercentRegion = percentRegionMC4;
acRegion.hybrid.mc4.gridVolumeLoss = nDefectMC4 * voxelVolume;

acRegion.hybrid.mc5.defectGridIdx = idxDefectMC5;
acRegion.hybrid.mc5.nGridPoints = nDefectMC5;
acRegion.hybrid.mc5.gridLossPercentDefect = percentDefectMC5;
acRegion.hybrid.mc5.gridLossPercentRegion = percentRegionMC5;
acRegion.hybrid.mc5.gridVolumeLoss = nDefectMC5 * voxelVolume;

acRegion.hybrid.outside.defectGridIdx = idxDefectOutside;
acRegion.hybrid.outside.nGridPoints = nDefectOutside;
acRegion.hybrid.outside.gridLossPercentDefect = percentDefectOutside;
acRegion.hybrid.outside.gridVolumeLoss = nDefectOutside * voxelVolume;

disp(['defect grid points assigned to hybrid acetabulum regions: pelvis defect ', num2str(pelvisNum)]);
end
