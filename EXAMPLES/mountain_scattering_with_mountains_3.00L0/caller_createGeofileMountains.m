clear all;
close all;
clc;

%%%

% without_0__with033L_1__with300L_2 = 0;
% without_0__with033L_1__with300L_2 = 1;
% without_0__with033L_1__with300L_2 = 2;

zmax_m = 15e3;
xmax_peaks = 18e3;
xmax = 23e3;

% Model.
% % down to mantle
% hLayers_m = [3600, 12960, 14580, 12960, 4000];
% dxGrdInts = [516, 762, 810, 864, 1000];
% % down to lower crust
% hLayers_m = [3600, 12960, 14580, 10000];
% dxGrdInts = [516, 762, 810, 864];
% remove middle sediments, go down to lower crust
hLayers_m = [16560, 14580, 10000];
dxGrdInts = [762, 810*2, 864*3];

dxAir = [1, 1.2]*110;
valleysAltitude = 0;
meshAlgorithm = 5; % delaunay

% for without_0__with033L_1__with300L_2 = (0:2)
for without_0__with033L_1__with300L_2 = 0

  switch(without_0__with033L_1__with300L_2)
    case 0
      % output to "without" folder
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_without_mountains/EXTMSH/extMesh.geo';
      nPeaks = 0;
      peaksHalfWidth = 2000;
      peaksAltitude = 0;
    case 1
      % output to "with 0.33L" folder
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_0.33L0/EXTMSH/extMesh.geo';
      peaksHalfWidth = 1000/2;
      nPeaks = xmax_peaks/peaksHalfWidth;
      peaksAltitude = 1500;
      meshAlgorithm = 9; % packing of parallelograms
    case 2
      % output to "with 3.00L" folder
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_3.00L0/EXTMSH/extMesh.geo';
      peaksHalfWidth = 9000/2;
      nPeaks = xmax_peaks/peaksHalfWidth;
      peaksAltitude = 1500;
  end
  outputFile

  distanceBetweenPeaks = peaksHalfWidth*ones(1, nPeaks-1);
  lastPeakHalfWidth = peaksHalfWidth;
  mHalfWidthLeftPeaks = peaksHalfWidth*ones(1, nPeaks);

  mPeakHeight = peaksAltitude*ones(1, nPeaks); % height of the peaks, left to right
  mBaseLeftHeight = [0, valleysAltitude*ones(1, nPeaks-1)]; % height of the left base of the peak for each peak, left to right

  zb_m = [-sum(hLayers_m), 0, zmax_m];

  hLayersProportion = flip(hLayers_m)/abs(zb_m(1)); % bottom to top

  x_peakToPeak_right = 2*[distanceBetweenPeaks, lastPeakHalfWidth]; % left to right
  mwlp = 2*mHalfWidthLeftPeaks;
  mwrp = flip(x_peakToPeak_right);

  nLayersGround = numel(hLayers_m);

  if(xmax < 0.5*nPeaks*2*peaksHalfWidth)
    error('xmax too small to contain all peaks');
  end
  xb_m = [-1,1]*xmax;

  distRight = max(xb_m)-peaksHalfWidth*nPeaks;

  createGeofileMountains(outputFile, meshAlgorithm, xb_m, zb_m, nLayersGround, hLayersProportion, nPeaks, mwrp, mwlp, distRight, flip(mPeakHeight), flip(mBaseLeftHeight), dxAir, flip(dxGrdInts));
end