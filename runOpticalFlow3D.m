
function runOpticalFlow3D
  clear; close all;


  datacase = 5;
  x = loadData( datacase );
  data1Indx = 2;
  data2Indx = 3;

  data1 = squeeze(x(:,:,:,data1Indx));
  data2 = squeeze(x(:,:,:,data2Indx));

  [nRows,nCols,nSlices] = size(data1);
  midSlice = floor(nSlices/2);

  mask = loadMask( [nRows,nCols,nSlices], datacase );

  maskedIndxs = find( mask>0 );
  scales = data1(maskedIndxs) ./ data2(maskedIndxs);
  data2 = data2 .* median( scales );
  
  

  %profile on
  tic;
  [du,dv,dw] = opticalFlow3D( data1, data2, mask );
  timeTaken = toc;
  %profile off
  %profile viewer

  interped = ofInterp3D( data2, du, dv, dw );
  close all;
  img1 = squeeze(x(:,:,midSlice,data1Indx));
  img2 = squeeze(x(:,:,midSlice,data2Indx));
  figure, imshow( imresize(img1,3), [] );
  title('img1', 'FontSize', 20 );
  interpedImg = interped(:,:,midSlice);
  figure, imshow( imresize(interpedImg,3), [min(data1(:)) max(data1(:))] );
  title('interped', 'FontSize', 20 );
  figure, imshow( imresize(img2,3), [min(data1(:)) max(data1(:))] );
  title('img2', 'FontSize', 20 );

  disp([ 'Time taken (s): ', num2str(timeTaken) ]);
end
