
function runOpticalFlow3D
  clear; close all;

  %load '../esp3dn25148.mat';
  % The max dimensions of x are (128,128,64,508);

  %data1 = abs( x( :, :, :, 250 ) );
  %data2 = abs( x( :, :, :, 251 ) );
load 'junk.mat';

  img1 = data1(:,:,32);
  %figure;  imshow( imresize( abs(img1), 5 ), [] );
  %title('Data 1');
  img2 = data2(:,:,32);
  %figure;  imshow( imresize( abs(img2), 5 ), [] );
  %title('Data 2');

  
  tic;
  profile on
  [du,dv] = opticalFlow3D( data1, data2 );
  profile off
  timeTaken = toc;
  profile viewer

save( 'ofResults.mat' );
  
  interped = ofInterp3D( data2, du, dv, dw );
  close all;
  figure, imshow( data1, [] );
  title('img1', 'FontSize', 20 );
  figure, imshow( interped, [min(data1(:)) max(data1(:))] );
  title('interped2', 'FontSize', 20 );
  figure, imshow( data2, [min(data1(:)) max(data1(:))] );
  title('img2', 'FontSize', 20 );

  disp([ 'Time taken (s): ', num2str(timeTaken) ]);
end