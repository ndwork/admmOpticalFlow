
function runOpticalFlow3D

  %data1 = 
  %data2 = 

  
  tic;
  profile on
  [du,dv] = opticalFlow3D( data1, data2 );
  profile off
  timeTaken = toc;
  profile viewer

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