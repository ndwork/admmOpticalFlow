
function runOpticalFlow3D
  clear; close all;


  datacase = 4;
  x = loadData( datacase );
  data1Indx = 1;
  data2Indx = 2;
  
  data1 = squeeze(x(:,:,:,data1Indx));
  data2 = squeeze(x(:,:,:,data2Indx));

  [~,~,nPages] = size( data1 );
  midPage = floor(nPages/2);

  %figure;  imshow( imresize( abs(img1), 5 ), [] );
  %title('Img 1');
  %figure;  imshow( imresize( abs(img2), 5 ), [] );
  %title('Img 2');


  %profile on
  tic;
  [du,dv,dw] = opticalFlow3D( data1, data2 );
  timeTaken = toc;
  %profile off
  %profile viewer

  interped = ofInterp3D( data2, du, dv, dw );
  close all;
  img1 = squeeze(x(:,:,midPage,data1Indx));
  img2 = squeeze(x(:,:,midPage,data2Indx));
  figure, imshow( imresize(img1,3), [] );
  title('img1', 'FontSize', 20 );
  interpedImg = interped(:,:,midPage);
  figure, imshow( imresize(interpedImg,3), [min(data1(:)) max(data1(:))] );
  title('interped', 'FontSize', 20 );
  figure, imshow( imresize(img2,3), [min(data1(:)) max(data1(:))] );
  title('img2', 'FontSize', 20 );

  disp([ 'Time taken (s): ', num2str(timeTaken) ]);
end
