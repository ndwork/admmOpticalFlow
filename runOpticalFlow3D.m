
function runOpticalFlow3D
  clear; close all;

  
  datacase = 1;
  switch datacase
    case 1
      %load '/Volumes/Seagate2TB/Data/mriCardiacNavs/esp3dn25148.mat';
      % The max dimensions of x are (128,128,64,508);

      %data1 = abs( x( :, :, :, 250 ) ) / 500;
      %data2 = abs( x( :, :, :, 251 ) ) / 500;
      %save( 'junk.mat', 'data1', 'data2' );
      load 'junk.mat';

      img1 = data1(:,:,32);
      img2 = data2(:,:,32);
      
    case 2
      % Not very much of a heartbeat
      load '/Volumes/Seagate2TB/Data/cardiacRecons/res1.2mm/img32013.mat';
      % The max dimensions of m are (256,256,112,6);
      data1 = abs( m( :, :, :, 1 ) );
      data2 = abs( m( :, :, :, 6 ) );
      
      
    case 3
      load '/Volumes/Seagate2TB/Data/cardiacRecons/res2mm/img318143.mat';
      % The max dimensions of m are (120,120,80,20);
      data1 = abs( m( :, :, :, 10 ) );
      data2 = abs( m( :, :, :, 12 ) );
      
      img1 = data1(:,:,40);
      img2 = data2(:,:,40);
  end


  %figure;  imshow( imresize( abs(img1), 5 ), [] );
  %title('Img 1');
  %figure;  imshow( imresize( abs(img2), 5 ), [] );
  %title('Img 2');


  tic;
  profile on
  [du,dv,dw] = opticalFlow3D( data1, data2 );
  profile off
  timeTaken = toc;
  profile viewer

save( 'ofResults.mat' );

  interped = ofInterp3D( data2, du, dv, dw );
  close all;
  figure, imshow( imresize(img1,3), [] );
  title('img1', 'FontSize', 20 );
  interpedImg = interped(:,:,32);
  figure, imshow( imresize(interpedImg,3), [min(data1(:)) max(data1(:))] );
  title('interped', 'FontSize', 20 );
  figure, imshow( imresize(img2,3), [min(data1(:)) max(data1(:))] );
  title('img2', 'FontSize', 20 );

  disp([ 'Time taken (s): ', num2str(timeTaken) ]);
end
