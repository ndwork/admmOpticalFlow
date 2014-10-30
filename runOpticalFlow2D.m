

function runOpticalFlow2D

  %img1 = double( imread( '../data/basketball/frame13.png' ) ) / 255;
  %img2 = circshift( img1, [2 2] );
  %img2 = double( imread( '../data/basketball/frame14.png' ) ) / 255;
  %img1 = double( imread( '../data/Yosemite/frame10.png' ) ) / 255;
  %img2 = double( imread( '../data/Yosemite/frame11.png' ) ) / 255;
  img1 = double( imread( '../Dumptruck/frame10.png' ) ) / 255;
  img2 = double( imread( '../Dumptruck/frame11.png' ) ) / 255;

  tic;
  profile on
  [du,dv] = opticalFlow2D( img1, img2 );
  profile off
  timeTaken = toc;
  profile viewer

  interped2 = ofInterp( img2, du, dv );
  close all;
  figure, imshow( img1, [] );
  title('img1', 'FontSize', 20 );
  figure, imshow( interped2, [min(img1(:)) max(img1(:))] );
  title('interped2', 'FontSize', 20 );
  figure, imshow( img2, [] );
  title('img2', 'FontSize', 20 );

  disp([ 'Time taken (s): ', num2str(timeTaken) ]);
end

