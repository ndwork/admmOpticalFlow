
function runOpticalFlow3D
  clear; close all;

  datacase = 9;
  [x,eta,rho,offset] = loadData( datacase );
  data1Indx = offset+1;
  data2Indx = data1Indx+3;
%data1Indx = 4;
%data2Indx = 10;

  bound = 1.2;

  data1 = squeeze(x(:,:,:,data1Indx));
  data2 = squeeze(x(:,:,:,data2Indx));

  %mask = loadMask( [nRows,nCols,nSlices], datacase );

  profile on
  if exist('mask','var') > 0
    maskedIndxs = find( mask>0 );
    scales = data1(maskedIndxs) ./ data2(maskedIndxs);
    data2 = data2 .* median( scales(:) );
    tic;
    [du,dv,dw] = opticalFlow3D( data1, data2, eta, rho, mask );
    timeTaken = toc;
  else
    maskedIndxs = find( data2 ~= 0 );
    scales = data1(maskedIndxs) ./ data2(maskedIndxs);
    data2 = data2 .* median( scales(:) );
    tic;
    if exist( 'bound', 'var' )
      [du,dv,dw] = boundedOpticalFlow3D( data1, data2, eta, rho, bound );
    else
      [du,dv,dw] = opticalFlow3D( data1, data2, eta, rho );
    end
    timeTaken = toc;
  end
  profile off
  profile viewer

  [nRows,nCols,nSlices] = size(data1);
  midRow = ceil(nRows/2);
  midCol = ceil(nCols/2);
  midSlice = ceil(nSlices/2);

  imgScale = 10;
  decimate = 2;
  interped = ofInterp3D( data2, du, dv, dw );

  makeArrowImage2( imgScale, decimate, ...
    squeeze(x(midRow,20:60,20:60,data1Indx)), squeeze(x(midRow,20:60,20:60,data2Indx)), ...
    squeeze(dv(midRow,20:60,20:60)), squeeze(dw(midRow,20:60,20:60)), ...
    squeeze(interped(midRow,20:60,20:60)), 1 );

  makeArrowImage2( imgScale, decimate, ...
    squeeze(x(20:60,midCol,20:60,data1Indx)), squeeze(x(20:60,midCol,20:60,data2Indx)), ...
    squeeze(dv(20:60,midCol,20:60)), squeeze(dw(20:60,midCol,20:60)), ...
    squeeze(interped(20:60,midCol,20:60)), 1 );

  makeArrowImage2( imgScale, decimate, ...
    squeeze(x(20:60,20:60,midSlice,data1Indx)), squeeze(x(20:60,20:60,midSlice,data2Indx)), ...
    squeeze(dv(20:60,20:60,midSlice)), squeeze(dw(20:60,20:60,midSlice)), ...
    squeeze(interped(20:60,20:60,midSlice)), 0 );
  
  sliceImg1 = squeeze(x(:,:,midSlice,data1Indx) );
  sliceImg2 = squeeze(x(:,:,midSlice,data2Indx) );
  sliceInterped = squeeze( interped(:,:,midSlice) );
  figure, imshow( imresize(sliceImg1,5), [min(data1(:)) max(data1(:))] );
  title('img1', 'FontSize', 20 );
  figure, imshow( imresize(sliceInterped,5), [min(data1(:)) max(data1(:))] );
  title('interped', 'FontSize', 20 );
  figure, imshow( imresize(sliceImg2,5), [min(data1(:)) max(data1(:))] );
  title('img2', 'FontSize', 20 );

  disp([ 'Time taken (s): ', num2str(timeTaken) ]);
end


function makeArrowImage( imgScale, decimate, ...
  img1, img2, dx, dy, interped, rotateIt )
  [nRows, nCols] = size(img1);
  if nargin < 4, rotateIt=0; end;
  img1Scaled = imresize( img1, imgScale, 'bilinear' );
  img2Scaled = imresize( img2, imgScale, 'bilinear' );
  imgInterpedScaled = imresize( interped, imgScale, 'bilinear' );
  if rotateIt > 0
    img1Scaled = rot90( img1Scaled, rotateIt );
    img2Scaled = rot90( img2Scaled, rotateIt );
    imgInterpedScaled = rot90( imgInterpedScaled, rotateIt );
  end
  quivH = figure;
  of2Quiver( dx, dy, decimate, 0, 10 );
  axis([0 nRows 0 nCols]);
  set(gca,'position',[0 0 1 1],'units','normalized')
  set(gca,'xtick',[]);  set(gca,'ytick',[]);
  quivImg = getframe;
  quivImg = double( rgb2gray( quivImg.cdata ) ) / 255.0;
  quivImg = imresize( quivImg, ...
    [nCols*imgScale, nCols*imgScale], 'nearest' );
  quivImg(:,1)=1;  quivImg(1,:)=1;
  quivImg = imdilate( quivImg ~= 1, [ 1 1; 1 0] );
  if rotateIt > 0, quivImg = rot90( quivImg); end;
  img2Scaled( quivImg == 1 ) = max( img2Scaled(:) ) * 1.4;
  img2Scaled(:,1) = 0;  img2Scaled(1,:) = 0;
  arrowImg = [ img1Scaled img2Scaled imgInterpedScaled ];
  arrowImg = arrowImg / max( arrowImg(:) );
  close( quivH );
  %figure; imshow( rowImg, [min(data1(:)) max(data1(:))] );
  figure; imshow( arrowImg, [] );
end


function makeArrowImage2( imgScale, decimate, ...
  img1, img2, dx, dy, interped, rotateIt )
  [nRows, nCols] = size(img1);
  if nargin < 4, rotateIt=0; end;
  img1Scaled = imresize( img1, imgScale, 'bilinear' );
  img2Scaled = imresize( img2, imgScale, 'bilinear' );
  imgInterpedScaled = imresize( interped, imgScale, 'bilinear' );
  if rotateIt > 0
    img1Scaled = rot90( img1Scaled, rotateIt );
    img2Scaled = rot90( img2Scaled, rotateIt );
    imgInterpedScaled = rot90( imgInterpedScaled, rotateIt );
  end
  quivH = figure;
  of2Quiver( dx, dy, decimate, 0, 10 );
  axis([0 nRows 0 nCols]);
  set(gca,'position',[0 0 1 1],'units','normalized')
  set(gca,'xtick',[]);  set(gca,'ytick',[]);
  quivImg = getframe;
  quivImg = double( rgb2gray( quivImg.cdata ) ) / 255.0;
  quivImg = imresize( quivImg, ...
    [nCols*imgScale, nCols*imgScale], 'nearest' );
  quivImg(:,1)=1;  quivImg(1,:)=1;
  quivImg = imdilate( quivImg ~= 1, [ 1 1; 1 0] );
  if rotateIt > 0, quivImg = rot90( quivImg); end;
  img2Quiv = img2Scaled;
  img2Quiv( quivImg == 1 ) = max( img2Quiv(:) ) * 1.4;
  img2Quiv(:,1) = 0;  img2Quiv(1,:) = 0;
  arrowImg = [ img1Scaled img2Scaled imgInterpedScaled img2Quiv ];
  arrowImg = arrowImg / max( arrowImg(:) );
  close( quivH );
  %figure; imshow( rowImg, [min(data1(:)) max(data1(:))] );
  figure; imshow( arrowImg, [] );
end
  
  