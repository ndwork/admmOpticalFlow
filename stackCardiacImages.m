
function stackCardiacImages
  close all; clear;

  %-- Parameters
  datacase = 10;
  nStack = 3;
  maxNData = 7;

  [x,eta,rho] = loadData( datacase );
  [nRows,nCols,nSlices,nData] = size(x);
  nData = min([nData,maxNData]);


  %mask = loadMask( [nRows,nCols,nSlices], datacase );

  nStack = min([nStack,nData]);

  midData = ceil(nData/2);
  midSlice = ceil(nSlices/2);
  flow = struct( 'du', zeros(nRows,nCols,nSlices), ...
                 'dv', zeros(nRows,nCols,nSlices), ...
                 'dw', zeros(nRows,nCols,nSlices)  );
  flows = cell(nData-1,1);

  parfor i=1:nData-1
    disp([ 'Working on i=', num2str(i), ' of ', num2str(nData-1) ]);
    flows{i} = flow;

    data1 = squeeze( x(:,:,:,i) );
    data2 = squeeze( x(:,:,:,i+1) );

    if exist('mask','var')
      maskedIndxs = find( mask>0 & data2~=0 );
      scales = data1(maskedIndxs) ./ data2(maskedIndxs);
      data2 = data2 .* median( scales );
      [du,dv,dw] = opticalFlow3D( data1, data2, eta, rho, mask );
    else
      maskedIndxs = find( data2~=0 );
      scales = data1(maskedIndxs) ./ data2(maskedIndxs);
      data2 = data2 .* median( scales );
      [du,dv,dw] = opticalFlow3D( data1, data2, eta, rho );
    end
    flows{i}.du = du;
    flows{i}.dv = dv;
    flows{i}.dw = dw;
  end

save( 'stackState.mat' );
%load 'stackState_1dn3.mat';   % eta=1d-3
%load 'stackState_5dn4.mat';   % eta=5d-4
%load 'stackState_1d10scaleByN.mat';


  % Make a stack
  if mod(nStack,2)==0, error('nStack must be odd'); end;
  nHalfStack = floor( nStack / 2 );
  stack = x(:,:,:,midData);
  simpleStack = x(:,:,:,midData);
figure;
%subplot(1,2,1);  imshow( imresize( x(:,:,midPage,midData), 3), [] );
%subplot(1,2,2);
imshow( imresize( x(:,:,midSlice,midData), 3), [] );
title(num2str(midData), 'FontSize', 20);
  for i=1:nHalfStack
    fIndx = midData + i;
    bIndx = midData - i;

    if i==1
      fFlows = flows{fIndx-1};
      bFlows = flows{bIndx};
      bFlows.du = -bFlows.du;
      bFlows.dv = -bFlows.dv;
      bFlows.dw = -bFlows.dw;
    else
      fFlows.du = fFlows.du + flows{fIndx-1}.du;
      fFlows.dv = fFlows.dv + flows{fIndx-1}.dv;
      fFlows.dw = fFlows.dw + flows{fIndx-1}.dw;
      bFlows.du = bFlows.du - flows{bIndx}.du;
      bFlows.dv = bFlows.dv - flows{bIndx}.dv;
      bFlows.dw = bFlows.dw - flows{bIndx}.dw;
    end

    fInterpd = ofInterp3D( x(:,:,:,fIndx), fFlows.du, fFlows.dv, fFlows.dw );    
    bInterpd = ofInterp3D( x(:,:,:,bIndx), bFlows.du, bFlows.dv, bFlows.dw );

figure; imshow( imresize( squeeze(x(125,:,:,fIndx)), 3), [] );
title([num2str(fIndx), ' - Orig'], 'FontSize', 20);
figure; imshow( imresize( squeeze(fInterpd(125,:,:)), 3), [] );
title(num2str(fIndx), 'FontSize', 20);
figure; imshow( imresize( squeeze(x(125,:,:,midData)), 3), [] );
title('Mid Data', 'Fontsize', 20);
figure; imshow( imresize( squeeze(x(125,:,:,bIndx)), 3), [] );
title([num2str(bIndx), ' - Orig'], 'FontSize', 20);
figure; imshow( imresize( squeeze(bInterpd(125,:,:)), 3), [] );
title(num2str(bIndx), 'FontSize', 20);

magFFlows = sqrt( fFlows.du.^2 + fFlows.dv.^2 + fFlows.dw.^2 );
figure;  imshow( imresize( magFFlows(:,:,midSlice), 3), [] );
title(num2str(fIndx), 'FontSize', 20);
magBFlows = sqrt( bFlows.du.^2 + bFlows.dv.^2 + bFlows.dw.^2 );
figure;  imshow( imresize( magBFlows(:,:,midSlice), 3), [] );
title(num2str(bIndx), 'FontSize', 20);

    stack = stack + fInterpd;
    stack = stack + bInterpd;

    simpleStack = simpleStack + x(:,:,:,fIndx);
    simpleStack = simpleStack + x(:,:,:,bIndx);
  end

  figure;  imshow( imresize( abs(x(:,:,midSlice,midData)), 3), [] );
  title('Single Image', 'FontSize', 20);
  figure;  imshow( imresize( abs(simpleStack(:,:,midSlice)), 3), [] );
  title('Simple Stack', 'FontSize', 20);
  figure;  imshow( imresize( abs(stack(:,:,midSlice)), 3), [] );
  title('OF Stacked Image', 'FontSize', 20);

  figure;  imshow( imresize( abs(squeeze(x(125,:,:,midData+1)))', 3), [] );
  title('Single Image-Mid-1', 'FontSize', 20);
  figure;  imshow( imresize( abs(squeeze(x(125,:,:,midData)))', 3), [] );
  title('Single Image-Mid', 'FontSize', 20);
  figure;  imshow( imresize( abs(squeeze(x(125,:,:,midData-1)))', 3), [] );
  title('Single Image-Mid+1', 'FontSize', 20);
  figure;  imshow( imresize( abs(squeeze(simpleStack(125,:,:)))', 3), [] );
  title('Simple Stack', 'FontSize', 20);
  figure;  imshow( imresize( abs(squeeze(stack(125,:,:)))', 3), [] );
  title('OF Stacked Image', 'FontSize', 20);

  midCol = ceil(nCols/2);
  figure;  imshow( imresize( abs(squeeze(x(:,125,:,midData)))', 3), [] );
  title('Single Image', 'FontSize', 20);
  figure;  imshow( imresize( abs(squeeze(simpleStack(:,125,:)))', 3), [] );
  title('Simple Stack', 'FontSize', 20);
  figure;  imshow( imresize( abs(squeeze(stack(:,125,:)))', 3), [] );
  title('OF Stacked Image', 'FontSize', 20);
  
  %saveAllSlices( stack, 1 );

end

