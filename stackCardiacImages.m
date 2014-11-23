

function stackCardiacImages
  close all; clear;

  datacase = 1;
  switch datacase
    case 1
      load '/Volumes/Seagate2TB/Data/mriCardiacNavs/esp3dn25148.mat';
      % The max dimensions of x are (128,128,64,508);
      x = x / 500;
      [nRows,nCols,nPages,nData] = size(x);
    case 3
      load '/Volumes/Seagate2TB/Data/cardiacRecons/res2mm/img318143.mat';
      % The max dimensions of m are (120,120,80,20);  
      x = m;
      clearvars m;
      [nRows,nCols,nPages,nData] = size(x);
  end
nData = 60;

  midPage = ceil(nPages/2);
  flow = struct( 'du', zeros(nRows,nCols,nPages), ...
                 'dv', zeros(nRows,nCols,nPages), ...
                 'dw', zeros(nRows,nCols,nPages)  );
  flows = cell(nData-1,1);
  for i=1:nData-1
    disp([ 'Working on i=', num2str(i), ' of ', num2str(nData-1) ]);
    flows{i} = flow;

    data1 = abs( x(:,:,:,i) );
    data2 = abs( x(:,:,:,i+1) );

    [du,dv,dw] = opticalFlow3D( data1, data2 );
    flows{i}.du = du;
    flows{i}.dv = dv;
    flows{i}.dw = dw;

    %figure; imshow( imresize( data1(:,:,midPage), 3 ), [] );
    %drawnow;
    %interped = ofInterp3D( data2, du, dv, dw );
    %figure; imshow( imresize( interped(:,:,midPage), 3 ), [] );
    %drawnow;

  end

%clearvars x;
%save( 'stackState1_nIter500.mat', 'flows' );

%load 'stackState1.mat';
load 'stackState1_diffEta.mat';



%   % Make a gif
%   img = zeros(nRows,nCols);
%   midPage = ceil(nPages/2);
%   for i=1:nData-1
%     img(:,:,1) = flows{i}.du(:,:,midPage);
%     img(:,:,2) = flows{i}.dv(:,:,midPage);
%     img(:,:,3) = flows{i}.dw(:,:,midPage);
%     [imind,cm] = rgb2ind(img,256);
%     if i == 1;
%       imwrite(imind,cm,'flows.gif','gif', 'Loopcount',inf);
%     else
%       imwrite(imind,cm,'flows.gif','gif','WriteMode','append');
%     end
%   end

  % Make a stack
  nStack = 3;
  midIndx = 30;
  if mod(nStack,2)==0, error('nStack must be odd'); end;
  nHalfStack = floor( nStack / 2 );
  stack = x(:,:,:,midIndx);
  simpleStack = x(:,:,:,midIndx);
  for i=1:nHalfStack
    fIndx = midIndx + i;
    bIndx = midIndx - i;

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

%close all;
%figure; imshow( imresize( abs(x(:,:,32,midIndx)), 4), [] );
%title('Mid Img');
%figure; imshow( imresize( abs(x(:,:,32,fIndx)), 4), [] );
%title('f Uninterpd');
%figure; imshow( imresize( abs(fInterpd(:,:,32)), 4), [] );
%title('fInterpd');
%figure; imshow( imresize( abs(x(:,:,32,bIndx)), 4), [] );
%title('b Uninterpd');
%figure; imshow( imresize( abs(bInterpd(:,:,32)), 4), [] );
%title('bInterpd');

    stack = stack + fInterpd;
    stack = stack + bInterpd;

    simpleStack = simpleStack + x(:,:,:,fIndx);
    simpleStack = simpleStack + x(:,:,:,bIndx);
  end

  figure;  imshow( imresize( abs(x(:,:,32,midIndx)), 4), [] );
  title('Single Image');
  figure;  imshow( imresize( abs(simpleStack(:,:,32)), 4), [] );
  title('Simple Stack');
  figure;  imshow( imresize( abs(stack(:,:,32)), 4), [] );
  title('OF Stacked Image');

  midPage = ceil(nPages/2);
  midCol = ceil(nCols/2);
  midRow = ceil(nRows/2);


  % Make a difference image
  img30 = x(:,:,:,30);
  img31 = x(:,:,:,31);
  interped31 = ofInterp3D( img31, flows{30}.du, flows{30}.dv, flows{30}.dw );
  absMid30 = abs(img30(:,:,midPage));
  absMid31 = abs(img31(:,:,midPage));
  figure; imshow( imresize(absMid30,4), [] ); title('Img 30');
  figure; imshow( imresize(absMid31,4), [] ); title('Img 31');
  absInt31 = abs(interped31(:,:,midPage));
  figure; imshow( imresize(absInt31,4), [] ); title('Interped 31');
  diff30 = abs( img30 - img31 );
  diffInterped30 = abs( img30 - interped31 );
  axDiff30 = diff30(:,:,midPage);
  axDiffInt30 = diffInterped30(:,:,midPage);
  figure; imshow( imresize(axDiff30,4), [] );  title('Axial Diff 30');
  %figure; imshow( imresize(axDiffInt30,4), [min(axDiff30(:)) max(axDiff30(:))] );
  figure; imshow( imresize(axDiffInt30,4), [0 0.25] );
  sagDiff30 = rot90( squeeze(diff30(:,midCol,:)) );
  sagDiffInt30 = rot90( squeeze(diffInterped30(:,midCol,:)) );
  figure; imshow( imresize(sagDiff30,4), [] );  title('Sagittal Diff 30');
  %figure; imshow( imresize(sagDiffInt30,4), [min(sagDiff30(:)) max(sagDiff30(:))] );
  figure; imshow( imresize(sagDiffInt30,4), [0 0.25] );
  title('Sagittal Diff Warp 30')
  corDiff30 = rot90( squeeze( diff30(midRow,:,:) ) );
  corDiffInt30 = rot90( squeeze( diffInterped30(midRow,:,:) ) );
  figure; imshow( imresize(corDiff30,4), [] );  title('Coronal Diff 30');
  %figure; imshow( imresize(corDiffInt30,4), [min(corDiff30(:)) max(corDiff30(:))] );
  figure; imshow( imresize(corDiffInt30,4), [0 0.25] );
  title('Coronal Diff Warp 30');

  % Make difference images
  diffStep = 1;
  for diffIndx=1:5:21
    imgA = x(:,:,:,diffIndx);
    imgB = x(:,:,:,diffIndx+diffStep);
    interpedB = ofInterp3D( imgB, flows{diffIndx}.du, ...
      flows{diffIndx}.dv, flows{diffIndx}.dw );
    diff = abs( imgA - imgB);
    axDiff = diff(:,:,midPage);
    diffInterped = abs( imgA - interpedB );
    axDiffInterped = diffInterped(:,:,midPage);
    figure; imshow(imresize(axDiff,4),[0 0.25]);
    title(['Ax Diff ', num2str(diffIndx) ]);
    figure; imshow(imresize(axDiffInterped,4), [0 0.25 ]);
      title(['Ax Diff ', num2str(diffIndx), ' Interped' ]);
  end


  % % Select heart portion
  %uLows = zeros(59,1);
  %uHighs = zeros(59,1);
  %vLows = zeros(59,1);
  %vHighs = zeros(59,1);
  %wLows = zeros(59,1);
  %wHighs = zeros(59,1);
  %for i=1:59
  % close all;
  % img = abs(x(:,:,:,i));
  % axImg = squeeze(img(:,:,midPage));
  % figure; imshow( imresize(axImg,4), [] )
  % [u,v] = ginput(1);
  % uLows(i) = round(u/4);
  % vLows(i) = round(v/4);
  % [u,v] = ginput(1);
  % uHighs(i) = round(u/4);
  % vHighs(i) = round(v/4);
  % sagImg = squeeze(img(:,midCol,:));
  % close all;
  % sagImg(uLows(i),:) = max( sagImg(:) );
  % sagImg(uHighs(i),:) = max( sagImg(:) );
  % figure; imshow( imresize(sagImg,4), [] );
  % [w,~] = ginput(1);
  % wLows(i) = round(w/4);
  % [w,~] = ginput(1);
  % wHighs(i) = round(w/4);
  %end
  %close all;
  %save( 'heartRegions.mat', 'uLows', 'uHighs', ...
  % 'vLows', 'vHighs', 'wLows', 'wHighs');
load 'heartRegions.mat';


  % Find average improvement
  improvements = zeros(59,1);
  energyDiffs = zeros(59,1);
  energyInterpeds = zeros(59,1);
  for i=1:59
    imgA = x(:,:,:,i);
    imgB = x(:,:,:,i+1);
    interpedB = ofInterp3D( imgB, flows{i}.du, flows{i}.dv, flows{i}.dw );
    diff = abs( imgA(uLows(i):uHighs(i),vLows(i):vHighs(i),wLows(i):wHighs(i)) - ...
      imgB(uLows(i):uHighs(i),vLows(i):vHighs(i),wLows(i):wHighs(i)) );
    diffInterped = abs( imgA(uLows(i):uHighs(i),vLows(i):vHighs(i),wLows(i):wHighs(i)) - ...
      interpedB(uLows(i):uHighs(i),vLows(i):vHighs(i),wLows(i):wHighs(i)) );
    energyDiffs(i) = norm( diff(:), 2 );
    energyInterpeds(i) = norm( diffInterped(:), 2 );
    improvements(i) = energyDiffs(i) - energyInterpeds(i);
  end
  disp(['Mean Energy Diffs: ', num2str(mean(energyDiffs)) ]);
  disp(['Mean Energy Diff Interps: ', num2str(mean(energyInterpeds)) ]);
  disp(['Mean Improvement: ', num2str(mean(improvements)) ]);

  midImg = abs(x(:,:,:,midIndx));
  du = flows{midIndx}.du;
  dv = flows{midIndx}.dv;
  dw = flows{midIndx}.dw;
  magFlow = sqrt( du.*du + dv.*dv + dw.*dw );
  minMagFlow = min( magFlow(:) );
  maxMagFlow = max( magFlow(:) );

  axMidFlow = zeros(nRows,nCols,3);
  axMidFlow(:,:,1) = du(:,:,midPage);
  axMidFlow(:,:,2) = dv(:,:,midPage);
  axMidFlow(:,:,3) = dw(:,:,midPage);
  axialFlowMag = sqrt( ...
    du(:,:,midPage) .* du(:,:,midPage) + ...
    dv(:,:,midPage) .* dv(:,:,midPage) + ...
    dw(:,:,midPage) .* dw(:,:,midPage) );
  figure; imshow( imresize(midImg(:,:,midPage),4), [] );
  title('Axial Img');
  figure; imshow( imresize(axialFlowMag,4), [] );
  title('Axial Mag Flow');
  figure; imshow( imresize(axMidFlow,4), [minMagFlow maxMagFlow] );
  title('Axial Flow');

  sagFlow = zeros(nRows,nPages,3);
  sagFlow = imrotate( sagFlow, 90 );
  sagFlow(:,:,1) = rot90( squeeze( du(:,midCol,:) ));
  sagFlow(:,:,2) = rot90( squeeze( dv(:,midCol,:) ));
  sagFlow(:,:,3) = rot90( squeeze( dw(:,midCol,:) ));
  sagittalFlowMag = sqrt( ...
    du(:,midCol,:) .* du(:,midCol,:) + ...
    dv(:,midCol,:) .* dv(:,midCol,:) + ...
    dw(:,midCol,:) .* dw(:,midCol,:) );
  sagittalFlowMag = rot90( squeeze( sagittalFlowMag ) );
  midColImg = squeeze( midImg(:,midCol,:) );
  midColImg = rot90( midColImg );
  figure; imshow( imresize(midColImg,4), [] );
  title('Sagittal Img');
  figure; imshow( imresize(sagittalFlowMag,4), [minMagFlow maxMagFlow] );
  title('Sagittal Mag Flow');
  figure; imshow( imresize(sagFlow,4), [] );
  title('Sagittal Flow')

  corFlow = zeros(nCols,nPages,3);
  corFlow = imrotate( corFlow, 90 );
  corFlow(:,:,1) = rot90( squeeze( du(midRow,:,:) ));
  corFlow(:,:,2) = rot90( squeeze( dv(midRow,:,:) ));
  corFlow(:,:,3) = rot90( squeeze( dw(midRow,:,:) ));
  coronalFlowMag = sqrt( ...
    du(midRow,:,:) .* du(midRow,:,:) + ...
    dv(midRow,:,:) .* dv(midRow,:,:) + ...
    dw(midRow,:,:) .* dw(midRow,:,:) );
  coronalFlowMag = rot90( squeeze( coronalFlowMag ) );
  midRowImg = squeeze( midImg(midRow,:,:) );
  midRowImg = rot90( midRowImg );
  figure; imshow( imresize(midRowImg,4), [] );
  title('Coronal Img');
  figure; imshow( imresize(coronalFlowMag,4), [minMagFlow maxMagFlow] );
  title('Coronal Mag Flow');
  figure; imshow( imresize(corFlow,4), [] );
  title('Coronal Flow');

end











