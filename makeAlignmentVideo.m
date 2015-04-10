
function makeAlignmentVideo
  close all; clear;

  maxNData = 100;

  %datacase = 1;
  datacase = 9;
  %datacase = 11;

  orientation = 1;   % 0 - 3x2, 1 - 2x3

  outDir = ['alignImgs_',num2str(datacase)];

  [x,eta,rho,offset] = loadData( datacase );
  x = abs(x);
  [nRows,nCols,nPages,nData] = size(x);
  midRow = ceil(nRows/2);
  midCol = ceil(nCols/2);
  midSlice = ceil(nPages/2);
  nData = min([offset+nData,maxNData,nData]);

  bound = 1.2;

  mkdir( outDir );
  data1 = squeeze( x(:,:,:,offset+1) );
  parfor i=1:nData-offset
    disp(['Working on image ', num2str(i)]);
    data2 = squeeze( x(:,:,:,offset+i) );

    maskedIndxs = find( data2 ~= 0 );
    scales = data1(maskedIndxs) ./ data2(maskedIndxs);
    data2scaled = data2 .* median( scales );

    %[du,dv,dw] = opticalFlow3D( data1, data2scaled, eta, rho );
    [du,dv,dw] = boundedOpticalFlow3D( data1, data2, eta, rho, bound );
    interp2 = ofInterp3D( data2, du, dv, dw );

    axial2 = data2(:,:,midSlice);
    axialInterp = interp2(:,:,midSlice);
    if orientation==0
      axialImg2Save = [ axial2 axialInterp ];
    else
      axialImg2Save = [ axial2; axialInterp; ];
    end
    [nRowsAx, nColsAx] = size(axialImg2Save);

    cor2 = rot90( squeeze( data2(midRow,:,:) ) );
    corInterp = rot90( squeeze( interp2(midRow,:,:) ) );
    if orientation==0    
      corImg2Save = [ cor2 corInterp ];
      [~, nColsCor] = size(corImg2Save);
      corImg2Save = imresize( corImg2Save, nColsAx/nColsCor );
    else
      corImg2Save = [ cor2; corInterp; ];
      [nRowsCor, ~] = size(corImg2Save);
      corImg2Save = imresize( corImg2Save, nRowsAx/nRowsCor );
    end
    
    sag2 = rot90( squeeze( data2(:,midCol,:) ) );
    sagInterp = rot90( squeeze( interp2(:,midCol,:) ) );
    if orientation==0
      sagImg2Save = [ sag2 sagInterp ];
      [~, nColsSag] = size(sagImg2Save);
      sagImg2Save = imresize( sagImg2Save, nColsAx/nColsSag );
    else
      sagImg2Save = [ sag2; sagInterp; ];
      [nRowsSag, ~] = size(sagImg2Save);
      sagImg2Save = imresize( sagImg2Save, nRowsAx/nRowsSag );
    end

    if orientation==0
      img2Save = [ axialImg2Save; corImg2Save; sagImg2Save; ];
    else
      img2Save = [ axialImg2Save corImg2Save sagImg2Save ];
    end
    img2Save = imresize( img2Save, 3 );

    % % Add text onto image
    %textImg = 1 - text2im( ['Frame ', num2str(i,'%3.3i')] );
    %sImText = size(textImg);
    %img2Save(end-sImText(1)-9:end-10,20:19+sImText(2)) = ...
    %  textImg*max(img2Save(:));

    %imshow( img2Save, [] );
    imwrite( img2Save, ...
      ['./',outDir, '/frame_',num2str(offset+i,'%3.3i'), '.jpg'] );
  end

end

