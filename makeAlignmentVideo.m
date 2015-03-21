
function makeAlignmentVideo
  close all; clear;

  maxNData = 100;

  %datacase = 1;
  datacase = 9;
  %datacase = 11;
  
  outDir = ['alignImgs_',num2str(datacase)];
  
  [x,eta,rho,offset] = loadData( datacase );
  x = abs(x);
  [nRows,nCols,nPages,nData] = size(x);
  midRow = ceil(nRows/2);
  midCol = ceil(nCols/2);
  midSlice = ceil(nPages/2);
  nData = min([offset+nData,maxNData,nData]);


  mkdir( outDir );
  data1 = squeeze( x(:,:,:,offset+1) );
  parfor i=1:nData-offset 
    disp(['Working on image ', num2str(i)]);
    data2 = squeeze( x(:,:,:,offset+i) );

    maskedIndxs = find( data2 ~= 0 );
    scales = data1(maskedIndxs) ./ data2(maskedIndxs);
    data2scaled = data2 .* median( scales );

    [du,dv,dw] = opticalFlow3D( data1, data2scaled, eta, rho );
    interp2 = ofInterp3D( data2, du, dv, dw );

    axial2 = data2(:,:,midSlice);
    axialInterp = interp2(:,:,midSlice);
    axialImg2Save = [ axial2 axialInterp ];
    [~, nColsAx] = size(axialImg2Save);

    sag2 = rot90( squeeze( data2(:,midCol,:) ) );
    sagInterp = rot90( squeeze( interp2(:,midCol,:) ) );
    sagImg2Save = [ sag2 sagInterp ];
    [~, nColsSag] = size(sagImg2Save);
    sagImg2Save = imresize( sagImg2Save, nColsAx/nColsSag );

    cor2 = rot90( squeeze( data2(midRow,:,:) ) );
    corInterp = rot90( squeeze( interp2(midRow,:,:) ) );
    corImg2Save = [ cor2 corInterp ];
    [~, nColsCor] = size(corImg2Save);
    corImg2Save = imresize( corImg2Save, nColsAx/nColsCor );

    img2Save = [ axialImg2Save; sagImg2Save; corImg2Save ];
    img2Save = imresize( img2Save, 3 );

    % Add text onto image
    textImg = 1 - text2im( ['Frame ', num2str(i,'%3.3i')] );
    sImText = size(textImg);
    img2Save(end-sImText(1)-9:end-10,20:19+sImText(2)) = ...
      textImg*max(img2Save(:));

    %imshow( img2Save, [] );
    imwrite( img2Save, ...
      ['./',outDir, '/frame_',num2str(offset+i,'%3.3i'), '.jpg'] );
  end

end

