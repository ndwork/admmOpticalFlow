
function mask = loadMask( maskSize, datacase )

  mask = zeros( maskSize );

  files = dir( ['./masks_', num2str(datacase), '/mask_*'] );

  for i=1:numel(files)
    thisMask = imread( ['./masks_', num2str(datacase), ...
      '/', files(i).name ] );
    mask(:,:,i) = (thisMask>0);
  end

  %checkMask(mask)
end


function checkMask(mask)

  figure;
  [nRows,nCols,nSlices] = size(mask);
  for i=1:nSlices
    maskedImg = mask(:,:,i) * data(:,:,i);
    imshow( maskedImg, [] );
  end

end
