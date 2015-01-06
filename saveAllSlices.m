
function saveAllSlices( x )

  [nRows,nCols,nSlices,nData] = size(x);

  midDataIndx = ceil(nData/2);
  mkdir('./slices');
  for i=1:nSlices
   thisImg = squeeze( x(:,:,i,midDataIndx) );
   thisImg = thisImg / max( thisImg(:) );
   thisImgFilename = ['./slices/slice_', num2str(i,'%3.3i'), '.jpg'];
   imwrite( thisImg, thisImgFilename, 'jpg' );
  end

end

