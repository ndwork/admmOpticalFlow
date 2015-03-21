
function saveAllSlices( data, varargin )
  % saveAllSlices( data, dim, resize, outDir );
  % data is 3D MRI data
  % dim is the dimension to use as the slice index (1,2,3)
  % outDir is a character string to call the output directory

  defaultDimension = 1;
  defaultOutDir = 'slices';
  defaultResize = 1;

  p = inputParser;
  p.addRequired( 'data' );
  p.addOptional( 'dim', defaultDimension, @isnumeric );
  p.addOptional( 'resize', defaultResize, @isnumeric );
  p.addOptional( 'outDir', defaultOutDir );
  p.parse( data, varargin{:} );
  dim = p.Results.dim;
  p = inputParser;
  p.addRequired( 'data' );
  p.addOptional( 'dim', dim, @isnumeric );
  defaultOutDir = ['slices',num2str(dim)];
  p.addOptional( 'resize', defaultResize, @isnumeric );
  p.addOptional( 'outDir', defaultOutDir );
  p.parse( data, varargin{:} );
  outDir = p.Results.outDir;
  resize = p.Results.resize;

  [nRows,nCols,nPages] = size(data);
  mkdir(['./',outDir]);

  switch dim

    case 1
      for i=1:nRows
       thisImg = imresize( squeeze( data(i,:,:) ), resize );
       thisImg = thisImg / max( thisImg(:) );
       thisImgFilename = ['./',outDir,'/slice_', num2str(i,'%3.3i'), '.jpg'];
       imwrite( thisImg, thisImgFilename, 'jpg' );
      end

    case 2
      for i=1:nCols
       thisImg = imresize( squeeze( data(:,i,:) ), resize );
       thisImg = thisImg / max( thisImg(:) );
       thisImgFilename = ['./',outDir,'/slice_', num2str(i,'%3.3i'), '.jpg'];
       imwrite( thisImg, thisImgFilename, 'jpg' );
      end

    case 3
      for i=1:nPages
       thisImg = imresize( squeeze( data(:,:,i) ), resize );
       thisImg = thisImg / max( thisImg(:) );
       thisImgFilename = ['./',outDir,'/slice_', num2str(i,'%3.3i'), '.jpg'];
       imwrite( thisImg, thisImgFilename, 'jpg' );
      end

  end

end

