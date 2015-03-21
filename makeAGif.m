
function makeAGif( varargin )
  % makeAGif( inDir, outFile, delayTime )

  defaultInDir = './alignImgs';
  defaultOutFile = 'orig.gif';
  defaultDelayTime = 0.5;

  p = inputParser;
  p.addOptional('inDir',defaultInDir);
  p.addOptional('outFile',defaultOutFile);
  p.addOptional('delayTime',defaultDelayTime,@isnumeric);
  p.parse( varargin{:} );
  inDir = p.Results.inDir;
  outFile = p.Results.outFile;
  delayTime = p.Results.delayTime;


  imgs = dir( inDir );
  imgs = imgs(3:end);
  if numel(imgs)<1, disp('No images found'); end;

  for i=1:numel(imgs)
    img = imread( [inDir,'/',imgs(i).name] );
    if ndims(img)==2
      img = repmat(img,[1,1,3]);
    end
    [A,map] = rgb2ind(img,256);

    if i == 1;
      imwrite(A,map,outFile,'gif','LoopCount',Inf,'DelayTime',delayTime);
    else
      imwrite(A,map,outFile,'gif','WriteMode','append','DelayTime',delayTime);
    end
  end

end
