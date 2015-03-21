
function of2Quiver( hFlow, vFlow, varargin )
  %function of2Quiver( hFlow, vFlow, decimate, vThresh, scale )
  %  hFlow is a 2D array representing horizontal velocities in pixels
  %  vFlow is a 2D array representing vertical velocities in pixels
  %  decimate is how often the vectors are displayed
  %  vThresh: no vectors with magnitude less than this amount are displayed
  %  scale: vectors are scaled by this amount for display

  defaultDecimate = 1;
  defaultVThresh = -1;    % Show all velocities by default
  defaultScale = 1;

  p = inputParser;
  p.addRequired('hFlow');
  p.addRequired('vFlow');
  p.addOptional('decimate',defaultDecimate,@isnumeric);
  p.addOptional('vThresh',defaultVThresh,@isnumeric);
  p.addOptional('scale',defaultScale,@isnumeric);
  p.parse( hFlow, vFlow, varargin{:} );
  decimate = p.Results.decimate;
  vThresh = p.Results.vThresh;
  scale = p.Results.scale;

  [M,N] = size(hFlow);
  [Mv,Nv] = size(vFlow);
  if M~=Mv || N~=Nv, error('hFlow and vFlow must be the same size'); end;

  [xs,ys] = meshgrid( 1:N, 1:M );

  xs = xs(1:decimate:end,1:decimate:end);
  ys = ys(1:decimate:end,1:decimate:end);

  if decimate > 1
    hFlowSmooth = filter2(fspecial('average',decimate),hFlow);
    hFlow = hFlowSmooth(1:decimate:end,1:decimate:end);
    vFlowSmooth = filter2(fspecial('average',decimate),vFlow);
    vFlow = vFlowSmooth(1:decimate:end,1:decimate:end);
  end

  vs = sqrt( hFlow.*hFlow + vFlow.*vFlow );
  goodIndxs = find( vs > vThresh );

  if numel(goodIndxs) == 0, disp('of2Quiver: no relevant velocities'); end;

  hFlow = hFlow * scale;
  vFlow = vFlow * scale;
  
  quiver( xs(goodIndxs), ys(goodIndxs), ...
    hFlow(goodIndxs), vFlow(goodIndxs) );

end

