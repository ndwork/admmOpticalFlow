
function mask = loadMask( maskSize )

  mask = zeros( maskSize );

  files = dir( './masks/mask_*' );

  for i=1:numel(files)
    thisMask = imread( ['./masks/', files(i).name ] );
    mask(:,:,i) = (thisMask>0);
  end

end
