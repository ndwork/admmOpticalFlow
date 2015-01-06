
function x = loadData( datacase )
  switch datacase
    case 1
      load '/Volumes/Seagate2TB/Data/mriCardiacNavs/esp3dn25148.mat';
      % The max dimensions of x are (128,128,64,508);
      x = x / 500;
    case 2
      load '/Volumes/Seagate2TB/Data/cardiacRecons/res2mm/img318143.mat';
      % The max dimensions of m are (120,120,80,20);  
      x = m;  clearvars m;
    case 3
      % Not very much of a heartbeat
      load '/Volumes/Seagate2TB/Data/cardiacRecons/res1.2mm/img32013.mat';
      % The max dimensions of m are (256,256,112,6);
      x = m;  clearvars m;
    case 4
      x = loadSerenaData(1);
      x = x / max( x(:));
    case 5
      x = loadSerenaData(2);
      x = x / max( x(:));
    case 6
      x = loadSerenaData(3);
      x = x / max( x(:));
    case 7
      x = loadSerenaData(4);
      x = x / max( x(:));
  end
end


function imgs = loadSerenaData( fileIndx )
  imgFiles = cell(4,1);
  imgFiles{1} = '/Volumes/Seagate2TB/Data/cardiacImgs/P02560_img.mat';
  imgFiles{2} = '/Volumes/Seagate2TB/Data/cardiacImgs/P27648_img.mat';
  imgFiles{3} = '/Volumes/Seagate2TB/Data/cardiacImgs/P62976_img.mat';
  imgFiles{4} = '/Volumes/Seagate2TB/Data/cardiacImgs/P712133_img.mat';

  img = 0;
  load( imgFiles{fileIndx} );

  sImgFiles = size( img );
  nCoils = sImgFiles(4);
  nPhases = sImgFiles(5);
  imgs = zeros([ sImgFiles(1:3), nPhases ]);

  for phase = 1:nPhases
    thisImg = abs(img(:,:,:,1,phase)).^2;
    for coil=2:nCoils
      % Combine images from each coil
      thisImg = thisImg + abs(img(:,:,:,coil,phase)).^2;
    end
    imgs(:,:,:,phase) = sqrt( thisImg );
  end

end

