
function [du,dv] = opticalFlow3D( data1, data2 )

  % optical flow parameters
  %eta = 1d-2;
  eta = 1d-4;

  pyramid1 = makeDataPyramid3( data1, 3 );
  pyramid2 = makeDataPyramid3( data2, 3 );

  sLevel1 = size( pyramid1{end} );
  du = zeros( sLevel1(1), sLevel1(2), sLevel1(3) );
  dv = zeros( sLevel1(1), sLevel1(2), sLevel1(3) );
  dw = zeros( sLevel1(1), sLevel1(2), sLevel1(3) );

  for level=numel(pyramid1):-1:1
    disp(['Working on pyramid level ', num2str(level)]);
    tmp1 = pyramid1{level};
    [nRows nCols nPages] = size( pyramid1{level} );

    % Apply a median filter
    du = medfilt3( du, [5 5 5], 'symmetric' );
    dv = medfilt3( dv, [5 5 5], 'symmetric' );
    dw = medfilt3( dw, [5 5 5], 'symmetric' );

    scale = nRows / size(du,1);
    du = resize(du, [nRows nCols nPages], 'bilinear') * scale;
    dv = resize(dv, [nRows nCols nPages], 'bilinear') * scale;
    dw = resize(dw, [nRows nCols nPages], 'bilinear') * scale;

    tmp2 = ofInterp3D( pyramid2{level}, du, dv, dw );

    [newDu,  newDv, newDw] = ofADMM3D( tmp1, tmp2, eta );

    du = du + newDu;
    dv = dv + newDv;
    dw = dw + newDw;

    %showDiagnostics = 1;
    %if showDiagnostics==1
    %  close all;
    %  interped2 = ofInterp( pyramid2{level}, du, dv );
    %  figure, imshow( imresize( tmp1, ceil(512/nRows), 'nearest' ), [] );
    %  title('img1', 'FontSize', 20);
    %  figure, imshow( imresize( tmp2, ceil(512/nRows), 'nearest' ), [] );
    %  title('img2', 'FontSize', 20);
    %  figure, imshow( imresize( interped2, ceil(512/nRows), 'nearest' ), [] );
    %  title('interped2', 'FontSize', 20);
    %  figure, imshow( imresize( du, ceil(512/nRows), 'nearest' ), [] );
    %  title('du', 'FontSize', 20);
    %  figure, imshow( imresize( dv, ceil(512/nRows), 'nearest' ), [] );
    %  title('dv', 'FontSize', 20);
    %  drawnow;
    %end

  end

end


function [du,dv,dw] = ofADMM3D( data1, data2, eta )
  
  % data1 and data2 must be of type double
  data1 = double( data1 );
  data2 = double( data2 );

  % ADMM Parameters
  rho = 1.0;

  [nRows, nCols, nPages] = size( data1 );
  nPix = nRows*nCols*nPages;

  [Iu1, Iv1, Iw1] = imgDeriv3D( data1 );
  [Iu2, Iv2, Iw2] = imgDeriv3D( data2 );
  Iu = ( Iu1 + Iu2 ) / 2;
  Iv = ( Iv1 + Iv2 ) / 2;
  Iw = ( Iw1 + Iw2 ) / 2;
  It = data2 - data1;

  D1 = sparse(nPix,nPix);
  for k = 1:nPages
    for j = 1:(nCols-1)
      for i = 1:nRows
        ijkIdx = (k-1)*nRows*nPages + (j-1)*nRows + i;
        ijp1kIdx = (k-1)*nRows*nPages + (j+1-1)*nRows + i;
        D1(ijkIdx,ijkIdx) = -1;
        D1(ijkIdx,ijp1kIdx) = 1;
      end
    end
  end

  D2 = sparse(nPix,nPix);
  for k = 1:nPages
    for j = 1:nCols
      for i = 1:(nRows - 1)
        ijkIdx = (k-1)*nRows*nPages + (j-1)*nRows + i;
        ip1jkIdx = (k-1)*nRows*nPages + (j-1)*nRows + i + 1;
        D2(ijkIdx,ijkIdx) = -1;
        D2(ijkIdx,ip1jkIdx) = 1;
      end
    end
  end
  
  D3 = sparse(nPix,nPix);
  for k = 1:(nPages-1)
    for j = 1:nCols
      for i = 1:(nRows-1)
        ijkIdx = (k-1)*nRows*nPages + (j-1)*nRows + i;
        ijkp1Idx = (k-1+1)*nRows*nPages + (j-1)*nRows + i;
        D3(ijkIdx,ijkIdx) = -1;
        D3(ijkIdx,ijkp1Idx) = 1;
      end
    end
  end

  B = speye(nPix) + D1'*D1 + D2'*D2 + D3'*D3;

  % Store an LU decomposition of M for faster solving later
  % NOTE: PMQ = LU, and Py = y(p); Qy = y(qInv);
  [L,U,p,q] = lu(B,'vector');
  p = p';
  q = q';
  pInv(p) = 1:nPix;
  qInv(q) = 1:nPix;
  pInv = pInv';
  qInv = qInv';
  
  % Initializations
  x1 = zeros( nPix, 1 );
  x2 = zeros( nPix, 1 );
  x3 = zeros( nPix, 1 );
  y1 = zeros( nPix, 1 );
  y2 = zeros( nPix, 1 );
  y3 = zeros( nPix, 1 );
  z1_1 = zeros( nPix, 1 );
  z1_2 = zeros( nPix, 1 );
  z1_3 = zeros( nPix, 1 );
  z2_1 = zeros( nPix, 1 );
  z2_2 = zeros( nPix, 1 );
  z2_3 = zeros( nPix, 1 );
  z3_1 = zeros( nPix, 1 );
  z3_2 = zeros( nPix, 1 );
  z3_3 = zeros( nPix, 1 );

  lambda1_1 = zeros( nPix, 1 );
  lambda1_2 = zeros( nPix, 1 );
  lambda1_3 = zeros( nPix, 1 );
  lambda2_1 = zeros( nPix, 1 );
  lambda2_2 = zeros( nPix, 1 );
  lambda2_3 = zeros( nPix, 1 );
  lambda3_1 = zeros( nPix, 1 );
  lambda3_2 = zeros( nPix, 1 );
  lambda3_3 = zeros( nPix, 1 );
  lambda4_1 = zeros( nPix, 1 );
  lambda4_2 = zeros( nPix, 1 );
  lambda4_3 = zeros( nPix, 1 );

  b = -It(:);
  Au = Iu(:);
  Av = Iv(:);
  Aw = Iw(:);

  %cvx_begin
  %  variables cvxDu(nPix,1) cvxDv(nPix,1)
  %  minimize( 0.5*sum( (Au.*cvxDu + Av.*cvxDv - b).^2 ) + ...
  %    eta*norm(D1*cvxDu,1) + eta*norm(D2*cvxDu,1) + ...
  %    eta*norm(D1*cvxDv,1) + eta*norm(D2*cvxDv,1) );
  %cvx_end
  %cvxDu = reshape( cvxDu, [nRows nCols] );
  %cvxDv = reshape( cvxDv, [nRows nCols] );
  
  M11 = Au.*Au + 1;   M12 = Au.*Av;         M13 = Au.*Aw;
  M21 = M12;          M22 = Av.*Av + 1;     M23 = Av.*Aw;
  M31 = M13;          M32 = M23;            M33 = Aw.*Aw+1;

  K1 = M32 - M31./M11.*M12;
  K2 = M22 - M21./M11.*M12;
  K = K1 ./ K2;
  C = M31./M11 - K.*M21./M11;
  
  nIter = 1000;
  %objectives = zeros(nIter,1);
  for i=1:nIter
    if mod(i,50)==0
      disp([ 'Working on iteration ', num2str(i), ' of ', num2str(nIter) ]);
    end

    %objectives(i) = ofObjective( Au, Av, b, x1, x2, eta, D1, D2 );

    % Update x
    arg1 = y1 + D1'*z1_1 + D2'*z1_2 + D3'*z1_3 - lambda1_1/rho - ...
      (D1'*lambda2_1 + D2'*lambda2_2 + D3'*lambda2_3)/rho;
    x1 = U\(L\(arg1(p)));
    x1 = x1(qInv);
    %x1 = M \ arg1;
    arg2 = y2 + D1'*z2_1 + D2'*z2_2 + D3'*z2_3 - lambda1_2/rho - ...
      (D1'*lambda3_1 + D2'*lambda3_2 + D3'*lambda3_3)/rho;
    x2 = U\(L\(arg2(p)));
    x2 = x2(qInv);
    %x2 = M \ arg2;
    arg3 = y3 + D1'*z3_1 + D2'*z3_2 + D3'*z3_3 - lambda1_3/rho - ...
      (D1'*lambda4_1 + D2'*lambda4_2 + D3'*lambda4_3);
    x3 = U\(L\(arg3(p)));
    x3 = x3(qInv);

    % Update y
    nu1 = Au.*b + lambda1_1 + rho*x1;
    nu2 = Av.*b + lambda1_2 + rho*x2;
    nu3 = Aw.*b + lambda1_3 + rho*x3;
    y3 = (nu3 - K.*x2 - C.*x1 ) ./ ...
      (K.*M21./M11.*M13 - K.*M23 - M31./M11.*M13 + M33 );
    y2 = ( nu2 - M21./M11.*nu1 + M21./M11.*M13.*nu3 - M23.*nu3 ) ./ ...
      (M22 - M21./M11.*M12);
    y1 = ( nu1 - M12.*y2 - M13.*y3 ) ./ M11;

    % Update z
    z1_1 = softThresh( D1*x1 + lambda2_1/rho, eta/rho );
    z1_2 = softThresh( D2*x1 + lambda2_2/rho, eta/rho );
    z1_3 = softThresh( D3*x1 + lambda2_3/rho, eta/rho );
    z2_1 = softThresh( D1*x2 + lambda3_1/rho, eta/rho );
    z2_2 = softThresh( D2*x2 + lambda3_2/rho, eta/rho );
    z2_3 = softThresh( D3*x2 + lambda3_3/rho, eta/rho );
    z3_1 = softThresh( D1*x3 + lambda4_1/rho, eta/rho );
    z3_2 = softThresh( D2*x3 + lambda4_2/rho, eta/rho );
    z3_3 = softThresh( D3*x3 + lambda4_3/rho, eta/rho );

    % Update lambdas
    lambda1_1 = lambda1_1 + rho * ( x1 - y1 );
    lambda1_2 = lambda1_2 + rho * ( x2 - y2 );
    lambda1_3 = lambda1_3 + rho * ( x3 - y3 );
    lambda2_1 = lambda2_1 + rho * ( D1*x1 - z1_1 );
    lambda2_2 = lambda2_2 + rho * ( D2*x1 - z1_2 );
    lambda2_3 = lambda2_3 + rho * ( D3*x1 - z1_3 );
    lambda3_1 = lambda3_1 + rho * ( D1*x2 - z2_1 );
    lambda3_2 = lambda3_2 + rho * ( D2*x2 - z2_2 );
    lambda3_3 = lambda3_3 + rho * ( D3*x2 - z2_3 );
    lambda4_1 = lambda4_1 + rho * ( D1*x3 - z3_1 );
    lambda4_2 = lambda4_2 + rho * ( D2*x3 - z3_2 );
    lambda4_3 = lambda4_3 + rho * ( D3*x3 - z3_3 );

  end

  du = reshape( x1, [nRows nCols nPages] );
  dv = reshape( x2, [nRows nCols nPages] );
  dw = reshape( x3, [nRows nCols nPages] );

  % for diagnostics
save( 'state.mat' );
  showDiagnostics = 0;
  if showDiagnostics==1
    close all;
    admmOptVal = ofObjective( Au, Av, b, x1, x2, eta, D1, D2 );
    disp(['ADMM Optimal Value: ', num2str(admmOptVal) ] );
    ofRes = ofResidual3D( Iu, Iv, It, du, dv );
    figure, imshow( imresize( ofRes, ceil(512/nRows), 'nearest' ), [] );
    title('OF Residual', 'FontSize', 20);
    %load( 'star.mat' );
    %objStar = ofObjective( Au, Av, b, x1Star, x2Star, eta, D1, D2 );
    %figure, plot( ( objectives - objStar ) / objStar );
    %title('relative error vs iteration');
    drawnow;
  end
end


function [dx, dy, dz] = imgDeriv3D( data )
  [M, N, K] = size( data );
  dx = zeros(M,N,K);
  dy = zeros(M,N,K);
  dz = zeros(M,N,K);

  %dx(:,1:nCols-1) = img1(:,2:nCols) - img1(:,1:nCols-1);
  %dy(1:nRows-1,:) = img1(2:nRows,:) - img1(1:nRows-1,:);
  
  dx(:,2:N-1,:) = ( data(:,3:N,:) - data(:,1:N-2,:) ) / 2;
  dx(:,1,:) = data(:,2,:) - data(:,1,:);
  dx(:,N,:) = data(:,N,:) - data(:,N-1,:);
  
  dy(2:M-1,:,:) = ( data(3:M,:,:) - data(1:M-2,:,:) ) / 2;
  dy(1,:,:) = data(2,:,:) - data(1,:,:);
  dy(M,:,:) = data(M,:,:) - data(M-1,:,:);
  
  dz(:,:,2:K-1) = ( data(:,:,3:K) - data(:,:,1:K-2) ) / 2;
  dz(:,:,1) = data(:,:,2) - data(:,:,1);
  dz(:,:,K) = data(:,:,K) - data(:,:,K-1);
end


function ofRes = ofResidual3D( Iu, Iv, Iw, It, du, dv, dw )
  ofRes = Iu.*du + Iv.*dv + It + Iw.*dw;
end


function out = softThresh( in, thresh )
  out = sign(in) .* max( ( abs(in) - thresh ), 0 );
end


function out = ofObjective( Iu, Iv, b, du, dv, eta, D1, D2 )
  Agam = Iu.*du + Iv.*dv;
  Agamb = Agam - b;
  out = 0.5*sum(Agamb.*Agamb) + ...
    eta*norm(D1*du,1) + eta*norm(D2*du,1) + ...
    eta*norm(D1*dv,1) + eta*norm(D2*dv,1);
end


function pyramid = makeDataPyramid3( data, nLevels, spacing )
  if nargin < 2
    nLevels = 4;
  end
  if nargin < 3
    spacing = 2;
  end

  smooth_sigma = spacing/2;
  f = gaussFilt3( smooth_sigma, 2*round(1.5*smooth_sigma)+1 );
  ratio = 1. / spacing;

  pyramid = cell(nLevels,1);
  tmp = data;
  for m = 1:nLevels
    pyramid{m} = tmp;
    sTmp = size( tmp );
    newSize = floor( sTmp .* ratio );
    tmp = imfilter(tmp, f, 'corr', 'symmetric', 'same');  % Gauss filter
    tmp = resize(tmp, newSize, 'bilinear');  % Downsampling
  end
end


function h = gaussFilt3( sig, hsize )
  siz   = (hsize-1)/2;
  if ndims(siz) ~= 3 siz = ones(1,3)*siz; end;
  if ndims(sig) ~= 3 sig = ones(1,3)*sig; end;
  [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
  h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
  h = h/sum(h(:));
end

