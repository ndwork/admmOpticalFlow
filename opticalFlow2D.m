
function [du,dv] = opticalFlow2D( img1, img2 )

  % optical flow parameters
  %eta = 1d-2;
  eta = 1d-4;

  pyramid1 = makeImagePyramid( img1 );
  pyramid2 = makeImagePyramid( img2 );

  sLevel1 = size( pyramid1{end} );
  du = zeros( sLevel1(1), sLevel1(2) );
  dv = zeros( sLevel1(1), sLevel1(2) );

  for level=numel(pyramid1):-1:1
    disp(['Working on pyramid level ', num2str(level)]);
    tmp1 = pyramid1{level};
    [nRows, nCols] = size( pyramid1{level} );

    du = medfilt2( du, [3 3] );
    dv = medfilt2( dv, [3 3] );

    scale = nRows / size(du,1);
    du = imresize(du, [nRows nCols], 'bilinear') * scale;
    dv = imresize(dv, [nRows nCols], 'bilinear') * scale;

    tmp2 = ofInterp( pyramid2{level}, du, dv );

    [newDu,  newDv] = ofADMM( tmp1, tmp2, eta );

    du = du + newDu;
    dv = dv + newDv;

    showDiagnostics = 1;
    if showDiagnostics==1
      close all;
      interped2 = ofInterp( pyramid2{level}, du, dv );
      figure, imshow( imresize( tmp1, ceil(512/nRows), 'nearest' ), [] );
      title('img1', 'FontSize', 20);
      figure, imshow( imresize( tmp2, ceil(512/nRows), 'nearest' ), [] );
      title('img2', 'FontSize', 20);
      figure, imshow( imresize( interped2, ceil(512/nRows), 'nearest' ), [] );
      title('interped2', 'FontSize', 20);
      figure, imshow( imresize( du, ceil(512/nRows), 'nearest' ), [] );
      title('du', 'FontSize', 20);
      figure, imshow( imresize( dv, ceil(512/nRows), 'nearest' ), [] );
      title('dv', 'FontSize', 20);
      drawnow;
    end

  end

end


function [du,dv] = ofADMM( img1, img2, eta )

  % ADMM Parameters
  rho = 1.0;

  [nRows, nCols] = size( img1 );
  nPix = nRows*nCols;

  [Iu1, Iv1] = imgDeriv( img1 );
  [Iu2, Iv2] = imgDeriv( img2 );
  Iu = ( Iu1 + Iu2 ) / 2;
  Iv = ( Iv1 + Iv2 ) / 2;
  It = img2 - img1;

  D1 = makeD1( nRows, nCols );
  D2 = makeD2( nRows, nCols );

  B = speye(nPix) + D1'*D1 + D2'*D2;

  % Store an LU decomposition of M for faster solving later
  [L,U,p,q] = lu(B,'vector'); % NOTE: PMQ = LU, and Py = y(p); Qy = y(qInv);
  p = p';
  q = q';
  pInv(p) = 1:nPix;
  qInv(q) = 1:nPix;
  pInv = pInv';
  qInv = qInv';
  
  % Store for update X with DCT
  applyM = @(u) u + D1'*D1*u + D2'*D2*u;
  allOnes = ones(nRows,nCols);
  %idctAllOnes = idct2(allOnes);
  idctAllOnes = mirt_idctn(allOnes);
  idctAllOnes = idctAllOnes(:);
  MIdctAllOnes = applyM(idctAllOnes);
  %eigValsM = dct2( reshape( MIdctAllOnes, [nRows nCols] ) );
  eigValsM = mirt_dctn( reshape( MIdctAllOnes, [nRows nCols] ) );
  eigValsMInv = eigValsM.^(-1);


% Now we'll solve Mx = b .
b = randn(nRows,nCols);
xDCT = updateXwDCT( eigValsMInv, nRows, nCols, b );
xLU = updateXwLU(L, U, p, qInv, b );
diff = sum( abs( xDCT(:) - xLU(:) ) );
disp(['Diff is: ', num2str(diff)]);
  
  
  % Initializations
  x1 = zeros( nPix, 1 );
  x2 = zeros( nPix, 1 );
  y1 = zeros( nPix, 1 );
  y2 = zeros( nPix, 1 );
  z1_1 = zeros( nPix, 1 );
  z1_2 = zeros( nPix, 1 );
  z2_1 = zeros( nPix, 1 );
  z2_2 = zeros( nPix, 1 );

  lambda1_1 = zeros( nPix, 1 );
  lambda1_2 = zeros( nPix, 1 );
  lambda2_1 = zeros( nPix, 1 );
  lambda2_2 = zeros( nPix, 1 );
  lambda3_1 = zeros( nPix, 1 );
  lambda3_2 = zeros( nPix, 1 );

  b = -It(:);
  Au = Iu(:);
  Av = Iv(:);

  %cvx_begin
  %  variables cvxDu(nPix,1) cvxDv(nPix,1)
  %  minimize( 0.5*sum( (Au.*cvxDu + Av.*cvxDv - b).^2 ) + ...
  %    eta*norm(D1*cvxDu,1) + eta*norm(D2*cvxDu,1) + ...
  %    eta*norm(D1*cvxDv,1) + eta*norm(D2*cvxDv,1) );
  %cvx_end
  %cvxDu = reshape( cvxDu, [nRows nCols] );
  %cvxDv = reshape( cvxDv, [nRows nCols] );
  
  M11 = Au.*Au + 1;
  M12 = Au.*Av;
  M21 = Av.*Au;
  M22 = Av.*Av + 1;

  nIter = 1000;
  %objectives = zeros(nIter,1);
  for i=1:nIter

    %objectives(i) = ofObjective( Au, Av, b, x1, x2, eta, D1, D2 );

    arg1 = y1 + D1'*z1_1 + D2'*z1_2 - lambda1_1/rho - ...
      (D1'*lambda2_1 + D2'*lambda2_2)/rho;
    %x1 = B \ arg1;
    %x1 = updateXwLU( L, U, p, qInv, arg1);
    x1 = updateXwDCT( eigValsMInv, nRows, nCols, arg1 );

    arg2 = y2 + D1'*z2_1 + D2'*z2_2 - lambda1_2/rho - ...
      (D1'*lambda3_1 + D2'*lambda3_2)/rho;
    %x2 = B \ arg2;
    %x2 = updateXwLU( L, U, p, qInv, arg2 );
    x2 = updateXwDCT( eigValsMInv, nRows, nCols, arg2 );


    % Update y
    nu1 = Au.*b + lambda1_1 + rho*x1;
    nu2 = Av.*b + lambda1_2 + rho*x2;
    y2 = (nu2 - M21./M11.*nu1) ./ ( M22 - M21./M11.*M12 );
    y1 = 1./M11.*nu1 - 1./M11.*M12.*y2;

    % Update z
    z1_1 = softThresh( D1*x1 + lambda2_1/rho, eta/rho );
    z1_2 = softThresh( D2*x1 + lambda2_2/rho, eta/rho );
    z2_1 = softThresh( D1*x2 + lambda3_1/rho, eta/rho );
    z2_2 = softThresh( D2*x2 + lambda3_2/rho, eta/rho );

    % Update lambdas
    lambda1_1 = lambda1_1 + rho * ( x1 - y1 );
    lambda1_2 = lambda1_2 + rho * ( x2 - y2 );
    lambda2_1 = lambda2_1 + rho * ( D1*x1 - z1_1 );
    lambda2_2 = lambda2_2 + rho * ( D2*x1 - z1_2 );
    lambda3_1 = lambda3_1 + rho * ( D1*x2 - z2_1 );
    lambda3_2 = lambda3_2 + rho * ( D2*x2 - z2_2 );

  end

  du = reshape( x1, [nRows nCols] );
  dv = reshape( x2, [nRows nCols] );

  % for diagnostics
  showDiagnostics = 0;
  if showDiagnostics==1
    close all;
    admmOptVal = ofObjective( Au, Av, b, x1, x2, eta, D1, D2 );
    disp(['ADMM Optimal Value: ', num2str(admmOptVal) ] );
    ofRes = ofResidual( Iu, Iv, It, du, dv );
    figure, imshow( imresize( ofRes, ceil(512/nRows), 'nearest' ), [] );
    title('OF Residual', 'FontSize', 20);
    %load( 'star.mat' );
    %objStar = ofObjective( Au, Av, b, x1Star, x2Star, eta, D1, D2 );
    %figure, plot( ( objectives - objStar ) / objStar );
    %title('relative error vs iteration');
    drawnow;
  end
end


function D1 = makeD1( nRows, nCols )
  %D1 = sparse(nPix,nPix);
  %for j = 1:(nCols-1)
  %  for i = 1:nRows
  %    ijIdx = (j-1)*nRows + i;
  %    ijp1Idx = (j + 1 - 1)*nRows + i;
  %    D1(ijIdx,ijIdx) = -1;
  %    D1(ijIdx,ijp1Idx) = 1;
  %  end
  %end
  nPix = nRows * nCols;
  nData = 2* (nRows * (nCols-1));
  rows = zeros(1,nData);
  cols = zeros(1,nData);
  values = zeros(1,nData);
  idx = 1;
  for j = 1:(nCols-1)
    for i = 1:nRows
      tmp = (j-1)*nRows + i;
      rows(idx) = tmp;
      cols(idx) = tmp;
      rows(idx+1) = tmp;
      cols(idx+1) = (j+1-1)*nRows + i;
      values(idx:idx+1) = [-1 1];
      idx = idx + 2;
    end
  end
  D1 = sparse(rows,cols,values,nPix,nPix);
end

function D2 = makeD2( nRows, nCols )
  %D2 = sparse(nPix,nPix);
  %for j = 1:nCols
  %  for i = 1:(nRows - 1)
  %    ijIdx = (j-1)*nRows + i;
  %    ip1jIdx = (j-1)*nRows + i + 1;
  %    D2(ijIdx,ijIdx) = -1;
  %    D2(ijIdx,ip1jIdx) = 1;
  %  end
  %end
  nPix = nRows * nCols;
  nData = 2 * ( nCols * (nRows-1) );
  rows = zeros(1,nData);
  cols = zeros(1,nData);
  values = zeros(1,nData);
  idx = 1;
  for j = 1:nCols
    for i = 1:(nRows - 1)
      tmp = (j-1)*nRows + i;
      rows(idx) = tmp;
      cols(idx) = tmp;
      rows(idx+1) = tmp;
      cols(idx+1) = (j-1)*nRows + i + 1;
      values(idx:idx+1) = [-1 1];
      idx = idx + 2;
    end
  end
  D2 = sparse(rows,cols,values,nPix,nPix);
end

function x = updateXwDCT( eigValsMInv, nRows, nCols, b_in )
  b = reshape( b_in, [nRows nCols] );
  %x = idct2( eigValsMInv .* dct2(b) );
  x = mirt_idctn( eigValsMInv .* mirt_dctn(b) );
  x = x(:);
  %check = applyM(x) - b;
  %check = max(abs(check(:)));
end

function x = updateXwLU(L, U, p, qInv, arg );
  % % The inputs are defined as follows:
  % % Store an LU decomposition of B for faster solving later
  %[L,U,p,q] = lu(B,'vector'); % NOTE: PMQ = LU, and Py = y(p); Qy = y(qInv);
  %p = p';
  %q = q';
  %pInv(p) = 1:nPix;
  %qInv(q) = 1:nPix;
  %pInv = pInv';
  %qInv = qInv';

  x = U\(L\(arg(p)));
  x = x(qInv);
end

function [dx, dy] = imgDeriv( img )
  [M N] = size( img );
  dx = zeros(M,N);
  dy = zeros(M,N);

  %dx(:,1:nCols-1) = img1(:,2:nCols) - img1(:,1:nCols-1);
  %dy(1:nRows-1,:) = img1(2:nRows,:) - img1(1:nRows-1,:);
  
  dx(:,2:N-1) = ( img(:,3:N) - img(:,1:N-2) ) / 2;
  dx(:,1) = img(:,2) - img(:,1);
  dx(:,N) = img(:,N) - img(:,N-1);
  
  dy(2:M-1,:) = ( img(3:M,:) - img(1:M-2,:) ) / 2;
  dy(1,:) = img(2,:) - img(1,:);
  dy(M,:) = img(M,:) - img(M-1,:);
end


function ofRes = ofResidual( Iu, Iv, It, du, dv )
  ofRes = Iu.*du + Iv.*dv + It;
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


function pyramid = makeImagePyramid( img, nLevels, spacing )

  if nargin < 2
    nLevels = 4;
  end
  if nargin < 3
    spacing = 2;
  end

  smooth_sigma = spacing/2;
  f = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);
  ratio = 1. / spacing;

  pyramid = cell(nLevels,1);
  tmp = img;
  for m = 1:nLevels
    pyramid{m} = tmp;
    tmp = imfilter(tmp, f, 'corr', 'symmetric', 'same');  % Gauss filter
    tmp = imresize(tmp, ratio, 'bilinear');  % Downsampling
  end
end


