
c = .1;
numRows = 128;
numCols = 256;
numPix = numRows*numCols;

D1 = sparse([],[],2*numPix,numPix,numPix);
D2 = sparse([],[],2*numPix,numPix,numPix);

for j = 1:numCols
    for i = 1:(numRows - 1)
        ijIdx = (j-1)*numRows + i;
        ip1jIdx = (j-1)*numRows + i + 1;
        D1(ijIdx,ijIdx) = -1;
        D1(ijIdx,ip1jIdx) = 1;
    end
end

for j = 1:(numCols-1)
    for i = 1:numRows
        ijIdx = (j-1)*numRows + i;
        ijp1Idx = (j + 1 - 1)*numRows + i;
        
        D2(ijIdx,ijIdx) = -1;
        D2(ijIdx,ijp1Idx) = 1;
    end
end

applyD1 = @(u) [u(2:end,:) - u(1:end-1,:) ; zeros(1,size(u,2))];
applyD2 = @(u) [u(:,2:end) - u(:,1:end-1), zeros(size(u,1),1)];
applyD1Trans = @(u) [-u(1,:) ; u(1:end-2,:) - u(2:end-1,:); u(end-1,:)];
applyD2Trans = @(u) [-u(:,1) , u(:,1:end-2) - u(:,2:end-1), u(:,end-1)];

x = randn(numRows,numCols);
D1x = D1*x(:);
D1x = reshape(D1x,numRows,numCols);
D1xCheck = applyD1(x);
diff1 = D1x - D1xCheck;
diff1 = max(abs(diff1(:)));

D2x = D2*x(:);
D2x = reshape(D2x,numRows,numCols);
D2xCheck = applyD2(x);
diff2 = D2x - D2xCheck;
diff2 = max(abs(diff2(:)));

figure; spy(D1); title('D1')
figure; spy(D2); title('D2')

applyM = @(u) u + c*applyD1Trans(applyD1(u)) + c*applyD2Trans(applyD2(u));
M = speye(numPix) + c*D1'*D1 + c*D2'*D2;
figure; spy(M); title('M')

b = randn(numRows,numCols);
tic
xCheck = M\b(:);
timeBackslash = toc;
xCheck = reshape(xCheck,numRows,numCols);

allOnes = ones(numRows,numCols);
eigValsM = dct2(applyM(idct2(allOnes)));
eigValsMInv = eigValsM.^(-1);

% Now we'll solve Mx = b .
tic
x = idct2(eigValsMInv.*dct2(b));
timeDCT = toc
check1 = applyM(x) - b;
check1 = max(abs(check1(:)));

check2 = x - xCheck;
check2 = max(abs(check2(:)));

disp(['timeBackslash: ',num2str(timeBackslash)])
disp(['timeDCT: ',num2str(timeDCT)])
disp(['ratio:', num2str(timeBackslash/timeDCT)])



disp('finished')
