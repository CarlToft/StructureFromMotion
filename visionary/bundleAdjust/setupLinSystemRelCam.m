function [res, J] = setupLinSystemRelCam(P, relPoses, nPoints, res, J)
%w is square root of the inverse of the covariance matrix for 
%[x, y, z, a, b, c]^T



N = length(relPoses);
m = 3+9;
resAdd = zeros(m*N,1);

for i = 1:N
    I = (relPoses(i).M*[P{relPoses(i).i}; 0 0 0 1]) / ...
        [P{relPoses(i).j}; 0 0 0 1];
    resAdd((i-1)*m+(1:3)) = I(1:3,4).*reshape(relPoses(i).w(1:3),3,1);
    eAngle = (I(1:3,1:3)-eye(3));
    resAdd((i-1)*m+(4:12)) = eAngle(:)*mean(relPoses(i).w(4:end));
end
res = [res; resAdd];

if nargout >= 2
    %Need to output Jacobian as well
    nVars = size(J,2);
    nCams = numel(P);
    assert(nVars == nPoints*3 + nCams*6);

    %Rotation manifold bases
    Bs(:,:,1) = [0 1 0; -1 0 0; 0 0 0];
    Bs(:,:,2) = [0 0 1; 0 0 0; -1 0 0];
    Bs(:,:,3) = [0 0 0; 0 0 1; 0 -1 0];
    
    
    %Analytical Jacobian
    col = [];
    row = [];
    data = [];
    
    for k = 1:N
        %Shorter variable names
        i = relPoses(k).i;
        j = relPoses(k).j;
        Ri = P{i}(1:3,1:3);
        Rj = P{j}(1:3,1:3);
        Rm = relPoses(k).M(1:3,1:3);
        tj = P{j}(1:3,4);
        wR = mean(relPoses(k).w(4:end));
        wt = reshape(relPoses(k).w(1:3),3,1);
        
        %dRotation/d(a,b,c)
        dRRi = zeros(3,3,3);
        dRRj = zeros(3,3,3);
        for l = 1:3
            dRRi(:,:,l) = Rm*Bs(:,:,l)*Ri*Rj';
            dRRj(:,:,l) = Rm*Ri*(Bs(:,:,l)*Rj)';
        end
        dRRi = wR*dRRi;
        dRRj = wR*dRRj;
        %rows and cols
        rowsRR = repmat((k-1)*m+(4:12)',6,1);
        colsi = 3*nPoints+(i-1)*6+(1:3);
        colsj = 3*nPoints+(j-1)*6+(1:3);
        colsRR = [repmat(colsi,9,1), repmat(colsj,9,1)];
        
        %dTranslation/d(a,b,c)
        dtRi = zeros(3,3);
        dtRj = zeros(3,3);
        for l = 1:3
            dtRi(:,l) = -Rm*Bs(:,:,l)*Ri*Rj'*tj;
            dtRj(:,l) = -Rm*Ri*(Bs(:,:,l)*Rj)'*tj;
        end
        dtRi = wt.*dtRi;
        dtRj = wt.*dtRj;
        %rows and cols
        rowstR = repmat((k-1)*m+(1:3)',6,1);
        colsi = 3*nPoints+(i-1)*6+(1:3);
        colsj = 3*nPoints+(j-1)*6+(1:3);
        colstR = [repmat(colsi,3,1), repmat(colsj,3,1)];
        
        %dTranslation/d(t1,t2,t3)
        dtti = wt.*Rm;
        dttj = wt.*(-Rm*Ri*Rj');
        %rows and cols
        rowstt = repmat((k-1)*m+(1:3)',6,1);
        colsi = 3*nPoints+(i-1)*6+(4:6);
        colsj = 3*nPoints+(j-1)*6+(4:6);
        colstt = [repmat(colsi,3,1), repmat(colsj,3,1)];

        %Append to row/col/data
        row = [row; rowsRR; rowstR; rowstt];
        col = [col; colsRR(:); colstR(:); colstt(:)];
        data = [data; dRRi(:); dRRj(:); dtRi(:); dtRj(:); dtti(:); dttj(:)];
        
    end
    ix = abs(data) > 1e-10;
    J = [J; sparse(row(ix),col(ix),data(ix),N*m, nVars)];

end


end