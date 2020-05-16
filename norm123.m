%gplData = readtable('GPL1261-56135.xlsx'); 
%  gseData.Data.RowNames
%  gplProbesetIDs = gplData.Data(:, strcmp(gplData.ColumnNames, 'ID'));
%  geneSymbols = gplData.Data(:, strcmp(gplData.ColumnNames, 'Gene Symbol'));
%  gseData.Data = rownames(gseData.Data, ':', geneSymbols);
% gseData = geoseriesread('GSE23006_series_matrix.txt');
% mairplot(DataX, DataY)
% scatter(ak,mk,1,'b')
%------------------------------------------------------------------------------



%Qx = quan(gseData1);
%Rx = cloess(gseData1);
%Lx = linearnorm(gseData1);
%mairplot(Lx(:,4), Lx(:,6))
var = zeros(4,48);
bias = zeros(4,48);
for i = 1:3:46
    [var(1,i),bias(1,i)] = mse1(Nx,gseData1,i,i+1);
    [var(1,i+1),bias(1,i+1)] = mse1(Nx,gseData1,i,i+2);
    [var(1,i+2),bias(1,i+2)] = mse1(Nx,gseData1,i+1,i+2);
end
for i = 1:3:46
    [var(2,i),bias(2,i)] = mse1(Lx,gseData1,i,i+1);
    [var(2,i+1),bias(2,i+1)] = mse1(Lx,gseData1,i,i+2);
    [var(2,i+2),bias(2,i+2)] = mse1(Lx,gseData1,i+1,i+2);
end
for i = 1:3:46
    [var(3,i),bias(3,i)] = mse1(Qx,gseData1,i,i+1);
    [var(3,i+1),bias(3,i+1)] = mse1(Qx,gseData1,i,i+2);
    [var(3,i+2),bias(3,i+2)] = mse1(Qx,gseData1,i+1,i+2);
end
for i = 1:3:46
    [var(4,i),bias(4,i)] = mse1(Rx,gseData1,i,i+1);
    [var(4,i+1),bias(4,i+1)] = mse1(Rx,gseData1,i,i+2);
    [var(4,i+2),bias(4,i+2)] = mse1(Rx,gseData1,i+1,i+2);
end
xlswrite('var', var);
xlswrite('bias', bias);
%---------------------------------------
function Rx = cloess(x)
m = zeros(1,45101)';
a = zeros(1,45101)';
for j = 1:3:46
    %1-2
    for k = 1:45101
     m(k) = log2(x(k,j)/x(k,j+1));
     a(k) = 0.5*log2(x(k,j)*x(k,j+1));
    end
    msmooth = malowess(a,m);
    mnorm = m - msmooth;
    for k = 1:45101
        x(k,j) = 2^(a(k)+0.5*mnorm(k));
        x(k,j+1) = 2^(a(k)-0.5*mnorm(k));
    end
    %2-3
    for k = 1:45101
     m(k) = log2(x(k,j+1)/x(k,j+2));
     a(k) = 0.5*log2(x(k,j+1)*x(k,j+2));
    end
    msmooth = malowess(a,m);
    mnorm = m - msmooth;
    for k = 1:45101
        x(k,j+1) = 2^(a(k)+0.5*mnorm(k));
        x(k,j+2) = 2^(a(k)-0.5*mnorm(k));
    end
    %1-3
    for k = 1:45101
     m(k) = log2(x(k,j)/x(k,j+2));
     a(k) = 0.5*log2(x(k,j)*x(k,j+2));
    end
    msmooth = malowess(a,m);
    mnorm = m - msmooth;
    for k = 1:45101
        x(k,j) = 2^(a(k)+0.5*mnorm(k));
        x(k,j+2) = 2^(a(k)-0.5*mnorm(k));
    end
end
Rx = x;
end


function Qx = quan(x)
Qx = zeros(45101,1);
for i = 1:3:46;
    Q = quantilenorm([x(:,i),x(:,i+1),x(:,i+2)]);
    Qx = [Qx,Q];
end
Qx = Qx(:,2:49);
end

function Lx = linearnorm(x)
for i = 1:3:46
    Lx(:,i) = x(:,i);
    [Lx(:,i+1)]=x(:,i+1)/( mean(x(:,i+1)) / mean(x(:,i)) );
    [Lx(:,i+2)]=x(:,i+2)/( mean(x(:,i+2)) / mean(x(:,i)) );
end
end

function [var,bias] = mse1(esm,truem,j,i)
m = zeros(45101,2);
for k = 1:45101
     m(k,1) = log2(esm(k,j)/esm(k,i));
end
for k = 1:45101
    m(k,2) = log2(truem(k,j)/truem(k,i));
end

var = mean(  (m(:,1) - mean(m(:,1))).^2 );
bias = mean(  (m(:,2) - mean(m(:,1))).^2 );
bias = sqrt(bias);
end
