Nx = zeros(45101,48);
for i = 1:3:46
    Nx(:,i) = gseData1(:,i);
    [Nx(:,i+1)]=nonlinearnorm(gseData1,i,i+1);
    [Nx(:,i+2)]=nonlinearnorm(gseData1,i,i+2);
end

function [finalnormgene] = nonlinearnorm(gseData1,col1,col2)
exprvalue = zeros(45101,2);
exprvalue1 = zeros(45101,2);
genename = zeros(45101,1);
genename1 = zeros(45101,1);
% #1-2
[~,idx1] = sort(gseData1(:,col1),'descend');
[~,idx1] = sort(idx1);
[~,idx2] = sort(gseData1(:,col2),'descend');
[~,idx2] = sort(idx2);
k = 1;
for i = 1:45101   
    if (idx1(i,1)-idx2(i,1))/45101<0.001
        genename(k,1) = i;
        exprvalue(k,1) = gseData1(i,col1);
        exprvalue(k,2) = gseData1(i,col2);
        k = k+1;
    end
end

while 1

        [~,idx1] = sort(exprvalue(:,1),'descend');
        [~,idx1] = sort(idx1);
        [~,idx2] = sort(exprvalue(:,2),'descend');
        [~,idx2] = sort(idx2);
        
        i = 1;
        t = 1;
        while genename(i,1) ~= 0
            if (idx1(i,1)-idx2(i,1))/nnz(genename(:,1)) < 0.001
                genename1(t,1) = genename(i,1);
                exprvalue1(t,1) = exprvalue(i,1);
                exprvalue1(t,2) = exprvalue(i,2);
                t = t+1;
            end
            i = i+1;
        end
        if nnz(genename(:,1)) == nnz(genename1(:,1))
            break
        else 
            exprvalue = exprvalue1;
            genename = genename1;
            exprvalue1 = zeros(45101,2);
            genename1 = zeros(45101,1);
        end
end
%----------------------------------------------
[val1,idx1] = sort(gseData1(:,col1),'descend');
[val2,idx2] = sort(gseData1(:,col2),'descend');

pace1 = zeros(nnz(genename),1);
finalnormgene = zeros(45101,1);
for i = 1:nnz(genename)
    for j = 1:45101
        if genename(i,1) == idx1(j,1)
            pace1(i) = j;
            break
        end
    end
end
[sortpace1,~]=sort(pace1,'descend');

pace2 = zeros(nnz(genename),1);
for i = 1:nnz(genename)
    for j = 1:45101
        if genename(i,1) == idx2(j,1)
            pace2(i) = j;
            break
        end
    end
end
[sortpace2,~]=sort(pace2,'descend');

normgene = zeros(45101,1);
for k = 1:sortpace2(1,1)
    normgene(k,1) = val1(k,1);
end
for k = sortpace2(nnz(genename),1):45101
    normgene(k,1) = val1(k,1);
end
for k = 2:nnz(genename)
    frac = mean(val2( (sortpace2(k-1,1):sortpace2(k,1)) ,1)) / mean(val1( (sortpace1(k-1,1):sortpace1(k,1)) ,1));
    for j = sortpace2(k-1,1):sortpace2(k,1)
        normgene = 1/frac * val2( (sortpace2(k-1,1):sortpace2(k,1)) ,1);
    end
end
for k = 1:45101
    finalnormgene(idx2(k,1),1) = normgene(k);
end
end
     
 


    