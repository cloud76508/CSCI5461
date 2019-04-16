clear all
clc
load('Phenotype.mat')
load('GeneExpression.mat')

% Data match phenotype (GeneExpression)
DataGE.x = [];
DataGE.y = [];
for n =1:size(phenotype,1)
   idx = [];
   idx = find(strcmp(SampleName,phenotype(n,1)) == 1);
   if isempty(idx)
       DataGE.x(n,:) = zeros(1,size(Data,1));
   elseif isnan(phenotype{n,2})
       DataGE.x(n,:) = zeros(1,size(Data,1));
   else
       DataGE.x(n,:) = Data(:,idx)'; 
   end
end

%delete unmatched data
phenotype(sum(DataGE.x,2)==0,:) = [];
DataGE.x(sum(DataGE.x,2)==0,:) = [];

%Y is the five-year recurrence status
for n =1:size(phenotype,1)
   if  phenotype{n,3} == 1 && phenotype{n,2} < 5*365
       DataGE.y(n,1) = 1;
   else
       DataGE.y(n,1) = 0;
   end
end

%Data normalization
for n =1:size(DataGE.x,2)
    MinV = min(DataGE.x(:,n));
    MaxV = max(DataGE.x(:,n));
    Diff = MaxV-MinV;
    for m = 1:size(DataGE.x,1)
        tempData(m,n) = (DataGE.x(m,n) - MinV)/Diff;
    end
end
DataGE.x = tempData;

i=1;
j=1;
for n =1:size(DataGE.y,1)
    if DataGE.y(n,1) == 1
       Recc.x(i,:) =  DataGE.x(n,:);
       i=i+1;
    else
       NonRecc.x(j,:) = DataGE.x(n,:);
       j=j+1;
    end
end

[~,p1] = ttest2(Recc.x,NonRecc.x);
[~,p_ind]=sort(p1,'ascend');
clearvars -except DataGE GeneName phenotype p1 p_ind
save('LearningSet_GE.mat')

