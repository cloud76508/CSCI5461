clear all
clc
load('Phenotype.mat')
load('Mutation.mat')

% Data match phenotype (GeneExpression)
DataMU.x = [];
DataMU.y = [];
for n =1:size(phenotype,1)
   idx = [];
   idx = find(strcmp(SampleName,phenotype(n,1)) == 1);
   if isempty(idx)
       DataMU.x(n,:) = zeros(1,size(Data,1));
   elseif isnan(phenotype{n,2})
       DataMU.x(n,:) = zeros(1,size(Data,1));
   else
       DataMU.x(n,:) = Data(:,idx)'; 
   end
end

%delete unmatched data
phenotype(sum(DataMU.x,2)==0,:) = [];
DataMU.x(sum(DataMU.x,2)==0,:) = [];

%Y is the five-year recurrence status
for n =1:size(phenotype,1)
   if  phenotype{n,3} == 1 && phenotype{n,2} < 5*365
       DataMU.y(n,1) = 1;
   else
       DataMU.y(n,1) = 0;
   end
end

i=1;
j=1;
for n =1:size(DataMU.y,1)
    if DataMU.y(n,1) == 1
       Recc.x(i,:) =  DataMU.x(n,:);
       i=i+1;
    else
       NonRecc.x(j,:) = DataMU.x(n,:);
       j=j+1;
    end
end

[~,p1] = ttest2(Recc.x,NonRecc.x);
[~,p_ind]=sort(p1,'ascend');
clearvars -except DataMU GeneName phenotype p1 p_ind
save('LearningSet_MU.mat')
