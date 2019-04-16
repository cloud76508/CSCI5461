clear all
%--------------------------------------------------------------------------
% Determine the data type (GE or MU) and the number of features for model
% estimation
%--------------------------------------------------------------------------
load('LearningSet_GE.mat') %for gene expression data
Data = DataGE; %for gene expression data
% load('LearningSet_MU.mat') %for somatic mutation data
% Data = DataMU; %for somatic mutation data

FeatureNumber = 200; %the number of features for model estimation
%--------------------------------------------------------------------------

for Exp = 1:10
    
clearvars -except Data p_ind Exp ExpResults FeatureNumber
clc

Recc = {};
NonRecc = {};
i=1;
j=1;

for n =1:size(Data.y,1)
    if Data.y(n,1) == 1
       Recc.x(i,:) =  Data.x(n,p_ind(1:FeatureNumber));
       i=i+1;
    else
       NonRecc.x(j,:) = Data.x(n,p_ind(1:FeatureNumber));
       j=j+1;
    end
end

Recc.y = zeros(size(Recc.x,1),1)+1;
NonRecc.y = zeros(size(NonRecc.x,1),1);

RandomIndex0 = randperm(size(Recc.y,1));
tempdRecc = Recc.x(RandomIndex0,:);
Recc.x = tempdRecc;

RandomIndex = randperm(size(NonRecc.y,1));

SubNonRecc.x = NonRecc.x(RandomIndex(1:size(Recc.y,1)),:);
SubNonRecc.y = zeros(size(SubNonRecc.x,1),1);

for  fold =1:3
    test.x = [Recc.x((fold-1)*size(Recc.y,1)/3+1 : (fold-1)*size(Recc.y,1)/3+size(Recc.y,1)/3,:) ...
        ; SubNonRecc.x((fold-1)*size(Recc.y,1)/3+1 : (fold-1)*size(Recc.y,1)/3+size(Recc.y,1)/3,:)];
    test.y = [zeros(size(Recc.y,1)/3,1)+1;zeros(size(Recc.y,1)/3,1)];
    Temp1 = Recc.x;
    Temp2 = SubNonRecc.x;
    Temp1((fold-1)*size(Recc.y,1)/3+1 : (fold-1)*size(Recc.y,1)/3+size(Recc.y,1)/3,:) = [];
    Temp2((fold-1)*size(Recc.y,1)/3+1 : (fold-1)*size(Recc.y,1)/3+size(Recc.y,1)/3,:) = [];
    train.x = [Temp1;Temp2];
    train.y = [zeros(size(Recc.y,1)/3*2,1)+1;zeros(size(Recc.y,1)/3*2,1)];
    
    k = 1:20;
    for idx = 1:length(k)
        crossAcc = [];
        crossLabel = [];
        for foldInner = 1:2
            TempAcc = [];
            TempLabel = [];
            Temp1 = train.x(1:size(train.x,1)/2,:);
            Temp2 = train.x(size(train.x,1)/2+1:end,:);
            Temp1 = Temp1((foldInner-1)*size(Temp1,1)/2+1 : (foldInner)*size(Temp1,1)/2, :);
            Temp2 = Temp2((foldInner-1)*size(Temp2,1)/2+1 : (foldInner)*size(Temp2,1)/2, :);
            learn.x = [Temp1; Temp2];
            learn.y = [zeros(size(Recc.y,1)/3,1)+1;zeros(size(Recc.y,1)/3,1)];

            Temp1 = train.x(1:size(train.x,1)/2,:);
            Temp2 = train.x(size(train.x,1)/2+1:end,:);
            Temp1((foldInner-1)*size(Temp1,1)/2+1 : (foldInner)*size(Temp1,1)/2, :) = [];
            Temp2((foldInner-1)*size(Temp2,1)/2+1 : (foldInner)*size(Temp2,1)/2, :) = [];
            valid.x = [Temp1;Temp2];
            valid.y = [zeros(size(Recc.y,1)/3,1)+1;zeros(size(Recc.y,1)/3,1)];
            
            model = fitcknn(learn.x,learn.y,'NumNeighbors',k(idx),'Standardize',1);
            label = predict(model,valid.x);
            
            crossAcc(foldInner,1) = sum(label==valid.y)/size(valid.y,1);
        end
        ValidAcc(idx,1) =   mean(crossAcc);
    end
    [~, optimal_idx] = max(ValidAcc);
    
    optimal_model = fitcknn(train.x,train.y,'NumNeighbors',k(optimal_idx),'Standardize',1);
    label = predict(optimal_model,test.x);
    acc = sum(label==test.y)/size(test.y,1);
    FinalAcc(fold) = acc(1,1);
end
%mean(FinalAcc)
%k(optimal_idx)

ExpResults(Exp,1) = mean(FinalAcc);

end

fileID = fopen('Accuracy.txt','w');
for Exp = 1:10
    fprintf(fileID,'%.4f\n',ExpResults(Exp,1));
end
fclose(fileID);
