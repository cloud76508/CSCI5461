%phenotype
%C1: smaple ID
%C2: first recurrence day (or last followup, if non-recurrence )
%C3: reccurence status (1: reccurence, 0: non-reccurence)

clear all
clc
tdfread('PRAD_clinicalMatrix') %phenotype data

clearvars -except days_to_first_biochemical_recurrence sample_type_id sampleID days_to_last_followup ...
    biochemical_recurrence

phenotype = {};
for n =1:size(sampleID,1)
    phenotype(n,1) = cellstr(sampleID(n,:));
    phenotype(n,2) = num2cell(sample_type_id(n,:));
    phenotype(n,3) = num2cell(days_to_first_biochemical_recurrence(n,:));
    phenotype(n,4) = num2cell(days_to_last_followup(n,:));
end

%delete normal sample or without recurrence status
phenotype(sample_type_id == 11,:) = [];
phenotype(:,2) = [];
for n =1:size(phenotype,1)
    if isnan(phenotype{n,2}) == 0 
       phenotype{n,3} = phenotype{n,2};
       phenotype{n,4} = 1;
    else
       phenotype{n,4} = 0;
    end
end

phenotype(:,2) = [];

% %delete censored sample or without recurrence status(extra test)
% for n =1:size(phenotype,1)
%     if phenotype{n,3} == 0 && phenotype{n,2} <= 5*365
%        phenotype{n,3} = 3;
%     end
% end
% i=1;
% for n=1:size(phenotype,1)
%     if phenotype{n,3} ~=3
%        tempphenotype(i,:) = phenotype(n,:) ;
%        i=i+1;
%     end
% end
% phenotype = tempphenotype;

       
clearvars -except phenotype
save('Phenotype.mat')

