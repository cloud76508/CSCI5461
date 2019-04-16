clear all;clc
% load data into MATLAB
fileID = fopen('mutation_broad_gene');
tline = fgets(fileID);
SampleName = strsplit(tline);
SampleName = SampleName(2:end)';
SampleName(end) = [];
fmt = ['%s' repmat('%f', 1, size(SampleName,1))];
C = textscan(fileID,fmt);
GeneName = C{1};
Data = cell2mat(C(2:end)); % Data = sparse(Data), whos()
fclose(fileID);

GeneName(1,:) = [];
Data(1,:) = [];

clearvars -except Data GeneName SampleName
save('Mutation.mat')