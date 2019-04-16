clear all;clc
% load data into MATLAB
fileID = fopen('HiSeqV2');
tline = fgets(fileID);
SampleName = strsplit(tline);
SampleName = SampleName(2:end)';
SampleName(end)=[];
fmt = ['%s' repmat('%f', 1, size(SampleName,1))];
C = textscan(fileID,fmt);
GeneName = C{1};
Data = cell2mat(C(2:end));
fclose(fileID);

clearvars -except Data GeneName SampleName
save('GeneExpression.mat')
