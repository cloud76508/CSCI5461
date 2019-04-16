1. Download Data from TCGA database:

Gene expression data:
https://tcga.xenahubs.net/download/TCGA.PRAD.sampleMap/HiSeqV2.gz

Somatic mutation data:
https://tcga.xenahubs.net/download/TCGA.PRAD.sampleMap/mutation_broad_gene.gz

phenotype data:
https://tcga.xenahubs.net/download/TCGA.PRAD.sampleMap/PRAD_clinicalMatrix.gz


2. After all data sets are prepared, please follow the workflow to execute the scripts:

Wrokflow(for gene expression):

Phenotype.m => gene_expression.m => DataPreprocess_expression => knnCode.m (KNN) or svmCode.m (linear SVM) or svmCode_nonlinear.m (non-linear SVM)

check the prediction accuracy (10 repetitions) in Accuracy.txt

Wrokflow(for gene expression):

Phenotype.m => mutation.m => DataPreprocess_mutation => knnCode.m or svmCode.m or svmCode_nonlinear.m

check the prediction accuracy (10 repetitions) in Accuracy.txt


3. Determine the genetic type and the number of features for model estimation in knnCode.m / svmCode.m / svmCode_nonlinear.m (line 2 - 12)

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