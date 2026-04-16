function out=calcRDMs(Y,SPM,condition,distType)
% Y is usually full(X(T.LI{i},:)') from rsa.runSearchlight.m: raw data,
% only withing the current searchlight.

% spatial prewhitening: 
% Do GLM and Get prewhitened beta weights
B=real(rsa.spm.noiseNormalizeBeta(Y,SPM));
%take only betas of relevant conditions
B=B(logical(condition),:);

% added bit: make sure that empty regressors have betas that are exactly 0
B(abs(B)<10^-8)=0;
warning('off', 'all');     % Turn off all warnings
% calculate distance
RDM = pdist(B,distType);    
warning('on', 'all');     % Turn off all warnings
% meanDist=mean(RDM,2);           % Calulate the mean distance
RDM=RDM';
% out = [RDM(:);meanDist(:)];     % Arrange the outputs in a vector, as they are written to files
out = RDM(:);     % Arrange the outputs in a vector, as they are written to files
