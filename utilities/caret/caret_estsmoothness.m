function [FWHM,resel]=caret_estsmoothness(cSPM,S)
% function [FWHM,resel]=caret_estsmoothness(cSPM,S)
% Estimates the smoothness in the principle directions  
% either on a flatmap (2D) or on Fiducial map (along 3D)
% The effective smoothness is F=prod(FWHM).^(1/D)
% In estimating it restrticts itself to pairs of Nodes,
% that have positive area-surface to them 
% ---------------------------------------------------
% v.1.0 Joern Diedrichsen jdiedric@bme.jhu.edu
% jdiedric@bme.jhu.edu

% ---------------------------------------------------
% get data and remove dimesions that have zero variance (e.g. Z in
% flatmaps)
X=S.Nodes.data;
a=find(var(X)>0);
X=X(:,a);
[N,D]=size(X);

% ---------------------------------------------------
% calculate the residuals of the linear model and 
% standardize them 
Z=cSPM.data;
[N,n]=size(Z);
Z=(Z-cSPM.b*cSPM.X');
STD=sqrt(nansum(Z.^2,2));
Z=Z./repmat(STD,1,n);

if (~isfield(S.Tiles,'good'))
    S.Tiles.good=true(S.num_tiles,1); 
end;

% ---------------------------------------------------
% Generate all possible M pairs of nodes 

Z1=Z(S.Edges.data(:,1),:); 
Z2=Z(S.Edges.data(:,2),:); 
X1=X(S.Edges.data(:,1),:); 
X2=X(S.Edges.data(:,2),:); 

Zd=[Z1-Z2];
Xd=[X1-X2];

j=find(~isnan(sum(Zd,2)) & sqrt(sum(Xd.^2,2))>0.1); 
Zd=[Zd(j,:);-Zd(j,:)]; 
Xd=[Xd(j,:);-Xd(j,:)]; 
Z=[Z1(j,:);Z2(j,:)];
d=sqrt(sum(Xd.^2,2));       


Zd=Zd./repmat(d,1,n);
Xd=Xd./repmat(d,1,D);       % Notice that if d==0 Xd can become NaN 
% d=mean(sqrt(sum(Xd.^2,2)));       
% Zd=Zd./d;
% Xd=Xd./d;       % Notice that if d==0 Xd can become NaN 


% ---------------------------------------------------
% generate vech(X'*X) 
if (D==2) 
    H=[Xd.^2 Xd(:,1).*Xd(:,2)];
elseif (D==3)
    H=[Xd.^2 Xd(:,1).*Xd(:,2) Xd(:,2).*Xd(:,3) Xd(:,1).*Xd(:,3)];    
else
    error('Map must be in 2D or 3D');
end;

% ---------------------------------------------------
% Only use in the calculation where we have complete data 
indx=find(~any(isnan(Zd)') & ~any(isnan(Xd)'));
DE=(n-2)./(n-1).*sum(Zd(indx,:).^2,2);
H=H(indx,:);
L=inv(H'*H)*H'*DE; 

% ---------------------------------------------------
% Reassemble the Smoothness matrix and take eigenvalues 
if (D==2)
    L=[L(1) L(3);L(3) L(2)];
else
    L=[L(1) L(4) L(5);L(4) L(2) L(6);L(5) L(6) L(3)];
end;
FWHM=sqrt(4*log(2)./eig(L));
resel=prod(1./FWHM);
