function []=JointUnwrap(maxIter,tolFun)
%=========================================================================
%
%   Fei Liu, Mar 2020
%
%   The main function to do phase unwrapping. And unwrap all the 
%   interferograms jointly
%   
%   Feb 2021, change the way to construct parameter matrix to run faster
%=========================================================================

%===================
%   Load data
%===================
if nargin<1
    maxIter=50;
    tolFun=1e-6;
end

data1=load('Tmpdata');
data2=load('TIN');

ph=data1.ph;
dk=data1.dk;
dk_s=data1.dk_s;
weight=data1.weight.^2;
res=data1.res;
edges=data2.edges_nz;
eles=data2.eles_nz;

n_ele=size(eles,1);
n_edge=size(edges,1);
n_ifg=size(ph,2);
n_pixel=size(ph,1);
n_closure=size(dk_s,2);

clear data1 data2;


%===================
%   Form Matrix
%===================

%   Object function C
w1=weight(edges(:,1),:);
w2=weight(edges(:,2),:);
corr=min(w1,w2);
C=[corr;corr];
C=double(C(:));
C(C<0.1)=0.1;
clear w1 w2 corr;

%   Matrix A
% for i=1:n_ele
%     index=eles(i,4:6);
%     for j=1:3
%         if index(j)>0
%             a11(i,index(j))=1;
%             a11(i,n_edge+index(j))=-1;
%         else
%             a11(i,-index(j))=-1;
%             a11(i,n_edge-index(j))=1;
%         end
%     end
% end
loci=zeros(1,6*n_ele);
for i=1:6
    loci(1,i:6:end)=1:n_ele;
end
locj=zeros(1,6*n_ele);
locj(1,1:2:end)=abs(reshape(eles(:,4:6)',3*n_ele,1));
locj(1,2:2:end)=n_edge+locj(1,1:2:end);
locv=zeros(1,6*n_ele);
idx1=find(reshape(eles(:,4:6)',3*n_ele,1)>0);
idx2=find(reshape(eles(:,4:6)',3*n_ele,1)<0);
locv(1,2*idx1-1)=1;
locv(1,2*idx1)=-1;
locv(1,2*idx2-1)=-1;
locv(1,2*idx2)=1;
a11=sparse(loci,locj,locv,n_ele,2*n_edge);

t1=repmat({a11},n_ifg,1);
Au=blkdiag(t1{:});

% Ad=sparse(n_edge*n_closure,2*n_ifg*n_edge);
% t2=speye(n_edge,n_edge);
% a22=[t2 -t2];
% for i=1:n_closure
%     index=dk_s{1,i};
%     for j=1:size(index,2)
%         if index(j)>0
%             Ad((i-1)*n_edge+1:n_edge*i,(index(j)-1)*2*n_edge+1:index(j)*2*n_edge)=a22;
%         else
%             Ad((i-1)*n_edge+1:n_edge*i,(-index(j)-1)*2*n_edge+1:-index(j)*2*n_edge)=-a22;
%         end
%     end
% end
locidx=1;
for i=1:n_closure
    index=dk_s{1,i};
    for j=1:size(index,2)
        starti=(i-1)*n_edge+1;
        startj=(abs(index(j))-1)*2*n_edge+1;
        loci(1,locidx:locidx+n_edge-1)=starti:starti+n_edge-1;
        locj(1,locidx:locidx+n_edge-1)=startj:startj+n_edge-1;
        locv(1,locidx:locidx+n_edge-1)=abs(index(j))/index(j);
        locidx=locidx+n_edge;
        
        startj2=startj+n_edge;
        loci(1,locidx:locidx+n_edge-1)=starti:starti+n_edge-1;
        locj(1,locidx:locidx+n_edge-1)=startj2:startj2+n_edge-1;
        locv(1,locidx:locidx+n_edge-1)=-abs(index(j))/index(j);
        locidx=locidx+n_edge;
    end
end
Ad=sparse(loci,locj,locv,n_edge*n_closure,2*n_ifg*n_edge);

A=[Au;Ad];
clear a11 a22 t1 t2 Au Ad;

%   Matrix b
b=zeros(n_ele*n_ifg+n_edge*n_closure,1);
b(1:n_ele*n_ifg)=double(res(:));
for i=1:n_closure
    index=dk_s{1,i};
    esum=zeros(n_edge,1);
    for j=1:size(index,2)
        if index(j)>0
            esum=esum+dk(:,index(j));
        else
            esum=esum-dk(:,-index(j));
        end
    end
    b(n_ele*n_ifg+(i-1)*n_edge+1:n_ele*n_ifg+i*n_edge)=-dk_s{2,i}(:)+esum(:);
end

%   Lower limitation
lb=sparse(2*n_ifg*n_edge,1);

%   Clear some data which no longer need
clear dk_s esum index res weight closure

%===================
%   Optimization
%===================
disp('Start Optimization...');
options=optimset('Display','iter','Algorithm','interior-point','MaxIter',maxIter,'TolFun',tolFun);

[x]=linprog(C,[],[],A,b,lb,[],[],options);
if isempty(x)
	msg='The algorithm of linprog  cannot output its optimal solution.';
	error(msg);
else
    disp('Optimization Done!');
end
x=round(x);
save('x','x');
clear A b C lb;

%==========================
%   Get unwrapped result
%==========================
K=zeros(n_edge,n_ifg);
for i=1:n_ifg
    K(:,i)=x((i-1)*2*n_edge+1:(i-1)*2*n_edge+n_edge)-x((i-1)*2*n_edge+n_edge+1:i*2*n_edge);
end
K=double(-K+dk);

res_K=zeros(n_ele,n_ifg);
for i=1:n_ele
    for j=1:3
        if eles(i,3+j)>0
            res_K(i,:)=res_K(i,:)+K(eles(i,3+j),:);
        else
            res_K(i,:)=res_K(i,:)-K(-eles(i,3+j),:);
        end
    end
end

A=sparse(n_edge,n_pixel);
for i=1:n_edge
	A(i,edges(i,1))=1;
	A(i,edges(i,2))=-1;
end

A=A(:,2:end);
k_map=lscov(A,K);

k_map=round([zeros(1,n_ifg);k_map]);
uph=(ph+k_map)*2*pi;

save('uph','uph');
clear;

end