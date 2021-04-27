function [Res,dk]=CalRes(ph,DATATYPE)
%=======================================================================
%   Calculate the estimated UPGs and residual matrix in each interferograms   
%
%   Fei Liu, Nov 2019
%
%   Input:
%       ph: wrapped phases, whose size is M*N*P in continuous case and
%           S*P in discrete case.
%
%   Output:
%       Res: the residual matrix, whose size is (M-1)*(N-1)*P in continous
%           case and TriNum*P (TriNum is the number of triangle in TIN) in 
%           discrete case.
%       dk: the estimated UPGs in space. In continuous case, it is a 2*1
%           cell. The first row store the UPGs in x direction, whose size
%           is M*(N-1)*P, and the second row store the UPGs in y direction,
%           whose size is (M-1)*N*P. And in discrete case, it is a matrix,
%           whose size is EdgesNum*P (EdgesNum is the number of Edges in 
%           TIN).
%=======================================================================

%Estimate dk_x and dk_y

if DATATYPE=='C'
    %   The estimated UPGs in space is acquired by assuming that phase 
    %   differences between nearby pixel is less than pi
    dk_x=diff(ph,1,2);
    index0= dk_x<0.5&dk_x>=-0.5;
    index1= dk_x<-0.5;
    indexm1= dk_x>=0.5;
    dk_x(index0)=0;
    dk_x(index1)=1;
    dk_x(indexm1)=-1;

    dk_y=diff(ph,1,1);
    index0= dk_y<0.5&dk_y>=-0.5;
    index1= dk_y<-0.5;
    indexm1= dk_y>=0.5;
    dk_y(index0)=0;
    dk_y(index1)=1;
    dk_y(indexm1)=-1;
    
    dk=cell(2,1);
    dk{1,1}=dk_x;
    dk{2,1}=dk_y;
    
    %   Get residual matrix according to dk_x and dk_y
    [m,n,p]=size(ph);
    Res=zeros(m,n,p,'single');
    Res(1:m-1,1:n-1,:)=dk_x(1:m-1,:,:)-dk_x(2:m,:,:)+dk_y(:,2:n,:)-dk_y(:,1:n-1,:);
else
    tdata=load('TIN.mat');
    edges=tdata.edges_nz;
    eles=tdata.eles_nz;
    
    dk=ph(edges(:,1),:)-ph(edges(:,2),:);
    index0= dk<0.5&dk>=-0.5;
    index1= dk<-0.5;
    indexm1= dk>=0.5;
    dk(index0)=0;
    dk(index1)=1;
    dk(indexm1)=-1;
    
    [~,p]=size(ph);
    [n_ele,~]=size(eles);
    Res=zeros(n_ele,p);
    for i=1:n_ele
        for j=1:3
            if eles(i,3+j)>0
                Res(i,:)=Res(i,:)+dk(eles(i,3+j),:);
            else
                Res(i,:)=Res(i,:)-dk(-eles(i,3+j),:);
            end
        end
    end
end

end


