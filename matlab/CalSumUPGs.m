function [dk_s]=CalSumUPGs(ph,closure,DATATYPE)
%=======================================================================
%   Calculate the sum of UPGs in closure phase.
%
%   Fei Liu, Oct 2019
%
%
%   Input:
%       ph: wrapped phases, whose size is M*N*P in continuous data and S*P
%           in discrete data.
%       closure: closure phase, which is a 1*X cell. And X stands for the
%           total number of closures (include the redundant).
%       DATATYPE: continous or discrete data.
%
%   Output:
%       dk_s: sum of UPGs.In continous case, it is a 3*X cell. The first 
%           row stores the closures information. The second row stores
%           the sum of UPGs in horizontal direction. And the third row
%           stores the sum of UPGs in veritcal direction. In discrete case,
%           it is a 2*X cell. The second row stores the UPGs. All the 
%           elements in dk_s are integers.
%=======================================================================


[~,n_closure]=size(closure);
warning('off');

%   Calculate sum of wrapped phases
if DATATYPE=='C'
    [m,n,n_ifg]=size(ph);
    A=zeros(n_closure,n_ifg);
    sph=zeros(m,n,n_closure);
    for i=1:n_closure
        closure_i=closure{1,i};
        ph_sum=zeros(m,n);
        for i=1:size(closure_i,2)
            if closure_i(i)<0
                ph_sum=ph_sum-ph(:,:,-closure_i(i));
                A(i,-closure_i(i))=-1;
            else
                ph_sum=ph_sum+ph(:,:,closure_i(i));
                A(i,closure_i(i))=1;
            end
        end
        sph(:,:,i)=ph_sum;
    end
else
    [n_pixel,n_ifg]=size(ph);
    A=zeros(n_closure,n_ifg);
    sph=zeros(n_pixel,n_closure);
    for i=1:n_closure
        closure_i=closure{1,i};
        ph_sum=zeros(n_pixel,1);
        for j=1:size(closure_i,2)
            if closure_i(j)<0
                ph_sum=ph_sum-ph(:,-closure_i(j));
                A(i,-closure_i(j))=-1;
            else
                ph_sum=ph_sum+ph(:,closure_i(j));
                A(i,closure_i(j))=1;
            end
        end
        sph(:,i)=ph_sum;
    end
end

%	Correct pixels for multiple closures
rA=rank(A);
[~,idx]=licols(A');
tn=0;

if DATATYPE=='C'
    new_sph=zeros(m,n,rA);
    for i=1:m
        for j=1:n
            b=round(sph(i,j,:));
            if rank([A b(:)])~=rA
                tn=tn+1;
                new_sph(i,j,:)=FindOptSol(A,b,sph(i,j,:));
            else
                new_sph(i,j,:)=b(idx);
            end
        end
    end
else
    new_sph=zeros(n_pixel,rA);
    for i=1:n_pixel
        b=round(sph(i,:));
        if rank([A b(:)])~=rA
            tn=tn+1;
            new_sph(i,:)=FindOptSlt(A,b,sph(i,:));
        else
            new_sph(i,:)=b(idx);
        end
    end
end
msg=['Total correction pixels number are ' num2str(tn) '.'];
disp(msg);

%   Remove the redundant closures and modify the value of X
n_closure=rA;
sph=new_sph;
new_closure=cell(1,n_closure);
for i=1:n_closure
    new_closure{1,i}=closure{1,idx(i)};
end
closure=new_closure;
clear new_sph new_closure;

%   Acquire the sum of UPGs
if DATATYPE=='C'
    dk_s=cell(3,n_closure);
    for i=1:n_closure
        dk_s{1,i}=closure{1,i};
        dk_s{2,i}=-diff(sph(:,:,i),1,2);
        dk_s{3,i}=-diff(sph(:,:,i),1,1);    
    end
else
    dk_s=cell(2,n_closure);
    tdata=load('TIN.mat');
    edges=tdata.edges_nz;
    for i=1:n_closure
        dk_s{1,i}=closure{1,i};
        dk_s{2,i}=-sph(edges(:,1),i)+sph(edges(:,2),i);
    end
end

end