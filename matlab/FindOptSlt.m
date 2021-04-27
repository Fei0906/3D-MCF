function [slt]=FindOptSlt(A,b,sph)
%=======================================================================
%   Find the optimal solution of multiple closures.
%
%   Fei Liu, Oct 2019
%
%
%   Input:
%       A: closure matrix. It is a X*p matrix, where X stands for the total
%           number of closures (include the redundant) and p stands for the
%           number of interferograms.
%       b: round values of sum of wrapped phases. It is a 1*X vector.
%       sph: absolute values of sum of wrapped phases. It is a 1*X vector.
%
%   Output:
%       slt: the solution of multiple closures. It is a 1*rA vector, where
%           rA is the rank of A.
%=======================================================================

[X,~]=size(A);
rA=rank(A);

[A2,idx]=licols(A');
A2=A2';

A3=zeros(X-rA,rA);
k=1;
for i=1:X
    if isempty(find(idx==i,1))
        C=A(i,:);
        A3(k,:)=round(C/A2);
        k=k+1;
    end
end
        
A4=zeros(X,rA);
A4(1:rA,1:rA)=eye(rA);
C=zeros(X,1);
C(1:rA)=b(idx);
A4(rA+1:end,:)=A3;
C(rA+1:end)=b(setdiff(1:X,idx));
        
Pr=1./((b-sph).*(b-sph));
Pr2=zeros(X,1);
Pr2(1:rA)=Pr(idx);
Pr2(rA+1:end)=Pr(setdiff(1:X,idx));
        
slt=FindMinRes(A4,Pr2,C);

end






