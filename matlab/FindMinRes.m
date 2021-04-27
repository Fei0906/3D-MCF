function [x]=FindMinRes(A,P,b)
%=======================================================================
%   Find the optimal solution of Ax=b by brute force
%
%
%   Input:
%       A: The input matrix in left
%       P: The weight matrix
%       b: The input matrix in right
%
%   Output:
%       x: The optimal result
%=======================================================================
[~,n]=size(A);

x=b(1:n);
mx=repmat(x,1,3*n);

for i=1:n
    mx(i,3*i-2:3*i)=b(i)-1:b(i)+1;
end

min=1000;

for i=1:3*n
    res=P'*abs(A*mx(:,i)-b);
    if res<min
        min=res;
        x=mx(:,i);
    end
end
end