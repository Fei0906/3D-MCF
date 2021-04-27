function H = sfd(x,grad,H,group,alpha,DiffMinChange,DiffMaxChange, ...
                  funfcn,varargin)
%

%SFD Sparse Hessian via finite gradient differences
%
% H = sfd(x,grad,H,group,alpha,DiffMinChange,DiffMaxChange,funfcn,varargin)
% returns the sparse finite difference approximation H of the Hessian 
% matrix of the function 'funfcn' at the current point x. The vector group 
% indicates how to use sparse finite differencing: group(i) = j means that 
% column i belongs to group (or color) j. Each group (or color) corresponds
% to a finite gradient difference.
%
% A non-empty input alpha overrides the default finite differencing 
% stepsize.
% 

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/09 01:06:10 $

xcurr = x(:);  % Preserve x so we know what funfcn expects
scalealpha = false;
[m,n] = size(H); 
v = zeros(n,1);
ncol = max(group); 
if isempty(alpha)
    scalealpha = true;
    alpha = repmat(sqrt(eps),ncol,1);
end
H = spones(H); 
for k = 1:ncol
    d = (group == k);
    if scalealpha 
        xnrm = norm(xcurr(d));
        xnrm = max(xnrm,1);
        alpha(k) = alpha(k)*xnrm;
    end

    % Ensure magnitude of step-size lies within interval [DiffMinChange, DiffMaxChange]
    alpha(k) = sign(alpha(k))*min(max(abs(alpha(k)),DiffMinChange), ...
                                  DiffMaxChange);
    
    y = xcurr + alpha(k)*d;    
    
    % Make x conform to user-x
    x(:) = y;

    switch funfcn{1}
    case 'fun'
        error(message('optim:sfd:ShouldNotReachThis'))
    case 'fungrad'
        [dummy,v(:)] = feval(funfcn{3},x,varargin{:});
    case 'fun_then_grad'
        v(:) = feval(funfcn{4},x,varargin{:});
    otherwise
        error(message('optim:sfd:UndefCallType'));
    end
    
    
    w = (v-grad)/alpha(k);
    cols = find(d); 
    A = sparse(m,n);
    A(:,cols) = H(:,cols);
    H(:,cols) = H(:,cols) - A(:,cols);
    [i,j,val] = find(A);
    [p,ind] = sort(i);
    val(ind) = w(p);
    A = sparse(i,j,full(val),m,n);
    H = H + A;
end
H = (H+H')/2;   % symmetricize



