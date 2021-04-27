function [gt,d]=findmax(data, flag)
%FINDMAX Interpolates the maxima in a vector of data.
%
%   Function used for returning all the maxima in a set of
%   data (vector or matrices). The maxima are calculated by
%   cubic or quadratic interpolation.

%   Copyright 1990-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2009/01/23 21:40:28 $

[m,n]=size(data);
if nargin < 2, flag = 1; end

if isempty(data)
    gt = [];
    d = [0.5 0.5];
    return;
elseif m==1 || n==1
    % 1-D
    d = [2, 2];
    data=data(:);
    n=max([n,m]);
    ind=find( (data>[data(1)-1;data((1:n-1)')]) ...
           & (data>=[data((2:n)');data(n)-1]) );
    s=length(ind);
    gt=zeros(s,1);
    for i=1:s
        ix=ind(i);
        factor = 1;
        err = 0;  % Smoothness factor
        if (ix==1 || ix==n);
            gt(i)=data(ix);
        elseif ix>n-3 || ix<4
            if data(ix-1)==data(ix) && data(ix+1)==data(ix)
                gt(i)=data(ix);
            else
                gt(i)=quad2(data(ix-1:ix+1));
            end
        elseif data(ix+2)<data(ix+1) 
            [gt(i), err] = cubic(data(ix-1:ix+2), data(ix-2), 0);
        elseif data(ix-2)<data(ix-1) 
            gt(i)=cubic(data(ix-2:ix+1));
        elseif data(ix)>data(ix+1)
            gt(i)=quad2(data(ix-1:ix+1));
            factor = 0.75; 
        else
            gt(i)=data(ix);
            factor = 0.5;
        end

% Discretization interval estimation
% Aim: keep error less than 1/10 of abs(constraint) greater than 0
%   or  1/50 of smoothness factor.
% error = data(ix) - gt(i)
        if flag 
            if ~err
                err = data(ix) - gt(i);
            end
            f = abs((gt(i) + 1.1e-5) /(10 * err + 1e-5));
            d(1,1) = min(d(1,1), factor * ((f>=0.3) + (0.7 + f)*(f< 0.3) +(f>1)*(log(f)/100)) );
        end
    end
else
% 2D
    [gt, d] = findmax2(data, flag); 
end
if (m*n)<4 
    d=[0.5,0.5]; 
end

%======================== SUBFUNCTIONS ====================================

function [gt,d]=findmax2(data, flag)
%FINDMAX2 Interpolates the maxima in a matrix of data.
%
%   Function used for returning all the maxima in a set of
%   data (matrices). The maxima are calculated by
%   quadratic interpolation.

if nargin < 2, flag = 0; end
% 2-D
d = [1.5, 1.5];
[m,n]=size(data);
% Factor for error control:    
factor = 5;  % Keep error five times less than nearest constraint

% Use point [0,0], [1,0], [0,1], [0, -1], [-1 0] and [1,1] on grid:
diagix3 = [0 1 0 0 -1 1 ];
diagix4 = [0 0 1 -1 0 1 ];
% Point [-1, -1] used to check goodness of quadratic fit. 
xerr = [-1; -1];

invH = [  ... 
-1.0000    0.5000         0         0    0.5000         0
0.5000   -0.5000   -0.5000         0         0    0.5000
-1.0000         0    0.5000    0.5000         0         0
0    0.5000         0         0   -0.5000         0
0         0    0.5000   -0.5000         0         0
1.0000         0         0         0         0         0
];

%%% Start of peaks calculation %%%
% [1] Find the peaks. 
crmax = (data>[data(1,:)-1;data(1:m-1,:)]) & (data>=[data(2:m,:);data(m,:)-1]) & ...
    (data>[data(:,1)-1,data(:,1:n-1)]) & (data>=[data(:, 2:n),data(:,n)-1]);

f = find(crmax);

% [2] Check the peaks. The calculation above finds elements of "data", P,
% (possible peaks), such that P is greater than each of the elements, S,
% shown in the schematic below.
%        S 
%      S P S
%        S
% The following calculation checks each P to ensure that the element of
% data at P is also larger than the data at the adjacent corner cells, C,
% that is
%      C S C  
%      S P S
%      C S C
for i = f'
    
    % Row index of i-th peak
    mi = 1 + rem(i-1, m);
    
    % Find the adjacent corner cells to the i-th peak.
    if mi == 1
        % Top row of matrix. Only consider corner cells in second row.
        ic = [i-m+1 i+m+1];
    elseif mi == m
        % Bottom row of matrix. Only consider corner cells in (m-1)-th row.
        ic = [i-m-1 i+m-1];
    else
        % Interior row of matrix. Can consider all corners.
        ic = [i-m+1 i+m+1 i+m-1 i-m-1];
    end
    % Filter out any corners that are in the zero-th or (n+1)-th column.
    ic = ic(ic > 0 & ic <= n*m);
    
    % If the data value of any of the adjacent corner cells is greater than
    % or equal to that of the i-th peak, then the i-th peak is not a
    % maximum.
    if any(data(i) <= data(ic))
        crmax(i) = 0;
    end
    
end
%%% End of peaks calculation %%%

% Find the indices of the remaining peaks
f = find(crmax)';

% Initialise return argument
gt = [];

% Fit a quadratic
H =zeros(2,2); 
cnt = 0; 
for i = f
    cnt = cnt + 1;
    ni = 1 + floor((i-0.5)/m);
    mi = i - (ni-1)*m;
    if (mi > 1 & mi < m & ni > 1 & ni < n) 
        ix = mi - diagix3;
        iy = ni - diagix4;
        for k = 1:6
            y(k,1) = data(ix(k),iy(k)); 
        end
        coeffs = invH*y; 
        H(1,1) = coeffs(1);
        H(1,2) = coeffs(2);
        H(2,1) = H(1,2);
        H(2,2) = coeffs(3); 
        b = coeffs(4:5);
        c = coeffs(6);
        x = -0.5 * (H\b); 
        gt(cnt,1) = x'*H*x + b'*x + c;
% Error control 
        if flag 
            err2 = abs(xerr'*H*xerr + b'*xerr + c - data(mi+1, ni+1 ));
            mult = 1; 
            if err2 < 1/20*abs(gt(cnt))
                mult = 1.25;
            elseif err2 > (5*abs(gt(cnt))+1e-5)
                mult = 0.75;
            end

            pky = quad2(data(mi-1:mi+1,ni));
            f = abs((pky + 1.1e-5) /(factor *(data(mi,ni)-pky) + 1e-5));
            d(1,2) = min(d(1,2), mult * ((f>=0.3) + (0.7 + f)*(f< 0.3) +(f>1)*(log(f)/100)));
            pkx =  quad2(data(mi, ni-1:ni+1)');
            f = abs((pkx + 1.1e-5) /(factor  *(data(mi,ni) - pkx) + 1e-5));
            d(1,1) = min(d(1,1), mult * ((f>=0.3) + (0.7 + f)*(f< 0.3) +(f>1)*(log(f)/100)));
        end
    elseif mi > 1 & mi < m
        gt(cnt,1) = quad2(data(mi-1:mi+1,ni));
        if flag 
            pky = gt(cnt);
            f = abs((pky + 1.1e-5) /(factor *(data(mi,ni) - pky) + 1e-5));
            d(1,2) = min(d(1,2), (f>=0.3) + (0.7 + f)*(f< 0.3) +(f>1)*(log(f)/100));
            if (pky > -1e-5), d(1,1) = min(1, d(1,1)); end
        end
    elseif ni > 1 & ni < n
        gt(cnt,1) = quad2(data(mi, ni-1:ni+1)');
        if flag 
            pkx = gt(cnt); 
            f = abs((pkx + 1.1e-5) /(factor  *(data(mi,ni) - pkx) + 1e-5));
            d(1,1) = min(d(1,1), (f>=0.3) + (0.7 + f)*(f< 0.3) +(f>1)*(log(f)/100));
            if (pkx > -1e-5), d(1,2) = min(1, d(1,2)); end
        end
    else
        gt(cnt,1) = data(mi, ni);
        if (gt(cnt) > -1e-5), 
            d(1,1) = min(1, d(1,1)); 
            d(1,2) = min(1, d(1,2)); 
        end
    end
end

% It is possible to have located no maxima at this point. If a non-empty
% data matrix has been specified, the data does have at least one maximum.
% We'll just return the raw maximum of the data elements. Note that NaN or
% inf will be allowed.
if isempty(gt) && ~isempty(data) 
    gt = max(data(:));
    % Use same code to update grid spacing as used if the maximum is
    % located in a corner of the data grid.
    if gt > -1e-5
        d(1,1) = min(1, d(1,1));
        d(1,2) = min(1, d(1,2));
    end
end

%--------------------------------------------------------------------------

function maximum=quad2(pts)
%QUAD2  Quadratically interpolates three points to find the maximum value.
%   QUAD2(PTS) interpolates the 3 element vector PTS to find the maximum.

c=pts(1);
ab=[-1 0.5; 2 -0.5]*(pts(2:3)-c*ones(2,1));
stepmin=-ab(2)/(2*ab(1));
maximum=ab(1)*stepmin^2+ab(2)*stepmin+c;

%--------------------------------------------------------------------------

function [maximum,err]=cubic(pts,checkpt,location)
%CUBIC Cubicly interpolates four points to find the maximum value.
%   The second argument is for estimation of the error in the 
%   interpolated maximum. 

d=pts(1);
abc=[0.5 -0.5 1/6 ; -2.5 2 -0.5; 3 -1.5 1/3]*[pts(2:4)-d*ones(3,1)];
root=real(sqrt(4*(abc(2)^2)-12*abc(1)*abc(3)));
x1=(-2*abc(2)+root)/(6*abc(1));
if 6*abc(1)*x1+2*abc(2)<0
    stepmin=x1; 
   else
    stepmin=(-2*abc(2)-root)/(6*abc(1));
end
maximum=abc(1)*stepmin^3+abc(2)*stepmin^2+abc(3)*stepmin+d;
if nargin>1
    if location==0
        checkpt2=-abc(1)+abc(2)-abc(3)+d;
    else
        checkpt2=64*abc(1)+16*abc(2)+4*abc(3)+d;
    end
    err=abs(checkpt-checkpt2);
end

%======================== END OF SUBFUNCTIONS =============================