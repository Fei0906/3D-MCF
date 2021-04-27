function [ph,weight,closure,datatype]=LoadData(INPUTFILE,DATATYPE,WRAPPEDPHASES, ...
    WEIGHTCOEFFICIENTS, INTERFEROGRAMS,COORDINATE)
%=======================================================================
%   Load data and check whether the input data is correct.
%
%   Fei Liu, Oct 2019
%=======================================================================


%   Load input file
msg=['The input file is ' INPUTFILE];
disp(msg);

load(INPUTFILE,WRAPPEDPHASES,WEIGHTCOEFFICIENTS,INTERFEROGRAMS,COORDINATE); 
eval(['ph=' WRAPPEDPHASES ';']);
eval(['weight=' WEIGHTCOEFFICIENTS ';']);
eval(['ifgs=' INTERFEROGRAMS ';']);
eval(['coor=' COORDINATE ';']);
datatype=DATATYPE;

ph=single(ph);
weight=single(weight);

%   Check input matrix
p1=0;
if DATATYPE=='D'
    [s1,p1]=size(ph);
    [s2,p2]=size(weight);
    [s3,p3]=size(coor);
    msg=['The size of input wrapped phases is '...
        num2str(s1) '*' num2str(p1) '.'  ];
    disp(msg);
    if s1~=s2 || p2~=p1
        error('Error occurred. The size of input weight matrix does not match phase matrix!');
    end
    if s3~=s2 || p3~=2
        error('Error occurred. The size of input coordinate matrix is not correct!');
    end
    
    if ~exist('TIN.mat')
        nodename=['uw.1.node'];
        fid=fopen(nodename,'w');
        fprintf(fid, '%d 2 0 0\n',s1);
        for i=1:s1
            fprintf(fid,'%d %f %f\n',i,coor(i,1),coor(i,2));
        end
        fclose(fid);

        !triangle -e uw.1.node

        fid=fopen('uw.2.edge','r');
        header=str2num(fgetl(fid));
        n_edge=header(1);
        edges_nz=zeros(n_edge,2);
        for i=1:n_edge
            linedata=str2num(fgetl(fid));
            edges_nz(i,:)=linedata(:,2:3);
        end
        fclose(fid);

        fid=fopen('uw.2.ele','r');
        header=str2num(fgetl(fid));
        n_ele=header(1);
        eles_nz=zeros(n_ele,6);
        for i=1:n_ele
            linedata=str2num(fgetl(fid));
            eles_nz(i,1:3)=linedata(:,2:4);
        end
        fclose(fid);

        for i=1:n_ele
            edgei=[eles_nz(i,1) eles_nz(i,2);eles_nz(i,2) eles_nz(i,3);eles_nz(i,3) eles_nz(i,1)];
            for j=1:3
                index1=find(edges_nz(:,1)==edgei(j,1));
                index2=find(edges_nz(:,2)==edgei(j,2));
                index=intersect(index1,index2);
                if ~isempty(index)
                    eles_nz(i,3+j)=index;
                else
                    index1=find(edges_nz(:,1)==edgei(j,2));
                    index2=find(edges_nz(:,2)==edgei(j,1));
                    index=intersect(index1,index2);
                    eles_nz(i,3+j)=-index;
                end
            end
        end
        save('TIN.mat','edges_nz','eles_nz');  
    end
else
    [m1,n1,p1]=size(ph);
    [m2,n2,p2]=size(weight);
    msg=['The size of input wrapped phases is '...
        num2str(m1) '*'  num2str(n1) '*' num2str(p1) '.'  ];
    disp(msg);
    
    if m2~=m1 || n2~=n1 || p2~=p1
        error('Error occurred. The size of input weight matrix does not match phase matrix!');
    end
end

if p1<3
	error('Error occurred. There should be at least three interferograms.');
end

%   Extract closure phases from ifgs.
[m3,n3]=size(ifgs);
if m3~=p1 || n3~=2
    error('Error occurred. The size of input interferogram matrix is not correct!');
end
if ~isequal(logical(~rem(ifgs,1)),ones(m3,n3)) || ~isempty(find(ifgs<=0, 1))
    error('Error occurred. The elements of interferogram matrix should all be positive integer!');
end

closure=cell(1,3*m3);
idx=1;
NK=nchoosek(1:m3,3);
for i=1:size(NK,1)
    un=union(union(ifgs(NK(i,1),:),ifgs(NK(i,2),:)),ifgs(NK(i,3),:));
    if size(un,2)==3
        a=ifgs(NK(i,1),1)-ifgs(NK(i,1),2);
        b=ifgs(NK(i,2),1)-ifgs(NK(i,2),2);
        c=ifgs(NK(i,3),1)-ifgs(NK(i,3),2);
        if a+b+c==0
            closure{idx}=[NK(i,1) NK(i,2) NK(i,3)];
        elseif a+b-c==0
            closure{idx}=[NK(i,1) NK(i,2) -NK(i,3)];
        elseif a-b+c==0
            closure{idx}=[NK(i,1) -NK(i,2) NK(i,3)];
        elseif a-b-c==0
            closure{idx}=[NK(i,1) -NK(i,2) -NK(i,3)];
        end
        idx=idx+1;
    end
end
closure=closure(~cellfun('isempty',closure));

end
