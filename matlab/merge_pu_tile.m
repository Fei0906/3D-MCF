clear;

lon1=-70;
lon2=-66;
lat1=-20;
lat2=-26;
polygon=[lon1 lon1 lon2 lon2;lat1 lat2 lat2 lat1];

n_tile=43;
n_ifg=516;
n_pixel=148947;
dratio=10;

ifgs=zeros(n_ifg,2);
uph=zeros(n_pixel,n_ifg);

for i=1:n_tile
    eval(['cd tile_' num2str(i)]);
    input=load('input_3DMCF_data.mat');
    output=load('uph.mat');
    coor=input.coor;
    
    n_ifgi=size(input.ifgs,1);
    if i==1
        ifgs(1:n_ifgi,:)=input.ifgs;
        uph(:,1:n_ifgi)=output.uph;
    else
        for j=1:n_ifg
            if ifgs(j,:)==input.ifgs(1,:)
                ifgs(j:j+n_ifgi-1,:)=input.ifgs;
                uph(:,j:j+n_ifgi-1)=output.uph;
                break;
            end
        end
    end
    
    cd ..;
end

%load('phuw_3d','ifgs','uph','coor');

uph_data=load('../phuw_sb2.mat');
aps_data=load('../tca_sb2.mat');
ifg_data=load('../ps2.mat');

idx=inpolygon(ifg_data.lonlat(:,1),ifg_data.lonlat(:,2),polygon(1,:),polygon(2,:));

uph2=uph_data.ph_uw(idx,:);
aps=aps_data.ph_tropo_gacos(idx,:);
aps=aps(1:dratio:end,:);
uph2=uph2(1:dratio:end,:);

uph=uph+aps;
dph=uph-uph2;

dk=mode(round(dph/2/pi),1);
uph=uph-repmat(dk*2*pi,size(uph,1),1);

save('phuw_sb2_3d','uph','aps');
quit;


