clear;

lon1=-68.75;
lon2=-67.5;
lat1=-22.5;
lat2=-24.25;
polygon=[lon1 lon1 lon2 lon2;lat1 lat2 lat2 lat1];

tile_ifg=15;

uph_data=load('phuw_sb2.mat');
aps_data=load('tca_sb2.mat');
ifg_data=load('ps2.mat');

idx=inpolygon(ifg_data.lonlat(:,1),ifg_data.lonlat(:,2),polygon(1,:),polygon(2,:));
lon=ifg_data.lonlat(idx,1);
lat=ifg_data.lonlat(idx,2);
ij=ifg_data.ij(idx,2:3);
coor=ifg_data.xy(idx,2:3);
    
uph=uph_data.ph_uw(idx,:);
aps=aps_data.ph_tropo_gacos(idx,:);
ifg_ix=ifg_data.ifgday_ix;

clear uph_data aps_data ifg_data;

uph=uph-aps;
ph=uph-round(uph/2/pi)*2*pi;
[n_pixel,n_ifg]=size(uph);

%   Read lines and samples
parlist=dir('Export/rslc/*.par');
fid=fopen(['Export/rslc/' parlist(1).name],'r');
par=textscan(fid,'%s%s%s%s');
fclose(fid);
samples=str2num(cell2mat(par{2}(15)));
lines=str2num(cell2mat(par{2}(16)));

%   Get coherence data as weight
ifglist=dir('Export/SMALL_BASELINES');
weight=zeros(size(ph));
%   Skip the first two '.' and '..' floder
for i=1:n_ifg
    fid=fopen(['./Export/SMALL_BASELINES/' ifglist(i+2).name ...
        '/' ifglist(i+2).name '.cc'],'r');
    cc=fread(fid,[samples,lines],'float',0,'b');
    
    for j=1:n_pixel
        weight(j,i)=cc(ij(j,2),ij(j,1));
    end
end

mkdir ThreeD_MCF
cd ThreeD_MCF

idx=1:tile_ifg:n_ifg;
for i=1:size(idx,2)
    eval(['mkdir tile_' num2str(i)]);
    eval(['cd tile_' num2str(i)]);
    if i==size(idx,2)
        wPh=ph(:,idx(i):end);
        wgt=weight(:,idx(i):end);
        ifgs=ifg_ix(idx(i):end,:);
        save('input_3DMCF_data','wPh','wgt','coor','ifgs');
        break;
    end
    
    wPh=ph(:,idx(i):idx(i+1));
    wgt=weight(:,idx(i):idx(i+1));
    ifgs=ifg_ix(idx(i):idx(i+1),:);
    save('input_3DMCF_data','wPh','wgt','coor','ifgs');
    
    cd ..;
end

cd ..;


