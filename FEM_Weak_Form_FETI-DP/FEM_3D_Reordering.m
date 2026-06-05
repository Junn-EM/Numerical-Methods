function SD=FEM_3D_Reordering(N_sd,SD)

%% node
all_nod_key=[];
all_nod_sd=[];
all_nod_loc=[];

for s = 1:N_sd
    nod_key_s=SD(s).nod_g(:); %각 서브도메인을 구성하는 node들의 global index
    n_nod_s=SD(s).N_0;

    %각 서브도메인들의 node들의 global index 모음
    all_nod_key=[all_nod_key;nod_key_s];

    %각 노드들이 어떤 서브도메인에서 왔는지 기록
    all_nod_sd=[all_nod_sd;s*ones(n_nod_s,1)];

    %각 노드들이 해당 서브도메인의 몇번째 local node인지 기록
    all_nod_loc=[all_nod_loc;(1:n_nod_s).'];
end

[unique_nod_key,~,nod_ic]=unique(all_nod_key,'rows');

%같은 노드가 몇번이나 나왔는지 카운트
nod_count=accumarray(nod_ic,1);

%% edge
all_edg_key=[];
all_edg_sd=[];
all_edg_loc=[];

for s = 1:N_sd
    edg_key_s=sort(SD(s).edg_nod_g,2);
    n_edg_s=size(edg_key_s,1);

    all_edg_key = [all_edg_key; edg_key_s];
    all_edg_sd = [all_edg_sd; s * ones(n_edg_s, 1)];
    all_edg_loc = [all_edg_loc; (1:n_edg_s).'];
end

[unique_edg_key,~,edg_ic]=unique(all_edg_key,'rows');
edg_count=accumarray(edg_ic,1);

%% face
all_fac_key=[];
all_fac_sd=[];
all_fac_loc=[];

for s = 1:N_sd
   fac_key_s=sort(SD(s).fac_nod_g,2);
   n_fac_s=size(fac_key_s,1);

   all_fac_key=[all_fac_key;fac_key_s];
   all_fac_sd=[all_fac_sd;s*ones(n_fac_s,1)];
   all_fac_loc=[all_fac_loc;(1:n_fac_s).'];
end

[unique_fac_key,~,fac_ic]=unique(all_fac_key,'rows');
fac_count=accumarray(fac_ic,1);

%

for s = 1:N_sd
    SD(s).nod_share_count=zeros(SD(s).N_0,1);
    SD(s).edg_share_count=zeros(SD(s).N_1,1);
    SD(s).fac_share_count=zeros(SD(s).N_2,1);
end

%모든 서브도메인의 노드 개수만큼 돌면서 
for k = 1:length(all_nod_loc)
    s=all_nod_sd(k); %현재 몇번 서브도메인인지
    ln=all_nod_loc(k); % 현재 서브도메인의 몇번째 로컬 노드 인지

    %nod_ic(k) -> 현재 로컬 인덱스를 갖는 노드가 unique의 몇번인지
    %nod_count(unique index) -> 해당 unique index가 몇번 count 되었는지 나타냄.
    SD(s).nod_share_count(ln)=nod_count(nod_ic(k)); 
end

for k = 1:length(all_edg_loc)
    s=all_edg_sd(k);
    le=all_edg_loc(k);

    SD(s).edg_share_count(le)=edg_count(edg_ic(k));
end

for k = 1:length(all_fac_loc)

    s=all_fac_sd(k);
    lf=all_fac_loc(k);

    SD(s).fac_share_count(lf)=fac_count(fac_ic(k));
end

%% Classify node/edge/face to inner,remainder ,corner
for s = 1:N_sd

    nc=SD(s).nod_share_count;
    SD(s).nod_inner = find(nc==1);
    SD(s).nod_remainder=find(nc==2);
    SD(s).nod_corner=find(nc>=3);

    SD(s).nod_inner=SD(s).nod_inner(:);
    SD(s).nod_remainder=SD(s).nod_remainder(:);
    SD(s).nod_corner=SD(s).nod_corner(:);

    ec=SD(s).edg_share_count;

    SD(s).edg_inner=find(ec==1);
    SD(s).edg_remainder = find(ec == 2);
    SD(s).edg_corner = find(ec >= 3);
    
    SD(s).edg_inner=SD(s).edg_inner(:);
    SD(s).edg_remainder = SD(s).edg_remainder(:);
    SD(s).edg_corner = SD(s).edg_corner(:);

    fc=SD(s).fac_share_count;

    SD(s).fac_inner=find(fc==1);
    SD(s).fac_remainder = find(fc == 2);
    
    SD(s).fac_inner = SD(s).fac_inner(:);
    SD(s).fac_remainder = SD(s).fac_remainder(:);

    SD(s).node=[SD(s).nod_inner;SD(s).nod_remainder;SD(s).nod_corner];
    SD(s).edge=[SD(s).edg_inner;SD(s).edg_remainder;SD(s).edg_corner];
    SD(s).face=[SD(s).fac_inner;SD(s).fac_remainder];
    
    SD(s).nni=length(SD(s).nod_inner);
    SD(s).nnr=length(SD(s).nod_remainder);
    SD(s).nnc=length(SD(s).nod_corner);

    SD(s).nei=length(SD(s).edg_inner);
    SD(s).ner=length(SD(s).edg_remainder);
    SD(s).nec=length(SD(s).edg_corner);

    SD(s).nfi=length(SD(s).fac_inner);
    SD(s).nfr=length(SD(s).fac_remainder);


end

GlobalCount=struct();
GlobalCount.nod_key=unique_nod_key;
GlobalCount.nod_count=nod_count;

GlobalCount.edg_key=unique_edg_key;
GlobalCount.edg_count=edg_count;

GlobalCount.fac_key=unique_fac_key;
GlobalCount.fac_count=fac_count;

SD(1).GlobalCount=GlobalCount;

end