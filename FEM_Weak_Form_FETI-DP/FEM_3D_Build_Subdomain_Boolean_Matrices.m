function SD=FEM_3D_Build_Subdomain_Boolean_Matrices(SD)

N_sd=length(SD);


rem_key_all=[];
rem_sd_all=[];
rem_k_all=[];
cor_key_all=[];
cor_sd_all=[];
cor_k_all=[];
%각 서브도메인의 edge_remainder 들을 가져오고(subdomain 내에서 local indexing 되어있음),
%각 서브도메인의 edg_nod_global을 가져와서 global하게 각 서브도메인의 remainder edge index를 비교함
%더 작은 인덱스의 서브도메인에 존재하는 remainder edge에 +1을 부여.
%더 큰 인덱스의 서브도메인에 존재하는 remainder edge에는 -1을 부여.
for s = 1 :N_sd

    r_edges=SD(s).edg_remainder(:);

    for k = 1:length(r_edges)
        le=r_edges(k);
        key=sort(SD(s).edg_nod_g(le,:));

        rem_key_all = [rem_key_all; key];
        rem_sd_all = [rem_sd_all; s];

        %k는 해당 서브도메인에서 몇번째 remainder edge인지 알려줌.
        rem_k_all = [rem_k_all; k];
    end
    
    % Bc를 만들기 위해 global corner edge들을 모음.
    c_edges=SD(s).edg_corner(:);

    for k = 1:length(c_edges)
        ce=c_edges(k);
        key_c=sort(SD(s).edg_nod_g(ce,:));
        cor_key_all=[cor_key_all;key_c];
        cor_sd_all=[cor_sd_all;s];
        cor_k_all=[cor_k_all;k];
    end
end

[~,~,rem_ic]=unique(rem_key_all,'rows');

%subdomain 내의 local remainder index가 global lambda matrix중 어디로 들어갈지 정해야함.
%즉 subdomain의 local remainder dof가 global remainder dof 중 몇번째인지 체크

N_lambda=max(rem_ic);

for s = 1:N_sd
    SD(s).lambda_idx=zeros(SD(s).ner,1);
end

%rem_k_all=[1번 서브도메인의 remainder ; 2번 서브도메인의 remainder ...]
for q=1:length(rem_k_all)
    s=rem_sd_all(q);

    %현재 다루는 dof가 서브도메인 내부에서 몇번째 local dof인지.
    k=rem_k_all(q);
    
    % local dof 자리에 global index를 집어넣음.
    % 예를 들어 local 32번 remainder dof 라면 k=32, rem_ic(q)는 이 q번이 몇번 global unique index로 매핑되는지를 나타냄.
    SD(s).lambda_idx(k)=rem_ic(q);
end
SD(1).N_lambda=N_lambda;

[unique_global_corner,~,cor_ic]=unique(cor_key_all,'rows');
N_global_corner=size(unique_global_corner,1);
for s = 1:N_sd
    SD(s).rem_sign=ones(SD(s).ner,1);
end

%unique한 전체 remainder edge 개수.
N_rem_group=max(rem_ic);

for g = 1:N_rem_group
    
    %remainder edge라면 2개의 edge가 unique한 edge 하나로 될것 
    %따라서 g와 같은 rem_ic는 항상 2개가 나올것.
    copies=find(rem_ic==g);

    if length(copies)~=2
        warning('Remainder edge group %d has %d copies. Expected 2.',g,length(copies));
        continue;
    end
    
    %c1과 c2는 각각 물리적으로 하나의 edge를 가리킴.
    c1=copies(1);
    c2=copies(2);
    
    %c1,c2가 어느 서브도메인에 속하는지
    s1=rem_sd_all(c1);
    s2=rem_sd_all(c2);
    
    %c1,c2가 각각 서브도메인에서 몇번째 local edge인지.
    k1=rem_k_all(c1);
    k2=rem_k_all(c2);

    part1= SD(s1).part_id;
    part2=SD(s2).part_id;

    if part1<part2
        SD(s1).rem_sign(k1)=+1;
        SD(s2).rem_sign(k2)=-1;
    else
        SD(s1).rem_sign(k1)=-1;
        SD(s2).rem_sign(k2)=+1;
    end

end


for s = 1:N_sd

    nei=SD(s).nei;
    ner=SD(s).ner;
    nec=SD(s).nec;

    n_local=nei+ner+nec;

    row_r=(1:ner).';
    col_r=nei+(1:ner).';
    val_r=SD(s).rem_sign(:);

    SD(s).Br=sparse(row_r,col_r,val_r,ner,nei+ner);

    %idx는 서브도메인 s의 모든 corner edge의 local nubmer (서브도메인내에서 edge index가 아님)
    % 예를 들어 서브도메인 1의 corner edge가 100개다 그러면 s=1일때 idx는 1:100 인셈.
    idx=find(cor_sd_all==s);

    row_c=cor_k_all(idx);
    col_c=cor_ic(idx);
    val_c=ones(length(idx),1);

    SD(s).Bc=sparse(row_c,col_c,val_c,nec,N_global_corner);
    
end

SD(1).unique_global_corner=unique_global_corner;
end

