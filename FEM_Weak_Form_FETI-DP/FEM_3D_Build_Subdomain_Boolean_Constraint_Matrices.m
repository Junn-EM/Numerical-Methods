function [SD,N_lambda,N_primal]=FEM_3D_Build_Subdomain_Boolean_Constraint_Matrices(SD)

%% edge dof에 대해서만 진행중. node dof에 대한것은 추후에 필요시 추가.

%% Remainder edge
rem_key_all=[];
rem_sd_all=[];
rem_le_all=[];
rem_pos_all=[];
N_sd=length(SD);
for s = 1:N_sd

    %서브도메인의 remainder edge
    r_edges=SD(s).edg_remainder;

    for k = 1 : length(r_edges)

        %서브도메인의 remainder edge 하나하나
        le = r_edges(k);

        %서브도메인의 remainder edge를 구성하는 global nodes index
        edge_nodes=SD(s).edg_nod_g(le,:);

        key=sort(edge_nodes);
        
        local_pos=SD(s).nei+k;
        
        %모든 서브도메인의 remainder edge_nodes
        rem_key_all=[rem_key_all;key];

        %모든 서브도메인의 index
        rem_sd_all=[rem_sd_all;s];
        
        %모든 서브도메인의 remainder edge
        rem_le_all=[rem_le_all;le];

        %모든 서브도메인의 local remainder edge index
        %즉 현재 서브도메인에서 edge는 [inner;remainder;corner] 순인데 여기서 몇번째인지.
        rem_pos_all=[rem_pos_all;local_pos];
    end
end

[unique_rem_key,~,rem_ic]=unique(rem_key_all,'rows');

%전체 unique remainder edge 개수
N_rem_group=size(unique_rem_key,1);

Ir=cell(N_sd,1);
Jr=cell(N_sd,1);
Vr=cell(N_sd,1);

for s= 1:N_sd
    Ir{s}=[];
    Jr{s}=[];
    Vr{s}=[];
end

row_id=0;

for g= 1:N_rem_group
    
    %g는 1부터 차례대로 올라가는중, g 값을 갖는 rem_ic행들을 찾음. 
    %즉 전체 remainder 목록중에서 같은 remainder edge를 가지는 행을 반환.
    %만약 copies가 [12,5234,38290] 이라면 이 세개의 행은 물리적으로 같은 edge를 가리킴.
    copies=find(rem_ic==g);
    
    % copies가 한개면 remainder도 아니고, corner도 아님
    if length(copies)<2
        continue;
    end
    
    ref=copies(1);
    
    for kk=2:length(copies)

        row_id=row_id+1;

        %copies의 1번을 reference로 쓰고 그 다음것부터 비교
        c1=ref;
        c2=copies(kk);
        
        %copies 1번이 속한 서브도메인
        s1=rem_sd_all(c1);

        %copies 1번이 속한 서브도메인에서의 local remainder edge index
        %즉 해당 서브도메인에서 몇번째 edge인지
        pos1=rem_pos_all(c1);

        s2=rem_sd_all(c2);
        pos2=rem_pos_all(c2);

        Ir{s1}=[Ir{s1};row_id];
        Jr{s1}=[Jr{s1};pos1];
        Vr{s1}=[Vr{s1};1];

        Ir{s2}=[Ir{s2};row_id];
        Jr{s2}=[Jr{s2};pos2];
        Vr{s2}=[Vr{s2};-1];
    end
end

N_lambda=row_id;

for s = 1:N_sd
    n_local=SD(s).nei+SD(s).ner+SD(s).nec;
    SD(s).Br=sparse(Ir{s},Jr{s},Vr{s},N_lambda,n_local);
end

N_primal=0;
end

%% Corner Edge

cor_key_all=[];
cor_sd_all=[];
cor_le_all=[];
cor_pos_all=[];

for s =1:N_sd
    c_edges=SD(s).edg_corner;

    for k = 1:length(c_edges)

        le=c_edges(k);
        edge_nodes=SD(s).edg_nod_g(le,:);
        key=sort(edge_nodes);
        local_pos=SD(s).nei+SD(s).ner+k;

        cor_key_all=[cor_key_all;key];
        cor_sd_all=[cor_sd_all;s];
        cor_le_all=[cor_le_all;le];
        cor_pos_all=[cor_pos_all;local_pos];

    end
end
[unique_cor_key,~,cor_ic]=unique(cor_key_all,'rows');
N_primal=size(unique_cor_key,1);

for s= 1:N_sd

