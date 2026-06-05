function [N_0,N_1,N_2,N_3,DT,...
    ele_fac,ele_edg,ele_nod,ele_edg_sign,...
    fac_edg,fac_nod,...
    edg_nod,...
    edg_length,...
    nod_crdn]=FEM_3D_Global_Mesh_Getter_gmsh(nod_crdn,ele_nod)

N_0=size(nod_crdn,1);
N_3=size(ele_nod,1);
DT=triangulation(ele_nod,nod_crdn);

fac_nod=zeros(4*N_3,3);
edg_nod=zeros(6*N_3,2);

loc_fac_nod=[1,2,3 ; 1,2,4 ; 1,3,4 ; 2,3,4];
loc_edg_nod=[1,2 ; 1,3 ; 1,4 ; 2,3 ; 2,4 ; 3,4];


cnt_f=0;
cnt_e=0;

all_fac_key=zeros(4*N_3,3);
all_edg_key=zeros(6*N_3,2);

for i = 1 :N_3
    
    tet=ele_nod(i,:); %원본 ele_nod를 복사
    
    local_faces=tet(loc_fac_nod); % 각 element를 이루는 face를 구성하는 노드를 우리 순서대로 정렬
    % 총 4N_3 개의 행이 생성
    % 원본 local 방향을 가진 face
    local_edges=tet(loc_edg_nod); % 각 element를 이루는 edge를 구성하는 노드를 우리 순서대로 정렬
    % 총 6N_3 개의 행이 생성
    % 원본 local 방향을 가진 edge

    all_fac_key(cnt_f+1:cnt_f+4,:)=sort(local_faces,2); % 행렬의 각 행 내부에서 정렬.
    %중복된 face들을 제거하기 위해서 행마다 정렬해줌.
    %정렬이 되었기 때문에 방향성은 소실.

    all_edg_key(cnt_e+1:cnt_e+6,:)=sort(local_edges,2);
    
    
    cnt_f=cnt_f+4;
    cnt_e=cnt_e+6;

end

[fac_nod,~,fac_label]=unique(all_fac_key,'rows'); 
%fac_nod는 all_fac_key에서 중복을 제거하고 정렬함. -> 원래의 local face 방향성, global face 행 순서도 사라짐
%fac_label에는 all_fac_key의 각 행이 global face list의 몇 번째인지 알려줌.

[edg_nod,~,edg_label]=unique(all_edg_key,'rows');

N_1=size(edg_nod,1);
N_2=size(fac_nod,1);

ele_fac=zeros(N_3,4);
ele_edg=zeros(N_3,6);
ele_edg_sign=zeros(N_3,6);

cnt_e=0;
cnt_f=0;

for i = 1 :N_3
    tet=ele_nod(i,:); %원본 element
    local_edges=tet(loc_edg_nod); %원본 local 방향성을 가진 edge. 6N_3 개의 행을 가짐
    
    for a=1:6
        cnt_e=cnt_e+1;
        gid=edg_label(cnt_e);
        %i번 element의 local 1번 edge부터 6번 edge까지 global edge list의 몇번째인지 확인.

        ele_edg(i,a)=gid; %i번 element의 local a번 edge에 찾아온 global edge 집어넣음
        %즉 element 내부의 local edge 순서는 보존

        global_edge=edg_nod(gid,:); % 방향성 없음.
        local_edge= local_edges(a,:); % 원본의 방향성이 존재

        if isequal(local_edge,global_edge) % 방향성이 같으면 +1, 다르면 -1
            ele_edg_sign(i,a)=1;
        else
            ele_edg_sign(i,a)=-1;
        end
    end

    for f=1:4
        cnt_f=cnt_f+1;
        ele_fac(i,f)=fac_label(cnt_f);
        %element 내부의 local face 순서는 보존하면서 
        %local f번째 element가 global face 몇번인지 적음.
    end
end
% -> ele_xxx 관계성은 1번 element부터 차례로 정렬

loc_fac_edg=[1,2,4 ; 1,3,5 ; 2,3,6 ; 4,5,6];

fac_edg=zeros(N_2,3);

for i = 1:N_3
    for f=1:4
        fid=ele_fac(i,f); %fid (face id)는 i번째 element의 f번째 face
        fac_edg(fid,:)=ele_edg(i,loc_fac_edg(f,:)); 
        %fid 번째 face를 구성하는 edge는 i번째 element를 구성하는 edge를 우리의 순서대로 정렬
        
    end
end
%-> fac_edg는 1번 face부터 face를 구성하는 edge로 정렬

edg_length = zeros(N_1,1);

for i = 1:N_1
    p1 = nod_crdn(edg_nod(i,1),:);
    p2 = nod_crdn(edg_nod(i,2),:);
    edg_length(i) = norm(p2-p1);
end

end
