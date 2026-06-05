function [N_0,N_1,N_2,N_3,DT,...
        ele_fac,ele_edg,ele_nod,ele_edg_sign,...
        fac_edg,fac_nod,...
        edg_nod,...
        edg_length,...
        nod_crdn]=FEM_3D_Mesh_Getter_gmsh(nod_crdn,ele_nod)

N_0=size(nod_crdn,1);
N_3=size(ele_nod,1);

DT=triangulation(ele_nod,nod_crdn);

loc_fac_nod=[1,2,3;1,2,4;1,3,4;2,3,4];
loc_edg_nod=[1,2;1,3;1,4;2,3;2,4;3,4];
loc_fac_edg=[1,2,4;1,3,5;2,3,6;4,5,6];

all_fac_key=zeros(4*N_3,3);
all_edg_key=zeros(6*N_3,2);

cnt_f=0;
cnt_e=0;

for e =1:N_3

    tet=ele_nod(e,:);
    
    local_faces=tet(loc_fac_nod); % 4 x 3
    local_edges=tet(loc_edg_nod); % 6 x 2

    all_fac_key(cnt_f+1:cnt_f+4,:)=sort(local_faces,2); %행 별로 오름차순으로 정렬
    all_edg_key(cnt_e+1:cnt_e+6,:)=sort(local_edges,2);

    cnt_f=cnt_f+4;
    cnt_e=cnt_e+6;
end

[fac_nod,~,fac_label]=unique(all_fac_key,'rows'); 
[edg_nod,~,edg_label]=unique(all_edg_key,'rows');

N_2=size(fac_nod,1);
N_1=size(edg_nod,1);

ele_fac = zeros(N_3,4);
ele_edg = zeros(N_3,6);
ele_edg_sign = zeros(N_3,6);

cnt_f = 0;
cnt_e = 0;

for e = 1:N_3

    tet = ele_nod(e,:);

    local_edges = tet(loc_edg_nod);

    for a = 1:6

        cnt_e = cnt_e + 1;

        gid = edg_label(cnt_e);

        ele_edg(e,a) = gid;

        global_edge = edg_nod(gid,:);
        local_edge  = local_edges(a,:);

        if isequal(local_edge,global_edge)
            ele_edg_sign(e,a) = 1;
        else
            ele_edg_sign(e,a) = -1;
        end

    end

    for f = 1:4

        cnt_f = cnt_f + 1;

        ele_fac(e,f) = fac_label(cnt_f);

    end

end

fac_edg=zeros(N_2,3);

for e = 1:N_3
    for f=1:4
        fid=ele_fac(e,f);

        fac_edg(fid,:)=ele_edg(e,loc_fac_edg(f,:));
    end
end
edg_length=zeros(N_1,1);

for i=1:N_1
    p1=nod_crdn(edg_nod(i,1),:);
    p2=nod_crdn(edg_nod(i,2),:);
    edg_length(i)=norm(p2-p1);
end
end
