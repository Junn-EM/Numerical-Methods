function Mesh = Build_Topology(Mesh)
    disp('Building Mesh Topology (Edges & Faces)...');
    
    % 노드 좌표 (Mesh.Nodes의 행 번호가 곧 노드 ID라고 가정)
    nod_crdn = Mesh.Nodes; 
    
    % 사면체(Type 4) 데이터 추출
    if isempty(Mesh.Elements{4})
        error('사면체(Tetrahedron) 엘리먼트가 존재하지 않습니다.');
    end
    ele_nod = Mesh.Elements{4}(:, 2:5);
    N_3 = size(ele_nod, 1);
    
    fac_nod = zeros(4*N_3, 3);
    edg_nod = zeros(6*N_3, 2);
    
    loc_fac_nod = [1,2,3; 1,2,4; 1,3,4; 2,3,4];
    loc_edg_nod = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
    
    cnt_1 = 0; cnt_2 = 0;
    for i = 1:N_3
        sorted_nodes = sort(ele_nod(i,:));        
        fac_nod(cnt_2+1:cnt_2+4, :) = sorted_nodes(loc_fac_nod);
        edg_nod(cnt_1+1:cnt_1+6, :) = sorted_nodes(loc_edg_nod);
        cnt_1 = cnt_1 + 6;        
        cnt_2 = cnt_2 + 4;    
    end
    
    % 고유한 면과 엣지 추출
    [fac_nod, ~, fac_label] = unique(fac_nod, 'rows');
    [edg_nod, ~, edg_label] = unique(edg_nod, 'rows');
    
    N_2 = size(fac_nod, 1);
    N_1 = size(edg_nod, 1);
    
    ele_fac = zeros(N_3, 4);
    ele_edg = zeros(N_3, 6);
    
    cnt_1 = 0; cnt_2 = 0;
    for i = 1:N_3   
        ele_edg(i,:) = edg_label(cnt_1+1:cnt_1+6);
        ele_fac(i,:) = fac_label(cnt_2+1:cnt_2+4);    
        cnt_1 = cnt_1 + 6;    
        cnt_2 = cnt_2 + 4;                    
    end
    
    % 🚀 [최적화] 엣지 길이를 for문 없이 완벽하게 벡터 연산으로 처리 (초고속)
    edg_crdn_1 = nod_crdn(edg_nod(:, 1), :);
    edg_crdn_2 = nod_crdn(edg_nod(:, 2), :);    
    edg_length = sqrt(sum((edg_crdn_1 - edg_crdn_2).^2, 2)); 
    
    % Mesh 구조체에 저장
    Mesh.Topology.edg_nod = edg_nod;
    Mesh.Topology.ele_edg = ele_edg;
    Mesh.Topology.fac_nod = fac_nod;
    Mesh.Topology.ele_fac = ele_fac;
    Mesh.Topology.edg_length = edg_length;
    
    disp(['Topology Built: ', num2str(N_1), ' Edges, ', num2str(N_2), ' Faces.']);
end