function Mesh = Extract_PEC(Mesh)
    disp('Extracting PEC nodes and edges...');
    Mesh.PEC = struct('ele', [], 'nod', [], 'edg', []);
    
    % 1. PEC 물리 태그 ID 찾기
    physNames = Mesh.Phys.NameToID.keys();
    pec_id = [];
    for i = 1:length(physNames)
        if strcmpi(physNames{i}, 'pec')
            pec_id = Mesh.Phys.NameToID(physNames{i});
            break;
        end
    end
    
    if isempty(pec_id)
        disp('Info: Mesh 내에 PEC 영역이 존재하지 않습니다.');
        return;
    end
    
    % PEC 엘리먼트들 찾기
    is_pec_ele = (Mesh.ele_physTag(:, 2) == pec_id);
    pec_ele_ids = Mesh.ele_physTag(is_pec_ele, 1);
    Mesh.PEC.ele = sort(pec_ele_ids);
    
    pec_nodes = [];
    pec_edges = [];
    
    % ---------------------------------------------------------
    % 케이스 A: PEC가 표면(Triangle)으로 정의된 경우
    % ---------------------------------------------------------
    if ~isempty(Mesh.Elements{2})
        tri_data = Mesh.Elements{2};
        [is_pec_tri, ~] = ismember(tri_data(:,1), pec_ele_ids);
        
        if any(is_pec_tri)
            pec_tris = tri_data(is_pec_tri, 2:4);
            pec_nodes = [pec_nodes; pec_tris(:)];
            
            % 삼각형의 3개 엣지 조합
            tri_edges = [pec_tris(:,[1,2]); pec_tris(:,[2,3]); pec_tris(:,[1,3])];
            tri_edges = sort(tri_edges, 2);
            
            % 전역 엣지 리스트와 비교하여 매칭되는 인덱스 추출
            [~, edge_ids] = ismember(tri_edges, Mesh.Topology.edg_nod, 'rows');
            pec_edges = [pec_edges; edge_ids(edge_ids > 0)];
        end
    end
    
    % ---------------------------------------------------------
    % 케이스 B: PEC가 체적(Tetrahedron)으로 정의된 경우
    % ---------------------------------------------------------
    if ~isempty(Mesh.Elements{4})
        tet_data = Mesh.Elements{4};
        [is_pec_tet, ~] = ismember(tet_data(:,1), pec_ele_ids);
        
        if any(is_pec_tet)
            pec_tets = tet_data(is_pec_tet, 2:5);
            pec_nodes = [pec_nodes; pec_tets(:)];
            
            % 🚀 꿀팁: ele_edg는 사면체 인덱스와 1:1 대응하므로 그대로 복사!
            pec_tet_edges = Mesh.Topology.ele_edg(is_pec_tet, :);
            pec_edges = [pec_edges; pec_tet_edges(:)];
        end
    end
    
    % 최종적으로 중복 제거
    Mesh.PEC.nod = unique(pec_nodes);
    Mesh.PEC.edg = unique(pec_edges);
    
    disp(['[PEC Result] Nodes: ', num2str(length(Mesh.PEC.nod)), ...
          ' / Edges: ', num2str(length(Mesh.PEC.edg))]);
end