function Mesh=ele_cent_crdn(Mesh)
if length(Mesh.Elements) >= 4 && ~isempty(Mesh.Elements{4})
    tets = Mesh.Elements{4}; % 구조: [요소ID, Node1, Node2, Node3, Node4]
    
    % 노드 ID 분리
    n1 = tets(:, 2);
    n2 = tets(:, 3);
    n3 = tets(:, 4);
    n4 = tets(:, 5);
    
    % 💡 [안전장치] 노드 ID가 중간에 비어있거나 1부터 시작하지 않을 경우를 대비해
    % 노드 ID -> 배열의 행(Row) 인덱스로 즉시 변환하는 매핑 테이블 생성
    maxNodeID = numel(Mesh.Nodes(:,1));
    id2row = zeros(maxNodeID, 1);
    id2row(1:maxNodeID) = 1:size(Mesh.Nodes, 1);
    
    % 벡터화 연산으로 X, Y, Z 중심 좌표를 한 방에 계산 (C++의 for문 대체)
    cx = (Mesh.Nodes(id2row(n1), 1) + Mesh.Nodes(id2row(n2), 1) + ...
          Mesh.Nodes(id2row(n3), 1) + Mesh.Nodes(id2row(n4), 1)) / 4.0;
          
    cy = (Mesh.Nodes(id2row(n1), 2) + Mesh.Nodes(id2row(n2), 2) + ...
          Mesh.Nodes(id2row(n3), 2) + Mesh.Nodes(id2row(n4), 2)) / 4.0;
          
    cz = (Mesh.Nodes(id2row(n1), 3) + Mesh.Nodes(id2row(n2), 3) + ...
          Mesh.Nodes(id2row(n3), 3) + Mesh.Nodes(id2row(n4), 3)) / 4.0;
    
    % [N x 3] 행렬로 병합
    Mesh.ele_cent_crdn = [cx, cy, cz];
    
    % 텍스트 파일로 추출
    writematrix(Mesh.ele_cent_crdn, 'ele_cent_crdn.txt', 'Delimiter', 'space');
    fprintf('✅ [Saved] ele_cent_crdn.txt (%d개 사면체 요소의 중심 좌표)\n', size(tets, 1));
else
    disp('⚠️ 사면체(Type 4) 요소가 존재하지 않아 중심 좌표를 추출할 수 없습니다.');
end