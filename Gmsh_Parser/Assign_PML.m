%%Natural Unit 전용 PML

function [pml_eps, pml_mu] = Assign_PML(target_ele_ids, Mesh, omega,bg_eps, bg_mu)
    num_targets = length(target_ele_ids);
    pml_eps = zeros(num_targets, 9);
    pml_mu  = zeros(num_targets, 9);
    
    % 1. 전체 도메인 Bounding Box 구하기
    % (Mesh.Nodes는 1열~3열이 각각 X, Y, Z 좌표라고 가정)
    all_coords = Mesh.Nodes; 
    X_min = min(all_coords(:,1)); X_max = max(all_coords(:,1));
    Y_min = min(all_coords(:,2)); Y_max = max(all_coords(:,2));
    Z_min = min(all_coords(:,3)); Z_max = max(all_coords(:,3));
    
    % PML 두께 설정 
    % (가장 좋은 건 GMSH의 pml_thickness 파라미터 값과 일치시키는 것. 여기선 근사치 사용)
    PML_thickness = (X_max - X_min) * 0.1; 
    
    % 2. 무차원(Unitless) 감쇠 팩터 계산 (단위계 충돌 방지 완벽 적용)
    m = 3; 
    R_0 = 1e-4; % 저주파수에서는 너무 작게 잡으면 행렬이 터지므로 1e-4 정도로 타협
    max_damping = -(m + 1) * log(R_0) / (2 * omega * PML_thickness); 
    
    % 3. 각 요소별 중심점 기준 PML 텐서 계산
    for i = 1:num_targets
        e_id = target_ele_ids(i);
        
        % --- 중심점(Centroid) 계산 ---
        % 현재 요소(사면체)를 구성하는 4개의 노드 ID를 가져옴
        node_idx = Mesh.Elements{4}(e_id, 2:5); 
        % 해당 4개 노드의 (X, Y, Z) 좌표들의 평균을 내어 중심점을 구함
        cent_coord = mean(Mesh.Nodes(node_idx, :), 1); 
        x = cent_coord(1); y = cent_coord(2); z = cent_coord(3);
        
        % --- 각 축 방향으로 PML 영역 깊이(d) 계산 (+방향, -방향 모두 고려) ---
        % X축
        if x > (X_max - PML_thickness)
            d_x = x - (X_max - PML_thickness);     % +X 방향 PML 진입
        elseif x < (X_min + PML_thickness)
            d_x = (X_min + PML_thickness) - x;     % -X 방향 PML 진입
        else
            d_x = 0;                               % 내부 영역
        end
        
        % Y축
        if y > (Y_max - PML_thickness)
            d_y = y - (Y_max - PML_thickness);
        elseif y < (Y_min + PML_thickness)
            d_y = (Y_min + PML_thickness) - y;
        else
            d_y = 0;
        end
        
        % Z축
        if z > (Z_max - PML_thickness)
            d_z = z - (Z_max - PML_thickness);
        elseif z < (Z_min + PML_thickness)
            d_z = (Z_min + PML_thickness) - z;
        else
            d_z = 0;
        end
        
        % --- Stretching variables 계산 ---
        % d_x, d_y, d_z가 0이면 허수부가 0이 되어 자연스럽게 1이 됨
        s_x = 1 - 1j * max_damping * (d_x / PML_thickness)^m;
        s_y = 1 - 1j * max_damping * (d_y / PML_thickness)^m;
        s_z = 1 - 1j * max_damping * (d_z / PML_thickness)^m;
        
        % --- 이방성(Anisotropic) 텐서 조립 ---
        % bg_eps(1), bg_eps(5), bg_eps(9)는 각각 배경 매질의 xx, yy, zz 성분
        eps_xx = bg_eps(1) * (s_y * s_z) / s_x;
        eps_yy = bg_eps(5) * (s_x * s_z) / s_y;
        eps_zz = bg_eps(9) * (s_x * s_y) / s_z;
        
        mu_xx  = bg_mu(1)  * (s_y * s_z) / s_x;
        mu_yy  = bg_mu(5)  * (s_x * s_z) / s_y;
        mu_zz  = bg_mu(9)  * (s_x * s_y) / s_z;
        
        pml_eps(i, :) = [eps_xx, 0, 0,  0, eps_yy, 0,  0, 0, eps_zz];
        pml_mu(i, :)  = [mu_xx,  0, 0,  0, mu_yy,  0,  0, 0, mu_zz];
    end
end