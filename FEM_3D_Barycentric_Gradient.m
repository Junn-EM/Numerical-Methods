function [ele_vol, grad_lambda,coef_lambda,curl_W]=FEM_3D_Barycentric_Gradient(nod_crdn,ele_nod,N_3)

ele_vol=zeros(N_3,1);
grad_lambda=zeros(N_3,4,3);
coef_lambda=zeros(N_3,4,4);
curl_W=zeros(N_3,6,3);

loc_edg_nod=[1,2;1,3;1,4;2,3;2,4;3,4];
for i = 1:N_3
    nodes=ele_nod(i,:); %i번 element를 구성하는 노드 index
    X=nod_crdn(nodes,:); %i번 element를 구성하는 노드들의 좌표 4x3 [x1 y1 z1
                                                                 % x2 y2 z2
                                                                 % x3 y3 z3
                                                                 % x4 y4 z4]

    A=[ones(4,1),X]; %[1 x1 y1 z1  [a1   [1
                     % 1 x2 y2 z2   b1    0  
                     % 1 x3 y3 z3   c1 =  0   이라는 조건을 이용하여 각 노드에서의 barycentric coordinate
                     % 1 x4 y4 z4 ] d1]   0]  의 coefficient 찾기.
    
    detA=det(A);
    if abs(detA)<1e-14
        error('Degenerate tetrahedron detected at element %d. det(A)=%.3e',k,detA);
    end

    coef = A\eye(4);
    coef_lambda(i,:,:)=coef; %i번째 element에 4 x 4 [a1 a2 a3 a4
                                                   % b1 b2 b3 b4
                                                   % c1 c2 c3 c4
                                                   % d1 d2 d3 d4] 가 들어감

    grad=coef(2:4,:).'; %[b1 c1 d1      [grad lambda1
                        % b2 c2 d2       grad lambda2  
                        % b3 c3 d3    =  grad lambda3
                        % b4 c4 d4]      grad lambda4]
    grad_lambda(i,:,:)=grad;
    
    J=[                %기준노드(1) 로부터 다른 노드까지의 벡터
        X(2,:)-X(1,:);
        X(3,:)-X(1,:);  
        X(4,:)-X(1,:)
        ];

    ele_vol(i)=abs(det(J))/6; %각 element의 부피는 각 벡터로 이루어진 행렬의 determinant/6

    for a=1:6
        k=loc_edg_nod(a,1);
        j=loc_edg_nod(a,2);

        curl_W(i,a,:)=2*cross(grad(k,:),grad(j,:));
    end
end
end
