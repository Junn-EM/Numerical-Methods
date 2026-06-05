%%
delete(gcp('nocreate'));
%%
parpool('Threads',48);
clear; clc; close all;
tstart=tic;

%% Constants - Natural Unit
[eps_0,mu_0,eta_0,c]=Constants();
[eps_0_SI,mu_0_SI,eta_0_SI,c_SI]=Constants_SI();

%% Input (SI unit)
load frequency_SI.txt;
load scale.txt;
load omega_tilde.txt;
load nod_crdn.txt;
load ele_nod.txt;
load ele_cent_crdn.txt;
load ele_part.txt;
load pec_edg_list.txt;
load pec_nod_list.txt;
eps_r = readmatrix('eps_r.txt');
mu_r  = readmatrix('mu_r.txt');

%% Global Domain Mesh Information Getter
disp('Building Global Mesh Connectivity...')
tic
[N_0,N_1,N_2,N_3,DT,...
    ele_fac,ele_edg,ele_nod,ele_edg_sign,...
    fac_edg,fac_nod,...
    edg_nod,...
    edg_length,...
    nod_crdn]=FEM_3D_Global_Mesh_Getter_gmsh(nod_crdn,ele_nod);
toc

%% Global Source Vector
f=complex(zeros(N_1,1));
src_edg=1;
f(src_edg)=1;

%% Getting Global PEC List
pec_edg=unique(pec_edg_list(:));
pec_nod=unique(pec_nod_list(:));

%% Subdomain Mesh Information getter
disp('Building SD(s).ele from ele_part.txt ...')
tic
part_ids=unique(ele_part(:,2));
N_sd=length(part_ids);
SD=struct();

for s=1:N_sd
    part_id=part_ids(s);
    
    %ele_part의 2열(partition index)과 part_id를 비교해서 같은것의 1열(element id)만 가져와라
    ele_s=ele_part(ele_part(:,2)==part_id,1); 

    SD(s).part_id=part_id;
    SD(s).ele=ele_s(:);
    SD(s).N_3=length(ele_s);
end
toc

disp('Building Subdomain Connectivity ...')
tic

SD_cell=cell(N_sd,1);

parfor s = 1:N_sd
    S=SD(s);
    
    %각 서브도메인의 global element index
    ele_g=S.ele;
    
    %Global table에서 subdomain에 해당하는 것들만 가져다가 unique 정렬.
    nod_g=unique(ele_nod(ele_g,:));
    edg_g=unique(ele_edg(ele_g,:));
    fac_g=unique(ele_fac(ele_g,:));

    %global node -> subdomain local node index로 매핑
    g2l_nod=zeros(N_0,1);
    g2l_edg=zeros(N_1,1);
    g2l_fac=zeros(N_2,1);

    % 이제 g2l_nod에 해당 서브도메인의 global node index를 넣으면 그 서브도메인에서의 local index가 나옴
    g2l_nod(nod_g)=1:length(nod_g);
    g2l_edg(edg_g)=1:length(edg_g);
    g2l_fac(fac_g)=1:length(fac_g);

    S.nod_g=nod_g(:);
    S.edg_g=edg_g(:);
    S.fac_g=fac_g(:);
    
    %서브도메인의 ele_nod, ele_edg, ele_fac local connectivity.
    %즉 connectivity를 서브도메인의 local index 가지고 나타냄.
    S.ele_nod=g2l_nod(ele_nod(ele_g,:));
    S.ele_edg=g2l_edg(ele_edg(ele_g,:));
    S.ele_fac=g2l_fac(ele_fac(ele_g,:));
    
    S.ele_edg_sign=ele_edg_sign(ele_g,:);
    
    %나머지 global connectivity도 subdomain의 local index를 가지고 재조립.
    S.edg_nod=g2l_nod(edg_nod(edg_g,:));
    S.fac_nod=g2l_nod(fac_nod(fac_g,:));
    S.fac_edg=g2l_edg(fac_edg(fac_g,:));

    S.edg_nod_g=edg_nod(edg_g,:);
    S.fac_nod_g=fac_nod(fac_g,:);
    
    S.edg_length=edg_length(edg_g);
    S.nod_crdn=nod_crdn(nod_g,:);
    
    S.N_0=length(S.nod_g);
    S.N_1=length(S.edg_g);
    S.N_2=length(S.fac_g);
    SD_cell{s}=S;
    
end
SD=[SD_cell{:}];
toc

%% Reordering 
disp('Reordering DoF (inner, remainder, corner)...')
tic
SD=FEM_3D_Reordering(N_sd,SD);
toc

%% Getting Barycentric Information
disp('Building Barycentric information...')
tic
parfor s = 1:N_sd

    S=SD(s);
    [S.ele_vol,S.grad_lambda,S.coef_lambda,S.curl_W]=FEM_3D_Barycentric_Gradient(S.nod_crdn,S.ele_nod,S.N_3);

    SD_cell{s}=S;
end
toc
SD=[SD_cell{:}];

%% Quadrature Information
disp('Building Information of 4-Quadrature Points in a Element ...')
tic
[lambda_q,weight_q]=Tet_Quad_4point();
toc

%% Getting Whitney Forms
disp('Building Whitney Forms...')
tic
SD_cell=cell(N_sd,1);
parfor s = 1:N_sd
    S=SD(s);
%각 element의 4개의 quadrature 지점에서의
%W_ij(x_q)=lambda_i(x_q)*grad_lambda_j-lambda_j(x_q)*grad_lambda_i 계산
%whitney 1 form은 edge가 만든 element 내부의 vector field라고 이해하고 넘어가자.
    [S.whitney_1]=Whitney_Forms_Quad(S.N_3,S.grad_lambda,lambda_q);  
    SD_cell{s}=S;
end
SD=[SD_cell{:}];
toc

%% Subdomain Matrices
disp('Building Subdomain Matrices...')
tic

SD_cell=cell(N_sd,1);

parfor s= 1:N_sd
    S=SD(s);

    eps_r_s=eps_r(S.ele,:);
    mu_r_s=mu_r(S.ele,:);

    I_all=zeros(36,S.N_3);
    J_all=zeros(36,S.N_3);
    V_all=complex(zeros(36,S.N_3));
    
    

    for e = 1:S.N_3
        [A_loc,Meps_loc, Mmu_loc, Kcurl_loc, Ccurl_loc]=FEM_3D_Local_Matrices(e,S.ele_vol,weight_q,S.whitney_1,S.curl_W,eps_r_s,mu_r_s,omega_tilde);

    I_e=zeros(36,1);
    J_e=zeros(36,1);
    V_e=complex(zeros(36,1));

    cnt=0;
    for a =1:6
        ia=S.ele_edg(e,a); %현재 element e의 a번째 edge의 글로벌 index
        sa=S.ele_edg_sign(e,a);

        for b=1:6
            ib=S.ele_edg(e,b); % 현재 element e의 b번째 edge의 글로벌 index
            sb=S.ele_edg_sign(e,b);
            
            cnt=cnt+1;
            I_e(cnt)=ia; %각 element의 행,열 index와 거기에 들어갈 value
            J_e(cnt)=ib;
            V_e(cnt)=sa*sb*A_loc(a,b);
        end
    end

    
    I_all(:,e)=I_e;
    J_all(:,e)=J_e;
    V_all(:,e)=V_e;

    end

    K_old=sparse(I_all(:),J_all(:),V_all(:),S.N_1,S.N_1);

    %subdomain source vector
    f_old=f(S.edg_g);
    
    is_pec_edge=ismember(S.edg_g,pec_edg);
    is_pec_node=ismember(S.nod_g,pec_nod);

    free_inner_e=S.edg_inner(~is_pec_edge(S.edg_inner));
    free_remainder_e=S.edg_remainder(~is_pec_edge(S.edg_remainder));
    free_corner_e=S.edg_corner(~is_pec_edge(S.edg_corner));
    
    free_inner_n=S.nod_inner(~is_pec_node(S.nod_inner));
    free_remainder_n=S.nod_remainder(~is_pec_node(S.nod_remainder));
    free_corner_n=S.nod_corner(~is_pec_node(S.nod_corner));

    
    %조립된 subdomain matrix를 [inner,remainder,corner] 순으로 재정렬
    %아직은 dof가 edge 요소뿐이니 edge요소 가지고만 재조립
    %node dof가 필요할땐 node 도 추가
    p=[free_inner_e;free_remainder_e;free_corner_e];
    S.edge=p(:);
   
    S.edg_inner=free_inner_e(:);
    S.edg_remainder=free_remainder_e(:);
    S.edg_corner=free_corner_e(:);

    S.nei=length(S.edg_inner);
    S.ner=length(S.edg_remainder);
    S.nec=length(S.edg_corner);
    
    S.N_1_free=length(S.edge);

    K_s=K_old(p,p);
    f_s=f_old(p);
    
    S.K=K_s;
    S.f=f_s;

    nei=S.nei;
    ner=S.ner;
    nec=S.nec;

    idx_ei=1:nei;
    idx_er=nei+(1:ner);
    idx_ec=nei+ner+(1:nec);

    S.K_ii = K_s(idx_ei,idx_ei);
    S.K_ir = K_s(idx_ei,idx_er);
    S.K_ic = K_s(idx_ei,idx_ec);

    S.K_ri = K_s(idx_er,idx_ei);
    S.K_rr = K_s(idx_er,idx_er);
    S.K_rc = K_s(idx_er,idx_ec);

    S.K_ci = K_s(idx_ec,idx_ei);
    S.K_cr = K_s(idx_ec,idx_er);
    S.K_cc = K_s(idx_ec,idx_ec);
    
    S.f_i=f_s(idx_ei);
    S.f_r=f_s(idx_er);
    S.f_c=f_s(idx_ec);
    
    idx_eIR=1:(S.nei+S.ner);
    
    %계산용 matrix form으로 재조립.
    S.K_RR=S.K(idx_eIR,idx_eIR);
    S.K_RC=S.K(idx_eIR,idx_ec);
    S.K_CR=S.K(idx_ec,idx_eIR);
    S.K_CC=S.K(idx_ec,idx_ec);
    
    S.f_R=S.f(idx_eIR);
    S.f_C=S.f(idx_ec);
    SD_cell{s} = S;
end
SD=[SD_cell{:}];
toc

%% Boolean Matrices
disp('Building Subdomain Boolean Matrices...')
tic
SD=FEM_3D_Build_Subdomain_Boolean_Matrices(SD);
toc

%% FETI-DP Interface Problem
disp('Building FETI-DP Interface Problem Matrices...')
tic

N_lambda=SD(1).N_lambda;
N_primal=size(SD(1).unique_global_corner,1);

F_Irr=sparse(N_lambda,N_lambda);
F_Irc=sparse(N_lambda,N_primal);
F_Icr=sparse(N_primal,N_lambda);
Kcc_star=sparse(N_primal,N_primal);

d_r=complex(zeros(N_lambda,1));
f_c_star=complex(zeros(N_primal,1));

for s = 1:N_sd
    S = SD(s);

    Br = S.Br;       % ner x (nei+ner)
    Bc = S.Bc;       % nec x N_primal

    KRR = S.K_RR;
    KRC = S.K_RC;
    KCR = S.K_CR;
    KCC = S.K_CC;

    fR = S.f_R;
    fC = S.f_C;

    rows = S.lambda_idx(:);   % global lambda index for local remainder edges

    KRR_fac = decomposition(KRR,'lu');

    KRR_inv_BrT   = KRR_fac \ Br.';
    KRR_inv_KRCBc = KRR_fac \ (KRC * Bc);
    KRR_inv_fR    = KRR_fac \ fR;

    F_Irr_s = Br * KRR_inv_BrT;
    F_Irc_s = Br * KRR_inv_KRCBc;
    F_Icr_s = Bc.' * KCR * KRR_inv_BrT;

    Kcc_star_s = Bc.' * KCC * Bc ...
               - Bc.' * KCR * KRR_inv_KRCBc;

    d_r_s = Br * KRR_inv_fR;

    f_c_star_s = Bc.' * fC ...
               - Bc.' * KCR * KRR_inv_fR;

    F_Irr(rows,rows) = F_Irr(rows,rows) + F_Irr_s;
    F_Irc(rows,:)    = F_Irc(rows,:)    + F_Irc_s;
    F_Icr(:,rows)    = F_Icr(:,rows)    + F_Icr_s;

    Kcc_star = Kcc_star + Kcc_star_s;

    d_r(rows) = d_r(rows) + d_r_s;
    f_c_star = f_c_star + f_c_star_s;
end
Kcc_fac = decomposition(Kcc_star,'lu');

Kcc_inv_FIcr = Kcc_fac \ F_Icr;
Kcc_inv_fc   = Kcc_fac \ f_c_star;
toc
%% Solve Interface Problem
disp('Solving Interface Problem...')
tic
A_lambda = F_Irr + F_Irc * Kcc_inv_FIcr;
b_lambda = d_r - F_Irc * Kcc_inv_fc;

lambda = A_lambda \ b_lambda;
toc
%% Recover global corner dofs
disp('Recovering global corner dofs...')
tic
u_c_global = Kcc_fac \ (F_Icr * lambda + f_c_star);
toc


%% Recover Subdomain Solutions
disp('Recovering Subdomain Solutions...')
tic

SD_cell = cell(N_sd,1);

parfor s = 1:N_sd
    S = SD(s);

    Br = S.Br;
    Bc = S.Bc;

    lambda_s = lambda(S.lambda_idx);

    u_C = Bc * u_c_global;

    KRR_fac = decomposition(S.K_RR,'lu');

    u_R = KRR_fac \ (S.f_R - Br.' * lambda_s - S.K_RC * u_C);

    S.u_R = u_R;
    S.u_C = u_C;
    S.u = [u_R; u_C];

    S.u_i = u_R(1:S.nei);
    S.u_r = u_R(S.nei+(1:S.ner));
    S.u_c = u_C;

    SD_cell{s} = S;
end

SD = [SD_cell{:}];

toc
%% Assemble Global Solution from subdomain solutions

disp('Assembling Global Solution from FETI-DP subdomain solutions...')
tic

x_sum = complex(zeros(N_1,1));
x_cnt = zeros(N_1,1);

for s = 1:N_sd
    S = SD(s);

    % S.edge: reordered free local edge indices
    % S.edg_g(S.edge): corresponding global edge indices
    g_edges = S.edg_g(S.edge);

    x_sum(g_edges) = x_sum(g_edges) + S.u;
    x_cnt(g_edges) = x_cnt(g_edges) + 1;
end

x = complex(zeros(N_1,1));

free_hit = x_cnt > 0;
x(free_hit) = x_sum(free_hit) ./ x_cnt(free_hit);

% PEC dofs are Dirichlet zero
x(pec_edg) = 0;

toc


%% Visualization
disp('Reconstructing E field at element centers...')
tic

%현재 우리가 구한것들은 dof (edge) 위에서의 적분 값에 가까움.
%Visualization 하기 위해서는 어떤 포인트에서의 E field 값이 필요함.
% Element 내부(center)의 E-field 값은 우리가 구한 edge에서의 value와 그 edge의 whitney vector value의 선형조합으로 구성됨.
E_cent=complex(zeros(N_3,3));
lambda_cent=[1/4,1/4,1/4,1/4];
for s = 1:N_sd
    S=SD(s);
    for e_loc = 1:S.N_3
        e_g=S.ele(e_loc);

        grad=squeeze(S.grad_lambda(e_loc,:,:));
    
        W_cent=whitney1_value(lambda_cent,grad);
    
        E_vec=complex(zeros(1,3));
    
        for a = 1:6
            gid=ele_edg(e_g,a);
            sgn=ele_edg_sign(e_g,a);
    
            coeff_local=sgn*x(gid);
    
            E_vec=E_vec+coeff_local*W_cent(a,:);
        end
    
        E_cent(e_g,:)=E_vec;
    end
end
toc

%% Magnitude of E
E_mag = sqrt(abs(E_cent(:,1)).^2 + abs(E_cent(:,2)).^2 + abs(E_cent(:,3)).^2);

x_c = ele_cent_crdn(:,1);
y_c = ele_cent_crdn(:,2);
z_c = ele_cent_crdn(:,3);

%% XY plane visualization at z = z0

z0 = 0;

F_mag = scatteredInterpolant(x_c, y_c, z_c, E_mag, 'linear', 'none');

Nx = 200;
Ny = 200;

xq = linspace(min(x_c), max(x_c), Nx);
yq = linspace(min(y_c), max(y_c), Ny);

[Xq, Yq] = meshgrid(xq, yq);
Zq = z0 * ones(size(Xq));

Eq_mag = F_mag(Xq, Yq, Zq);

figure;
surf(Xq, Yq, Zq, Eq_mag, 'EdgeColor', 'none');
view(2);
axis equal tight;
colorbar;
xlabel('x');
ylabel('y');
title(['|E| on XY plane, z = ', num2str(z0)]);

%% YZ plane visualization at x = x0

x0 = 0;

F_mag = scatteredInterpolant(x_c, y_c, z_c, E_mag, 'linear', 'none');

Ny = 200;
Nz = 200;

yq = linspace(min(y_c), max(y_c), Ny);
zq = linspace(min(z_c), max(z_c), Nz);

[Yq, Zq] = meshgrid(yq, zq);
Xq = x0 * ones(size(Yq));

Eq_mag = F_mag(Xq, Yq, Zq);

figure;
surf(Xq, Yq, Zq, Eq_mag, 'EdgeColor', 'none');
view(90,0);
axis equal tight;
colorbar;
xlabel('y');
ylabel('z');
title(['|E| on YZ plane, x = ', num2str(x0)]);


%% ZX plane visualization at y = y0

y0 =0;

F_mag = scatteredInterpolant(x_c, y_c, z_c, E_mag, 'linear', 'none');

Nz = 200;
Nx = 200;

zq = linspace(min(z_c), max(z_c), Nz);
xq = linspace(min(x_c), max(x_c), Nx);

[Zq, Xq] = meshgrid(zq, xq);
Yq = y0 * ones(size(Xq));

Eq_mag = F_mag(Xq, Yq, Zq);

figure;
surf(Xq, Yq, Zq, Eq_mag, 'EdgeColor', 'none');
view(0,0);
axis equal tight;
colorbar;
xlabel('z');
ylabel('x');
title(['|E| on ZX plane, y = ', num2str(y0)]);




%%

figure;
quiver3( ...
    ele_cent_crdn(:,1), ele_cent_crdn(:,2), ele_cent_crdn(:,3), ...
    real(E_cent(:,1)), real(E_cent(:,2)), real(E_cent(:,3)), ...
    'AutoScale', 'on' ...
);

axis equal;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('Reconstructed E field at element centers - Real part');

%% Visualize |E| at element centers
E_mag = sqrt(abs(E_cent(:,1)).^2 + abs(E_cent(:,2)).^2 + abs(E_cent(:,3)).^2);

figure;
scatter3( ...
    ele_cent_crdn(:,1), ele_cent_crdn(:,2), ele_cent_crdn(:,3), ...
    20, E_mag, 'filled' ...
);

axis equal;
grid on;
colorbar;
xlabel('x');
ylabel('y');
zlabel('z');
title('|E| at element centers');

%% Visualize edge DOF along global edge directions
edg_cent = zeros(N_1,3);
edg_tvec = zeros(N_1,3);

for i = 1:N_1
    p1 = nod_crdn(edg_nod(i,1),:);
    p2 = nod_crdn(edg_nod(i,2),:);

    edg_cent(i,:) = 0.5 * (p1 + p2);

    t = p2 - p1;
    edg_tvec(i,:) = t / norm(t);
end

E_edge_vec = real(x) .* edg_tvec;

figure;
quiver3( ...
    edg_cent(:,1), edg_cent(:,2), edg_cent(:,3), ...
    E_edge_vec(:,1), E_edge_vec(:,2), E_edge_vec(:,3), ...
    'AutoScale', 'on' ...
);

axis equal;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('Edge DOF visualization along global edge tangents');

tend=toc(tstart)