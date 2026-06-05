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

%% Mesh Information getter
disp('Building Mesh Information (Connectivity)...')
tic
[N_0,N_1,N_2,N_3,DT,...
    ele_fac,ele_edg,ele_nod,ele_edg_sign,...
    fac_edg,fac_nod,...
    edg_nod,...
    edg_length,...
    nod_crdn]=FEM_3D_Mesh_Getter_gmsh(nod_crdn,ele_nod);
toc
%% Getting Barycentric Information
disp('Building Barycentric information...')
tic
[ele_vol, grad_lambda, coef_lambda,curl_W]=FEM_3D_Barycentric_Gradient(nod_crdn,ele_nod,N_3);
toc

%% Quadrature Information
disp('Building Information of 4-Quadrature Points in a Element ...')
tic
[lambda_q,weight_q]=Tet_Quad_4point();
toc
%% Getting Whitney Forms
disp('Building Whitney Forms...')
tic
%각 element의 4개의 quadrature 지점에서의
%W_ij(x_q)=lambda_i(x_q)*grad_lambda_j-lambda_j(x_q)*grad_lambda_i 계산
%whitney 1 form은 edge가 만든 element 내부의 vector field라고 이해하고 넘어가자.
[whitney_1]=Whitney_Forms_Quad(N_3,grad_lambda,lambda_q);  
toc

%% Local Matrices
disp('Building Local Matrices...')
tic

I_all=zeros(36,N_3);
J_all=zeros(36,N_3);
V_all=complex(zeros(36,N_3));

parfor e = 1:N_3
    [A_loc,Meps_loc, Mmu_loc, Kcurl_loc, Ccurl_loc]=FEM_3D_Local_Matrices(e,ele_vol,weight_q,whitney_1,curl_W,eps_r,mu_r,omega_tilde);

    I_e=zeros(36,1);
    J_e=zeros(36,1);
    V_e=complex(zeros(36,1));

    cnt=0;
    for a =1:6
        ia=ele_edg(e,a); %현재 element e의 a번째 edge의 글로벌 index
        sa=ele_edg_sign(e,a);

        for b=1:6
            ib=ele_edg(e,b); % 현재 element e의 b번째 edge의 글로벌 index
            sb=ele_edg_sign(e,b);
            
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
toc

%% Global FEM Matrix Assemble
disp('Building Global FEM Matrices...');
tic

I_vec=I_all(:);
J_vec=J_all(:);
V_vec=V_all(:);
A=sparse(I_all,J_all,V_all,N_1,N_1);

toc

%% Source
b=complex(zeros(N_1,1));
src_edg=1;
b(src_edg)=1;

%% Apply PEC condition
pec_edg=unique(pec_edg_list(:));
pec_nod=unique(pec_nod_list(:));

all_edg=(1:N_1).';
free_edg=setdiff(all_edg,pec_edg); %pec가 아닌 dof
A_free=A(free_edg,free_edg);
b_free=b(free_edg);

%% Solve
disp('Solving Linear Equation...')
tic
x_free=A_free\b_free;

x=complex(zeros(N_1,1));
x(free_edg)= x_free;
x(pec_edg)=0;
toc

%% Visualization
disp('Reconstructing E field at element centers...')
tic

%현재 우리가 구한것들은 dof (edge) 위에서의 적분 값에 가까움.
%Visualization 하기 위해서는 어떤 포인트에서의 E field 값이 필요함.
% Element 내부(center)의 E-field 값은 우리가 구한 edge에서의 value와 그 edge의 whitney vector value의 선형조합으로 구성됨.
E_cent=complex(zeros(N_3,3));
lambda_cent=[1/4,1/4,1/4,1/4];
for e = 1:N_3

    grad=squeeze(grad_lambda(e,:,:));

    W_cent=whitney1_value(lambda_cent,grad);

    E_vec=complex(zeros(1,3));

    for a = 1:6
        gid=ele_edg(e,a);
        sgn=ele_edg_sign(e,a);

        coeff_local=sgn*x(gid);

        E_vec=E_vec+coeff_local*W_cent(a,:);
    end

    E_cent(e,:)=E_vec;
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

% Visualize |E| at element centers
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

% Visualize edge DOF along global edge directions
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

%%
tend=toc(tstart)