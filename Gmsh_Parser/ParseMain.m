clc;
clear;

%% 2026.04.27 
%  Msh file version 4.1 에서는 Partition 된 Msh file을 Parsing 했을때 
%  Physical Tag가 제대로 인덱싱 되지 않음.
%  근본적인 이유로는 Parsing을 하면서 Entity index가 새롭게 바뀌게 되는데 이를 제대로 parsing 하지 못하는듯.

%% README
% geo 파일에서 Mesh의 기하설정을 할때 PML 영역은 PML_background 이런식으로 써줘야함 (대소문자 상관 x)
% 예를 들어 내부 박스가 박스 1이고 외부박스가 박스 2일때
% 박스 1과 2 사이에 pml 설정을 해주고 싶다고 하자. 그리고 박스 1 내부는 air라고 할때
% 박스 1과 2 사이의 영역을 'PML_Air'로 physical tag를 붙여줘야 background가 air로 설정됨.
%PML 설정은 현재 Natural Unit Solver를 기준으로 되어있음.

%PEC는 Physical Tag가 PEC인 노드,에지를 추출(대소문자 구분 x). 
%PEC는 물질 유전율로 처리하는것이 아니라 노드와 에지번호를 솔버로 넘겨서
%해당 노드, 인덱스의 dof를 제거하는 방식으로 진행

%%
load PML_set.txt % 1이면 on, 2면 off 
%% User Define
mshfile='untitled3.msh';  %사용자는 mshfile 이름만 바꿔주면 됨.
f=3e9;
omega=2*pi*f;
eps_0=8.854e-12;
L_0=0.1; % Natural Unit Solver에 사용될 Scale 
c0=299792458;
omega_tilde=(omega*L_0)/c0;
%% MSH 파일 버전 확인

fid=fopen(mshfile,'r');
if fid==-1
    error('파일을 열 수 없습니다: %s',mshfile);
end
Version=0;
while ~feof(fid)
        line = strtrim(fgetl(fid));
        if strcmp(line,'$MeshFormat')
            fmt=fgetl(fid);
            v_info = sscanf(fmt, '%f %d %d');
            Version=v_info(1);
            break;
           
        end
end
fclose(fid);
if Version==0
    error('MSH 파일 포맷을 인식할 수 없습니다.');
end
disp(['감지된 MSH Version : ', num2str(Version)]);

%% Parser 실행 
if Version <4.0
    Mesh=parseMSHV2(mshfile);
else
    Mesh=parseMSHV4(mshfile);
end

Phys_nametoID_keys=Mesh.Phys.NameToID.keys();
Phys_nametoID_values=Mesh.Phys.NameToID.values();
%% 모든 노드 좌표를 scaling 함.
Mesh.Nodes=Mesh.Nodes/L_0;

%% 요소,페이스,엣지,노드 connectivity 추출
Mesh=Build_Topology(Mesh);

%% PEC 정보 추출 (Physical Tag가 PEC인것 추출 (대소문자 구분 없이))
Mesh=Extract_PEC(Mesh);

%% Element_Center_coordinate 저장
Mesh=ele_cent_crdn(Mesh);

%% Physical Tag에 따라 물성치 부여 (eps, mu)
Mesh=Assign_Material(Mesh,omega,eps_0,omega_tilde,PML_set);
%% Mesh 정보 추출
writematrix(Mesh.Nodes(:,:),'nod_crdn.txt','Delimiter','space');
writematrix(Mesh.Elements{4}(:,2:5),'ele_nod.txt','Delimiter','space');
writematrix(Mesh.ElementPartitions{4}(:,:),'ele_part.txt');
writematrix(Mesh.ele_physTag,'ele_physTag.txt');
writematrix(f,'frequency_SI.txt');
writematrix(L_0,'scale.txt');
writematrix(omega_tilde,'omega_tilde.txt');
writematrix(Mesh.Topology.edg_nod, 'edg_nod.txt', 'Delimiter', 'space');
writematrix(Mesh.Topology.ele_edg, 'ele_edg.txt', 'Delimiter', 'space');
writematrix(Mesh.Topology.edg_length, 'edg_length.txt');

writematrix(Mesh.PEC.nod, 'pec_nod_list.txt');
writematrix(Mesh.PEC.edg, 'pec_edg_list.txt');
% Exporting Material Property
disp('Exporting Material Properties to TXT files...');
eps_r_9D=Mesh.eps_r;
mu_r_9D=Mesh.mu_r;
eps_r_1D=Mesh.eps_r(:,1);
mu_r_1D=Mesh.mu_r(:,1);

write_9D_complex_array('eps_r.txt',eps_r_9D);
write_9D_complex_array('mu_r.txt',mu_r_9D);
% write_1D_complex_array('eps_r_vec.txt',eps_r_1D);
% write_1D_complex_array('mu_r_vec.txt',mu_r_1D);
 
disp(' Complete Data Parsing and Export');

