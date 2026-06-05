function [A_loc,Meps_loc, Mmu_loc, Kcurl_loc, Ccurl_loc]=FEM_3D_Local_Matrices(e,ele_vol,weight_q,whitney_1,curl_W,eps_r,mu_r,omega)

Nq=length(weight_q);

Meps_loc=zeros(6,6);
Mmu_loc=zeros(6,6);
Kcurl_loc=zeros(6,6);
Ccurl_loc=zeros(6,6);

vol=ele_vol(e);

eps_tensor=reshape(eps_r(e,:),[3,3]).'; %현재 
mu_tensor=reshape(mu_r(e,:),[3,3]).';

mu_inv_tensor=inv(mu_tensor);

curlW=squeeze(curl_W(e,:,:)); %e번 element의 curl whitney

for q = 1:Nq
    W=squeeze(whitney_1(e,q,:,:)); %6 x 3 %e번 element의 q번 point에서의 whitney 1 form 
    
    for a= 1:6
        Wa=W(a,:); % e번 element의 q번 point에서의 a번 edge에 대한 whitney 1 form vector 값.
        
        %하나의 edge마다 element 내부의 전체 edge와의 계산을 해야함.
        for b = 1:6
            Wb=W(b,:); 
            

            Meps_loc(a,b)=Meps_loc(a,b)+ vol*weight_q(q)*(Wa*eps_tensor*Wb.');
            % E-only formulation에서는 필요없음.
            % curlWb=curlW(b,:);
            % Mmu_loc(a,b)=Mmu_loc(a,b)+vol*weight_q(q)*(Wa*mu_tensor*Wb.');
            % 
            % Ccurl_loc(a,b)=Ccurl_loc(a,b)+vol*weight_q(q)*dot(Wa,curlWb);
        end
    end
end
for a = 1:6
    curlWa=curlW(a,:);

    for b= 1:6
        curlWb=curlW(b,:);

        Kcurl_loc(a,b)=vol*(curlWa*mu_inv_tensor*curlWb.');
    end
end
A_loc=Kcurl_loc-omega^2*Meps_loc;
end

