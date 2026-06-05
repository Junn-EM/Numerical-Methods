function [whitney_1]=Whitney_Forms_Quad(N_3,grad_lambda,lambda_q)

Nq=size(lambda_q,1); % 4 

loc_edg_nod=[1,2;1,3;1,4;2,3;2,4;3,4];

whitney_1=zeros(N_3,Nq,6,3); % 각 element는 4개의 quadrature point에서 6 x 3 짜리 whitney 1-form 값이 필요함.

for e=1:N_3
    grad=squeeze(grad_lambda(e,:,:));

    for q=1:Nq
        lambda=lambda_q(q,:);

        for a=1:6
            i=loc_edg_nod(a,1);
            j=loc_edg_nod(a,2);

            whitney_1(e,q,a,:)=lambda(i)*grad(j,:)-lambda(j)*grad(i,:);
        end
    end
end
end
