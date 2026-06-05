function W=whitney1_value(lambda,grad)

loc_edg_nod=[1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
W=zeros(6,3);

for a = 1:6
    i=loc_edg_nod(a,1);
    j=loc_edg_nod(a,2);

    W(a,:)=lambda(i)*grad(j,:)-lambda(j)*grad(i,:);
end
end
