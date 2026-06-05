function [lambda_q,weight_q]=Tet_Quad_4point()

alpha=0.5854101966249685; % 알려진 quadrautre variable
beta=0.1381966011250105;

lambda_q=[
    alpha beta beta beta;
    beta alpha beta beta;
    beta beta alpha beta;
    beta beta beta alpha];

weight_q=ones(4,1)/4;
end