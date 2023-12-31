% ��Legendre����ʽ��ͼ
% 
clear
% ����ʽ������0��N
N = 6;
syms P x
% ��ͼ������ɫ����
c={'b','r','g','k','m','y','c'};
for n = 0:N
    % Legendre����ʽ���޵���˱��ʽ��Ҳ���õ��ƹ�ϵ���㣩
    P = 1/(2^n*prod(1:n))*diff((x^2-1)^n,x,n);
    h = fplot(P,[-1 1],'LineWidth',1);
    set(h,'Color',c{n+1})
    hold on
end
legend('n=0','n=1','n=2','n=3','n=4','n=5','n=6','Location','NorthEastOutside')
title('Legendre polynomials');
xlabel('x')