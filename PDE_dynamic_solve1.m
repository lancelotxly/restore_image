function [fd,b,a,time]=PDE_dynamic_solve1(p,iter,eth,Nub,Rho,K,dt1,dt2)
[NY,NX]=size(p);
b=ones(NY,NX);
f=p;
time=0;
for i=1:iter
c=b.*b;
f_x_e = f(:,[2:NX,NX])-f;
   cE = 0.5*(c(:,[2:NX,NX])+c);
f_x_w  = f-f(:,[1 1:NX-1]);
   cW = 0.5*(c(:,[1 1:NX-1])+c);
f_y_s = f([2:NY,NY],:)-f;
   cS = 0.5*(c([2:NY,NY],:)+c);
f_y_n = f-f([1 1:NY-1],:);
   cN = 0.5*(c([1,1:NY-1],:)+c);
TermN=cN.*f_y_n;
TermS=cS.*f_y_s;
TermE=cE.*f_x_e;
TermW=cW.*f_x_w;
div1 = TermE - TermW + TermS - TermN;
f1=f;
f=f+dt1*(p-f+(Nub^2)*div1);         %%%%%%%%%%%最速下降迭代%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fN=[f(1,:); f(1:NY-1,:)]-f;      
    fS=[f(2:NY,:); f(NY,:)]-f;
    fE=[f(:,2:NX) f(:,NX)]-f;
    fW=[f(:,1) f(:,1:NX-1)]-f;
    
    psi_N=(abs(fN).^2)/K+1;
    psi_S=(abs(fS).^2)/K+1;
    psi_E=(abs(fE).^2)/K+1;
    psi_W=(abs(fW).^2)/K+1;
    part1=b.*(psi_N+psi_S+psi_E+psi_W);
    
    bN=[b(1,:); b(1:NY-1,:)]-b;      
    bS=[b(2:NY,:); b(NY,:)]-b;
    bE=[b(:,2:NX) b(:,NX)]-b;
    bW=[b(:,1) b(:,1:NX-1)]-b;
%     
    phi_N=1./((1+abs(bN).^2).^2);             %%%%%%%%%%%%%G-M%%%%%%
    phi_S=1./((1+abs(bS).^2).^2);
    phi_E=1./((1+abs(bE).^2).^2);
    phi_W=1./((1+abs(bW).^2).^2);
%     
%      phi_N=1;                               %%%%%%%%%%%%T-V%%%%%%%
%      phi_S=1;
%      phi_E=1;
%      phi_W=1;
% 
%     phi_N=tanh(abs(bN))./(2*abs(bN+1e-6));       %%%%%%%%%%%Green%%%%%%%%
%     phi_S=tanh(abs(bS))./(2*abs(bS+1e-6));
%     phi_E=tanh(abs(bE))./(2*abs(bE+1e-6));
%     phi_W=tanh(abs(bW))./(2*abs(bW+1e-6));
    
%     phi_N=exp(-abs(bN).^2);           %%%%%%%%%%%%%%P-M%%%%%%%%%%%%%%
%     phi_S=exp(-abs(bS).^2);
%     phi_E=exp(-abs(bE).^2);
%     phi_W=exp(-abs(bW).^2);
    g=(phi_N+phi_S+phi_E+phi_W)/4;
    

b_x_e = b(:,[2:NX,NX])-b;
   cE = 0.5*(phi_E+g);
b_x_w  = b-b(:,[1 1:NX-1]);
   cW = 0.5*(phi_W+g);
b_y_s = b([2:NY,NY],:)-b;
   cS = 0.5*(phi_S+g);
b_y_n = b-b([1 1:NY-1],:);
   cN = 0.5*(phi_N+g);
TermN=cN.*b_y_n;
TermS=cS.*b_y_s;
TermE=cE.*b_x_e;
TermW=cW.*b_x_w;
div2 = TermE - TermW + TermS - TermN;

b1=b;
b=b+dt2*(1-part1+(Rho^2)/(K^2)*div2);

a=sum(sum(abs(f1-f)));
time=time+1;
if(a<eth) break;
end
end
fd=f;
b=b;
time=time;
a=a;