function [fd,b,a,d]=PDE_dynamic_solve2(p,iter,eth,Nub,Rho,K,dt1,dt2)
[NY,NX]=size(p);
b=ones(NY,NX);
f=p;
d=0;
% laplacianOpreator=[0,1,0;1,-4,1;0,1,0];
for i=1:iter
    fN=[f(1,:); f(1:NY-1,:)]-f;      
    fS=[f(2:NY,:); f(NY,:)]-f;
    fE=[f(:,2:NX) f(:,NX)]-f;
    fW=[f(:,1) f(:,1:NX-1)]-f;
    
    psi_N=1./abs(fN+1e-6);
    psi_S=1./abs(fS+1e-6);
    psi_E=1./abs(fE+1e-6);
    psi_W=1./abs(fW+1e-6);
    
    TermN=psi_N.*fN;
    TermS=psi_S.*fS;
    TermE=psi_E.*fE;
    TermW=psi_W.*fW;
    div1=TermN+TermS+TermE+TermW;
    
    bN=[b(1,:); b(1:NY-1,:)]-b;      
    bS=[b(2:NY,:); b(NY,:)]-b;
    bE=[b(:,2:NX) b(:,NX)]-b;
    bW=[b(:,1) b(:,1:NX-1)]-b;
    part1=2*b.*(bN.*TermN+bS.*TermS+bE.*TermE+bW.*TermW);

    
    f1=f;
    f=f+dt1*((p-f)+(Nub^2)/2*((b.*b.*div1+part1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fN=[f(1,:); f(1:NY-1,:)]-f;      
    fS=[f(2:NY,:); f(NY,:)]-f;
    fE=[f(:,2:NX) f(:,NX)]-f;
    fW=[f(:,1) f(:,1:NX-1)]-f;
    
    phi_N=abs(fN)/(K^2)+1;
    phi_S=abs(fS)/(K^2)+1;
    phi_E=abs(fE)/(K^2)+1;
    phi_W=abs(fW)/(K^2)+1;
    part1=b.*(phi_N+phi_S+phi_E+phi_W);
    
    bN=[b(1,:); b(1:NY-1,:)]-b;      
    bS=[b(2:NY,:); b(NY,:)]-b;
    bE=[b(:,2:NX) b(:,NX)]-b;
    bW=[b(:,1) b(:,1:NX-1)]-b;
    
%     bhi_N=1;                        %%%%%%%%%%%%Tihonkov%%%%%%%%%%%%%%
%     bhi_S=1;
%     bhi_E=1;
%     bhi_W=1;
    
%       bhi_N=1./(2*abs(bN)+eth);      %%%%%%%%%%%%%%%TV%%%%%%%%%%%%%%
%       bhi_S=1./(2*abs(bS)+eth);
%       bhi_E=1./(2*abs(bE)+eth);
%       bhi_W=1./(2*abs(bW)+eth);

%     bhi_N=exp(-abs(bN).^2);           %%%%%%%%%%%%%%P-M%%%%%%%%%%%%%%
%     bhi_S=exp(-abs(bS).^2);
%     bhi_E=exp(-abs(bE).^2);
%     bhi_W=exp(-abs(bW).^2);

%     bhi_N=1./(2*sqrt(1+abs(bN).^2));     %%%%%%%%%%%%%%%Hyper%%%%%%%%
%     bhi_E=1./(2*sqrt(1+abs(bE).^2));
%     bhi_S=1./(2*sqrt(1+abs(bS).^2));
%     bhi_W=1./(2*sqrt(1+abs(bW).^2));
    
%     bhi_N=2./(1+abs(bN).^2);          %%%%%%%%%%%%%%H-L%%%%%%%%%%%%%%
%     bhi_S=2./(1+abs(bS).^2);
%     bhi_E=2./(1+abs(bE).^2);
%     bhi_W=2./(1+abs(bW).^2);
    
    bhi_N=1./((1+abs(bN).^2).^2);             %%%%%%%%%%%%%%%%%%G-M%%%%%%%%
    bhi_S=1./((1+abs(bS).^2).^2);
    bhi_E=1./((1+abs(bE).^2).^2);
    bhi_W=1./((1+abs(bW).^2).^2);
    
   
%     bhi_N=tanh(abs(fN))./(2*abs(fN+1e-6));       %%%%%%%%%%%Green%%%%%%%%
%     bhi_S=tanh(abs(fS))./(2*abs(fS+1e-6));
%     bhi_E=tanh(abs(fE))./(2*abs(fE+1e-6));
%     bhi_W=tanh(abs(fW))./(2*abs(fW+1e-6));
g=(bhi_N+bhi_S+bhi_E+bhi_W)/4;   
b_x_e = b(:,[2:NX,NX])-b;
   cE = 0.5*(bhi_E+g);
b_x_w  = b-b(:,[1 1:NX-1]);
   cW = 0.5*(bhi_W+g);
b_y_s = b([2:NY,NY],:)-b;
   cS = 0.5*(bhi_S+g);
b_y_n = b-b([1 1:NY-1],:);
   cN = 0.5*(bhi_N+g);
TermN=cN.*b_y_n;
TermS=cS.*b_y_s;
TermE=cE.*b_x_e;
TermW=cW.*b_x_w;
div2 = TermE - TermW + TermS - TermN;
  
b1=b;
b=b+dt2*(1-part1+(Rho^2)/K*div2);

a=sum(sum(abs(f1-f)));
if(a<eth) break;
end
d=d+1;
end
fd=f;
b=b;
a=a;
d=d;