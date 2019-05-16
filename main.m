clear all; clc; close all;
Original = imread('BoartsColor.bmp');
img1 =rgb2gray(Original);
figure(1) ; imshow(img1,[]);
figure(2);stairs(img1(100,:),'b');hold on;
p = imnoise(Original,'gaussian',0,0.025); 
p=rgb2gray(p);
stairs(p(100,:),'y');
legend('Orginal','noise');hold off;
ylabel('intensity');
figure(3) ; imshow(p,[]);
p=double(p);
snr1=SNR(p,img1);
psnr1=PSNR(p,img1);


iter =1000;    %迭代次数
eth=1.5;
dt1=1e-4;
dt2=1e-6;
r=2;

if(r==2)
    Nub=7.8;
    Rho=7.8;       %最佳50 
    K=5;          %最佳3
    [f,b,a,d]=PDE_dynamic_solve1(p,iter,eth,Nub,Rho,K,dt1,dt2);
    figure(4);imshow(f,[]);
    figure(5); imshow(b,[]);
end
if(r==1)
    Nub=40;
    Rho=50;       %最佳50 
    K=32;          %最佳3
    [f,b,a,d]=PDE_dynamic_solve2(p,iter,eth,Nub,Rho,K,dt1,dt2);
    f=uint8(f);
    figure(4); imshow(f,[]);
    figure(5); imshow(b,[]);
end
figure(6); stairs(f(100,:),'r');hold on;
stairs(img1(100,:),'b');
legend('restore','Orginal');hold off;
ylabel('intensity');
snr2=SNR(f,img1);
psnr2=PSNR(f,img1);