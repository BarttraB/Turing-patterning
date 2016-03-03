%code for finding the ranges of unstable wavenumbers as a function of
%parameters..NOTE: code will crash when you lose the last fixed point, so
%need to rerun from after the first loop (line 92), then rerun the final plotting stuff at the
%end once the 2nd loop crashes. 
tic

clear
% close all

numsamp=60;
numsamp2=20;

a1_vec=linspace(2.5,5.3,numsamp2);

Pmp11=zeros(1,numsamp2);

for jj=1:numsamp2

reset(symengine)

Pm_vec=linspace(0.2,15,numsamp);

Pmp1=zeros(1,numsamp); %rezero every outer loop

%parameters

d1=0.004;
d2=0.00202;

a1=a1_vec(jj)
a2=1.6667;
%a3=3.5;

a3=3.0;


g1=3.0;
g2=2.0;
g3=1.1;

c1=0.675;
%c2=1.0;
c3=0.01;
c4=0.675;
c5=0.675;

k1=100;
k2=248.15;
p_max=1;

D=100;

jj

for j=1:numsamp

p_max=Pm_vec(j);
j
        
syms L H A P positive

F(1)=(a1*(d1+P)/(c1+P))-g1*A*(L/(c3+L)) - k1*(p_max-P)*L + k2*P;
F(2)=(a2*(d2+P)/(c4+P))-g2*H;
F(3)=(a3*H)/(c5+H)-g3*A;
F(4)=k1*(p_max-P)*L - k2*P;

V=[L,H,A,P];

[FP.L,FP.H,FP.A,FP.P]=vpasolve([F==0], [V]);

 
% 
% %this keeps the lower (stable) fixed point when there are two
% FPL(j)=min(double(FP.L));
% FPH(j)=min(double(FP.H));
% FPA(j)=min(double(FP.A));
% FPP(j)=min(double(FP.P));

    FPL(j)=FP.L(1);
    FPH(j)=FP.H(1);
    FPA(j)=FP.A(1);
    FPP(j)=FP.P(1);
if (numel(FP.L)>1) 
    FPL(j)=FP.L(2);
    FPH(j)=FP.H(2);
    FPA(j)=FP.A(2);
    FPP(j)=FP.P(2);
end
    
J=jacobian(F,V);
% FP.L
Js=subs(J, [L,H,A,P], [FPL(j),FPH(j),FPA(j),FPP(j)]);
Js1=double(Js);

syms ks positive

Dmat=[ks 0 0 0; 0 ks*D 0 0; 0 0 0 0; 0 0 0 0];


%CEL0=vpasolve(det(Js-Dmat), ks);

Poly=collect(det(Js-Dmat), ks);
PolyCoef=coeffs(Poly, ks);
Discrim(j)=(PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1);
%h1(j)=double(det(Js1-Dmat)); 

%Polys=collect(det(J-Dmat), ks);
%PolyCoefs=coeffs(Polys, ks)
%Discrims=(PolyCoefs(2)).^2-4.*PolyCoefs(3).*PolyCoefs(1);
%Discrimss=subs(Discrims,[L,H,A,P], [FPL(j),FPH(j),FPA(j),FPP(j)])

ksH(j)=double((-PolyCoef(2)+sqrt((PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1)))/(2*PolyCoef(3)));
ksL(j)=double((-PolyCoef(2)-sqrt((PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1)))/(2*PolyCoef(3)));


        if(j>1 && Discrim(j)>0 && Discrim(j-1)<0 && ksH(j)>0 && ksL(j)>0)
            Pmp1(j)=Pm_vec(j)-Discrim(j)*(Pm_vec(j)-Pm_vec(j-1))/(Discrim(j)-Discrim(j-1));
%             kk=j
        end

        
syms Lam        
JsLam=vpasolve([det(Js1-Lam*eye(4))==0], [Lam]);
RLam5(j,1:4)=double(JsLam);     
        
end

Discrim_n=double(Discrim);

%restart here after error message

ns=size(nonzeros(double(Discrim)));
start=j-ns(1);




figure
hold on
%plot(a1_vec(start:j-1),Discrim(start:j-1))
subplot(2,1,1);plot(Pm_vec,Discrim);%plot(a3_vec,double(h1));
ylabel('Discriminant') 
subplot(2,1,2);plot(Pm_vec, ksH);hold; plot(Pm_vec,ksL,'r');
ylabel('k^2') 
xlabel('\alpha1') 

% figure
% subplot(3,1,1);plot(a3_vec,FPL)
% ylabel('L_*') 
% %subplot(2,1,2);plot3( real(RLam5), imag(RLam5), a3_vec)
% subplot(3,1,2);plot(a3_vec,real(RLam5))
% ylabel('Re(\lambda)') 
% 
% subplot(3,1,3);plot(a3_vec,imag(RLam5))
% ylabel('Im(\lambda)') 
% xlabel('\alpha3') 



dPmp1=double(Pmp1);
Pmp11(jj)=nonzeros(dPmp1)

end 

toc

% a1_vec1=a1_vec(1:11);
% Pmp111=Pmp11(1:11);
figure
plot(a1_vec,Pmp11,'.-')


