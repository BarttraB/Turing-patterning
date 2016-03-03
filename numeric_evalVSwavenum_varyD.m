
clear
reset(symengine)


%parameters

d1=0.004;
d2=0.00202;

a1=4.0;
a2=1.6667;
a3=2.5;


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


%loop parameter
numpts=100;
ks_vec=linspace(0,1,numpts);

syms L H A P positive

F4(1)=(a1*(d1+P)/(c1+P))-g1*A*(L/(c3+L)) - k1*(p_max-P)*L + k2*P;
F4(2)=(a2*(d2+P)/(c4+P))-g2*H;
F4(3)=(a3*H)/(c5+H)-g3*A;
F4(4)=k1*(p_max-P)*L - k2*P

V=[L,H,A,P];

%solve for fixed point

FP4=vpasolve([F4==0], [V]);
FP4.L

J=jacobian(F4,V);

%plugging in for the three fied points (FP)
%Js=struct('FP1',subs(J, [L,H,A], [FP.L(1),FP.H(1),FP.A(1)]), 'FP2', subs(J, [L,H,A], [FP.L(2),FP.H(2),FP.A(2)]), 'FP3', subs(J, [L,H,A], [FP.L(3),FP.H(3),FP.A(3)]))
Js=subs(J, [L,H,A,P], [FP4.L(1),FP4.H(1),FP4.A(1),FP4.P(1)]);

syms Lam

for i=1:numpts
    ks=ks_vec(i);
    Dmat1=[ks 0 0 0; 0 ks*D 0 0; 0 0 0 0; 0 0 0 0];
    JsLamD1=vpasolve([det(Js-Lam*eye(4)-Dmat1)==0], [Lam]);
%    RLam1(i,1)=double(real(JsLamD1.FP1));   
    RLam4(i,1:4)=double(real(JsLamD1));
%    RLam3(i,1:3)=double(real(JsLamD1.FP3));
end

figure
hold on


%plot(ks_vec,RLam4(:,4),'b*-');
%plot(ks_vec,RLam4(:,3),'b*-');



%a3 case 2

D=50;


F5(1)=(a1*(d1+P)/(c1+P))-g1*A*(L/(c3+L)) - k1*(p_max-P)*L + k2*P;
F5(2)=(a2*(d2+P)/(c4+P))-g2*H;
F5(3)=(a3*H)/(c5+H)-g3*A;
F5(4)=k1*(p_max-P)*L - k2*P

V=[L,H,A,P];

%solve for fixed point

FP5=vpasolve([F5==0], [V]);
FP5.L

J=jacobian(F5,V);

Js=subs(J, [L,H,A,P], [FP5.L(1),FP5.H(1),FP5.A(1),FP5.P(1)]);

for i=1:numpts
    ks=ks_vec(i);
    Dmat1=[ks 0 0 0; 0 ks*D 0 0; 0 0 0 0; 0 0 0 0];
    JsLamD1=vpasolve([det(Js-Lam*eye(4)-Dmat1)==0], [Lam]);
%    RLam1(i,1)=double(real(JsLamD1.FP1));   
    RLam5(i,1:4)=double(real(JsLamD1));
%    RLam3(i,1:3)=double(real(JsLamD1.FP3));
end


%a3 case 3

D=25;


F6(1)=(a1*(d1+P)/(c1+P))-g1*A*(L/(c3+L)) - k1*(p_max-P)*L + k2*P;
F6(2)=(a2*(d2+P)/(c4+P))-g2*H;
F6(3)=(a3*H)/(c5+H)-g3*A;
F6(4)=k1*(p_max-P)*L - k2*P

V=[L,H,A,P];

%solve for fixed point

FP6=vpasolve([F6==0], [V]);
FP6.L

J=jacobian(F6,V);

Js2=subs(J, [L,H,A,P], [FP6.L(1),FP6.H(1),FP6.A(1),FP6.P(1)]);
%Js=subs(J, [L,H,A,P], [FP6.L(2),FP6.H(2),FP6.A(2),FP6.P(2)]);

for i=1:numpts
    ks=ks_vec(i);
    Dmat1=[ks 0 0 0; 0 ks*D 0 0; 0 0 0 0; 0 0 0 0];
    JsLamD1=vpasolve([det(Js-Lam*eye(4)-Dmat1)==0], [Lam]);
    %JsLamD2=vpasolve([det(Js2-Lam*eye(4)-Dmat1)==0], [Lam]);

%    RLam1(i,1)=double(real(JsLamD1.FP1));   
    RLam6(i,1:4)=double(JsLamD1);
    %RLam6_2(i,1:4)=double(real(JsLamD2));

%    RLam3(i,1:3)=double(real(JsLamD1.FP3));
end


plot(ks_vec,RLam4(:,4),'bs-');
plot(ks_vec,RLam5(:,4),'b*-');
plot(ks_vec,real(RLam6(:,4)),'bd-');
%plot(ks_vec,RLam6_2(:,4),'rd-');


%plot(ks_vec,RLam4(:,3),'rs-');
%plot(ks_vec,RLam5(:,3),'r*-');
%plot(ks_vec,RLam6(:,3),'rd-');

%plot(ks_vec,RLam5(:,2),'gs-');
%plot(ks_vec,RLam5(:,1),'ks-');
%ylim([-1 1])

 xlabel('k^2','FontSize',14); ylabel('Re(\lambda)','FontSize',14);
hh= legend('\alpha_{3}=1.1','\alpha_{3}=0.9','\alpha_{3}=0.7, L_*=0.04128','\alpha_{3}=0.7, L_*=1.175', 'Location','NorthEast')
legend('boxoff')
%set(hh,'Visible', 'Off')

set(gca,'FontSize',14) 



%set(gcf, 'PaperPosition', [0.5 0.5 6.5 6.5]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
%set(gcf, 'PaperSize', [7 7]); %Keep the same paper size
%saveas(gcf, 'evals_models', 'pdf')

set(gcf, 'PaperSize', [7 4]); %Keep the same paper size
set(gcf, 'PaperPosition', [0.5 0.5 6.5 3.5]);
%saveas(gcf, 'evals_model_E4_c3', 'pdf')


%FP=vpasolve([f1==0, f2==0, f3==0], [L, H, A])