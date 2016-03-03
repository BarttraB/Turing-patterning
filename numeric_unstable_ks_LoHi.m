%code for finding the ranges of unstable wavenumbers as a function of
%parameters..NOTE: code will crash when you lose the last fixed point, so
%need to rerun from after the first loop (line 92), then rerun the final plotting stuff at the
%end once the 2nd loop crashes. 

clear
reset(symengine)


%parameters

d1=0.004;
d2=0.00202;

a1=4.5;
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

numsamp=500;
numsamp2=2;

a1_vec=linspace(0.1,10,numsamp);

for j=1:numsamp

a1=a1_vec(j);        
        
syms L H A P positive

F(1)=(a1*(d1+P)/(c1+P))-g1*A*(L/(c3+L)) - k1*(p_max-P)*L + k2*P;
F(2)=(a2*(d2+P)/(c4+P))-g2*H;
F(3)=(a3*H)/(c5+H)-g3*A;
F(4)=k1*(p_max-P)*L - k2*P;

V=[L,H,A,P];

[FP.L,FP.H,FP.A,FP.P]=vpasolve([F==0], [V]);

j

%this keeps the lower (stable) fixed point when there are two
FPL=min(double(FP.L));
FPH=min(double(FP.H));
FPA=min(double(FP.A));
FPP=min(double(FP.P));


%if

J=jacobian(F,V);
Js=subs(J, [L,H,A,P], [FPL,FPH,FPA,FPP]);

syms ks positive

Dmat=[ks 0 0 0; 0 ks*D 0 0; 0 0 0 0; 0 0 0 0];


%CEL0=vpasolve(det(Js-Dmat), ks);

Poly=collect(det(Js-Dmat), ks);
PolyCoef=coeffs(Poly, ks);

Discrim=(PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1);
if(Discrim<0)
    ksH(j)=0;
    ksL(j)=0;
else
    ksH(j)=(-PolyCoef(2)+sqrt((PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1)))/(2*PolyCoef(3));
    ksL(j)=(-PolyCoef(2)-sqrt((PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1)))/(2*PolyCoef(3));
end


end

ns=size(nonzeros(subplus(double(ksL))));
start=j-ns(1);

figure
hold on
%eh = errorbar(a1_vec(start:j-1),meank(start:j-1), ksL(start:j-1), ksH(start:j-1))
%set(eh2,'Color','red')

%plot(a1_vec(start:j-1),ksH(start:j-1), 'r');
%plot(a1_vec(start:j-1),ksL(start:j-1), 'r')

X=[a1_vec(start:j-1),fliplr(a1_vec(start:j-1))];                %#create continuous x value array for plotting
Y=[ksL(start:j-1),fliplr(ksH(start:j-1))];              %#create y values for out and then back
fill(X,Y,'r');    


%%%%****case 2

a3=1.0;


a1_vec=linspace(1,3,numsamp);

for j=1:numsamp

a1=a1_vec(j); 
j
        
syms L H A P positive

F(1)=(a1*(d1+P)/(c1+P))-g1*A*(L/(c3+L)) - k1*(p_max-P)*L + k2*P;
F(2)=(a2*(d2+P)/(c4+P))-g2*H;
F(3)=(a3*H)/(c5+H)-g3*A;
F(4)=k1*(p_max-P)*L - k2*P;

V=[L,H,A,P];

[FP.L,FP.H,FP.A,FP.P]=vpasolve([F==0], [V]);

%this keeps the lower (stable) fixed point when there are two
FPL=min(double(FP.L));
FPH=min(double(FP.H));
FPA=min(double(FP.A));
FPP=min(double(FP.P));

J=jacobian(F,V);
Js=subs(J, [L,H,A,P], [FPL,FPH,FPA,FPP]);
syms ks positive
Dmat=[ks 0 0 0; 0 ks*D 0 0; 0 0 0 0; 0 0 0 0];
Poly=collect(det(Js-Dmat), ks);
PolyCoef=coeffs(Poly, ks);

Discrim=(PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1);

if(Discrim<0)
    ksH2(j)=0;
    ksL2(j)=0;
else
    ksH2(j)=(-PolyCoef(2)+sqrt((PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1)))/(2*PolyCoef(3));
    ksL2(j)=(-PolyCoef(2)-sqrt((PolyCoef(2)).^2-4.*PolyCoef(3).*PolyCoef(1)))/(2*PolyCoef(3));
end



end

ns=size(nonzeros(subplus(double(ksH2))));
start=j-ns(1);
%meank=(ksH+ksL)/2;
%rangek=ksH-ksL/2;

%eh2 = errorbar(a1_vec(start:j-1),meank(start:j-1), ksL(start:j-1), ksH(start:j-1))
%plot(a1_vec(start:j-1),ksH2(start:j-1), 'b');
%plot(a1_vec(start:j-1),ksL2(start:j-1), 'b')

X2=[a1_vec(start:j-1),fliplr(a1_vec(start:j-1))];                %#create continuous x value array for plotting
Y2=[ksL2(start:j-1),fliplr(ksH2(start:j-1))];              %#create y values for out and then back
fill(X2,Y2,'b');    



xlim([0,9])
xlabel('\alpha_1 [min.^{-1}]','FontSize',16)
ylabel('unstable k^2','FontSize',16)
legend('\alpha_3=3 min.^{-1}','\alpha_3=1 min.^{-1}', 'location','northwest')
legend('boxoff')


set(gca,'FontSize',16)
%set(eh2,'Color','blue','linestyle','none')
set(gcf, 'PaperPosition', [0.5 0.5 6.5 6.5]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
set(gcf, 'PaperSize', [7 7]); %Keep the same paper size

saveas(gcf, 'wavenums_model_E4_a1_a3is1and3_1', 'pdf')