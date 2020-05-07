// Copyright (C) 2017 - AGH - Wojciech T. Witkowski
//
// Date of creation: 2017-05-15
//
function calculation(liczba)
clc;
messagebox("Computations in progress. Please wait...");
tic();
// retrive data
eksp = findobj("tag", "Aquifer:"); eksp = (eksp.string);
pkt = findobj("tag", "Profile:"); pkt = (pkt.string);
a = findobj("tag", "$C_{m} [{MPa}^{-1}]:$"); a = evstr(a.string);
tgb = findobj("tag", "$tg\beta [-]:$"); tgb = evstr(tgb.string);
il_iter = findobj("tag", "No. Iterations:"); il_iter = evstr(il_iter.string);

if liczba == 1 then il_iter=1 
end

u=file('open', eksp ,'unknown');
E=read(u,-1,2);
E1=size(E);
file('close',u);

v=file('open', pkt ,'unknown');
L=read(v,-1,4);
L1=size(L);
file('close',v);

my_plot_axes = gca();
delete(gca());
newaxes();
set(gca(),"auto_clear","off");
xlabel("Distance in profile [m]");
ylabel("Surface subsidence [mm]");

my_plot_axes = gca();
my_plot_axes.title.font_size = 2;
my_plot_axes.axes_bounds = [1/4,0,4/5,1];
my_plot_axes.title.text = "Theoretical and Empirical Surface Subsidence";

Hm=0;
gm=0;
Xtr=0;
Ytr=0;
value=1;
Mk=0.002;

Wo=ones(L1(1,1),2);
Wo(:,1)=L(:,1);
Wo(:,2)=-L(:,4);

x2=0;
x2(1,1)=0;
for i=1:L1(1,1)-1
    odl=sqrt((L(i,2)-L(i+1,2))^2+(L(i,3)-L(i+1,3))^2);
    x2(1,i+1)=odl+x2(1,i);
end;

plot2d(x2, Wo(:,2), 5, leg="Empirical Surface Subsidence");
plot(x2,Wo(:,2),'ro--');
mac_a=0;
mac_tgb=0;
V=0;
    deff('zp=W(x,y)','zp=exp((-%pi/H^2)*(tgb^2)*((x-s)^2+(y-t)^2))');
    deff('z2=W2(x,y)','z2=-2*%pi*tgb*((x-s)^2+(y-t)^2)/H^2*exp((-%pi/H^2)*(tgb^2)*((x-s)^2+(y-t)^2))');

kk=1;
while value>Mk
    Wt=zeros(L1(1),2);
    Wt(:,1)=L(:,1);
    A=zeros(L1(1),2);
    licz_w=0;
    for j=1:E(1,1)
            Wtt=0;
            Aa=0;
            Xtr=0;
            Ytr=0;
        il_w=E(licz_w+3,2);
        licz_tr=0;
        X0=E(licz_w+4,1);
        Y0=E(licz_w+4,2);
        pi=E(licz_w+2,1);
        H=E(licz_w+2,2);
        M=E(licz_w+3,1)/1000;
            for k=1:il_w-2
                Xtr(1,licz_tr+k)=X0;
                Xtr(2,licz_tr+k)=E(licz_w+3+k+1,1);
                Xtr(3,licz_tr+k)=E(licz_w+3+k+2,1);
        
                Ytr(1,licz_tr+k)=Y0;
                Ytr(2,licz_tr+k)=E(licz_w+3+k+1,2);
                Ytr(3,licz_tr+k)=E(licz_w+3+k+2,2);
            end;
            
        for i=1:L1(1,1)
            s=L(i,2);
            t=L(i,3);
            g=(pi)*M;
            
            if (sqrt((s-mean(Xtr))^2+(t-mean(Ytr))^2)<(3*H/tgb)) then 
                I=-((a*g)/(H^2))*(tgb^2)*1000*int2d(Xtr,Ytr,W, [1.d-10, 1, 50, 4000, 0]);
                clc;
                disp("No. of iterations: "+string(kk));disp("Progress: "+string(int(j/E(1,1)*100))+" [%]");
                Wtt(i,2)=I;
                
                Ix=2*Wtt(i,2)/(tgb*1000);
                Iy=int2d(Xtr,Ytr,W2, [1.d-10, 1, 50, 4000, 0]);
                clc;
                disp("No. of iterations: "+string(kk));disp("Progress: "+string(int(j/E(1,1)*100))+" [%]");
                Iy=-((a*g)/(H^2))*(tgb^2)*Iy;
                Aa(i,2)=Ix+Iy;
            else 
                Wtt(i,2)=0;
                Aa(i,2)=0;clc;disp("No. of iterations: "+string(kk));disp("Progress: "+string(int(j/E(1,1)*100))+" [%]");
            end;
        end;
        
    for i=1:L1(1,1)
        Wt(i,2)=Wt(i,2)+Wtt(i,2);
        A(i,2)=A(i,2)+Aa(i,2);
    end;
    licz_tr=licz_tr+k;
    licz_w=licz_w+il_w+2;
    end;

    A(:,1)=Wt(:,2)/a/1000;
    csvWrite(Wt, 'functions\Model_results.txt', ';');
    l=Wo-Wt;
    LL=l(:,2)/1000; 
    X=(inv(A'*A))*A'*LL;
    mac_a(kk)=a;
    mac_tgb(kk)=tgb;
    a=a+X(1);
    tgb=tgb+X(2);
    plot2d(x2, Wt(:,2), kk, leg="Theoretical Surface Subsidence");
    V=Wt(:,2)- Wo(:,2);
    VV=0;
    for i=1:L1(1,1)
        V(i,1)=V(i,1)^2;
        VV=VV+V(i,1);
    end;
    V(kk,1)=(sqrt(s/(L1(1,1)-1))/(-min(Wo)))*100;
    value=abs(X(2));
    if kk==il_iter then value=0.0000001;
    end
    kk=kk+1;
end;
kk=kk-1;
set(gca(),"auto_clear","on");

d=LL-A*X;
sigma_kw=((d'*d)/((L1(1,1)-2)));
Cov_X=sigma_kw*inv((A'*A));
od_st_a=sqrt(Cov_X(1,1));
od_st_tgb=sqrt(Cov_X(2,2));
r=Cov_X(1,2)/(od_st_a*od_st_tgb);
sig_kw_max=(Cov_X(1,1)+Cov_X(2,2))/2+sqrt(((Cov_X(1,1)-Cov_X(2,2))/2)^2+Cov_X(1,2));
sig_kw_min=(Cov_X(1,1)+Cov_X(2,2))/2-sqrt(((Cov_X(1,1)-Cov_X(2,2))/2)^2+Cov_X(1,2));
kat_skr_b=0.5*atan((2*Cov_X(1,2))/(Cov_X(1,1)-Cov_X(2,2)));

nic=0;
wyn=mopen('functions\value_in_iterations.txt','w+');
mfprintf(wyn, "%.2f %.2f\n", kk,nic);
for i=1:kk
    mfprintf(wyn, "%.2f %.2f\n", mac_a(i), mac_tgb(i));
end
mclose(wyn);

time=int(toc()*100/60)/100;
clc;
disp("Progress: "+string(j/E(1,1)*100)+" [%]");
disp("Time: "+string(time)+" [min]");

if liczba == 1 then messagebox("Calculations completed in "+string(time)+" [min]");
end
if liczba == 2 then 
wyn=mopen('Results.txt','w+');
mfprintf(wyn, "Value of compaction coefficient (Cm [MPa-1]) with standard deviation:\t %.10f   +/- %.10f \n", a, od_st_a);
mfprintf(wyn, "Value of influence coefficient (tgB [-]) with standard deviation:\t %.2f   +/- %.2f \n \n", tgb, od_st_tgb);
mfprintf(wyn, "Value of correlation (R [-]): \t %.2f \n", r);
mclose(wyn);

messagebox(["Calculations completed in "+string(time)+" [min]", "No. Iterations: "+string(kk)+"/"+string(il_iter), "Value of parameters in file: Results.txt"], "Results");
end

endfunction
