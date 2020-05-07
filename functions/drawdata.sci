// Copyright (C) 2017 - AGH - Wojciech T. Witkowski
//
// Date of creation: 2017-05-15
//
function drawdata()
clc;
xset("window",3);
xset("thickness",1);
f=get("current_figure");
f.figure_size=[1200 900];
p=gca();
p.isoview="on";

eksp = findobj("tag", "Aquifer:"); eksp = (eksp.string);
pkt = findobj("tag", "Profile:"); pkt = (pkt.string);

u=file('open', eksp ,'unknown');
E=read(u,-1,2);
E1=size(E);
file('close',u);

uu=file('open', pkt ,'unknown');
W=read(uu,-1,4);
W1=size(W);
file('close',uu);
Wx=W(:,2);
Wy=W(:,3);

ilosc_pol=E(1,1);
licznik_wierszy=3;
licznik=1;

for i=1:ilosc_pol
    ilosc_wierzcholkow=E(licznik_wierszy,2);
    for j=1:ilosc_wierzcholkow
        Xe(licznik)=E(licznik_wierszy+j,1);
        Ye(licznik)=E(licznik_wierszy+j,2);
        licznik=licznik+1;
    end
    licznik_wierszy=licznik_wierszy+j+2;
end;

plot2d(0,0,-1,"012"," ",[min(Xe)-100,min(Ye)-100,max(Xe)+100,max(Ye)+100]);
Xe=0;
Ye=0;

f=gcf();
f.color_map=jetcolormap(ilosc_pol);

licznik_wierszy=3;
licznik=1;
for i=1:ilosc_pol
    ilosc_wierzcholkow=E(licznik_wierszy,2);
    for j=1:ilosc_wierzcholkow
        Xe(licznik)=E(licznik_wierszy+j,1);
        Ye(licznik)=E(licznik_wierszy+j,2);
        licznik=licznik+1;
    end
    
    xfpoly(Xe, Ye );
   
    p=get("hdl"); 
    p.background=i;
    licznik=1;
    licznik_wierszy=licznik_wierszy+j+2;
    Xe=0;
    Ye=0;
    wysw=int(i/E(1)*100);
    disp("Progress:"+string(wysw)+" %");
end;
disp(i, "Number of elements in aquifer:");

xpoly(Wx, Wy, "lines");
p=get("hdl"); 
p.foreground=1;
p.thickness=1;
p.mark_style=9;
xpoly(Wx(1), Wy(1), "lines");
p=get("hdl"); 
p.thickness=3;
p.mark_style=3;
endfunction
