close all;
clc;
clear all;

[Aur,Alr,Aut,Alt] = GeomtoCST();
x0 = [5.6,0.9,14,25,25,0,0,Aur',Alr',Aut',Alt']';
lb = [0.5,0.1,5,0.5,0.5,-5,-5,0.5*Aur',0.5*Alr',0.5*Aut',0.5*Alt']';
ub = [20,20,30,60,60,5,5,2*Aur',2*Alr',2*Aut',2*Alt']';

x0n = (x0-lb)./(ub-lb);
lbn = 0.*lb;
ubn = ub./ub;
xinit = x0n.*(ub-lb)+lb;
kinkloc = 0.5;
[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(Aur,Alr,Aut,Alt,kinkloc);