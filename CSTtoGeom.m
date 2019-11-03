%Function D_airfoil2 to transform input Bernstein parameters (shape + class function method) to
%complete Airfoil coordinates

%output [upper surface y coord, lower surface y coord, Class function pos y
%coord, up surf thickness distribution, lw surf thickness distb, camber
%distb] = input (up surf Bernstein parameters, lw surf BS params,
%X-ordinates)
function[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(Aur,Alr,Aut,Alt,Xr,Xt,kinkloc)

xr = Xr(:,1);
xt = Xt(:,1);
N1 = 0.5;   %Class function N1
N2 = 1;     %Class function N2

zeta_u = 0.000;     %upper surface TE gap
zeta_l = -0.000;     %lower surface TE gap


nur = length(Aur)-1;
nlr = length(Alr)-1;
nut = length(Aut)-1;
nlt = length(Alt)-1;

for i = 1:length(xr)
    
    %calculate Class Function for x(i):
    C(i) = (xr(i)^N1)*(1-xr(i))^N2;
    
    %calculate Shape Functions for upper and lower surface at x(i)
    Su(i) = 0;  %Shape function initially zero
    for j = 0:nur
        Krnu = factorial(nur)/(factorial(j)*factorial(nur-j));
        Su(i) = Su(i) + Aur(j+1)*Krnu*(1-xr(i))^(nur-j)*xr(i)^(j);
    end
    Sl(i) = 0;  %Shape function initially zero
    for k = 0:nlr        
        Krnl = factorial(nlr)/(factorial(k)*factorial(nlr-k));
        Sl(i) = Sl(i) + Alr(k+1)*Krnl*(1-xr(i))^(nur-k)*xr(i)^(k);
    end
    
    %calculate upper and lower surface ordinates at x(i)
    Yur(i) = C(i)*Su(i) + xr(i)*zeta_u;
    Ylr(i) = C(i)*Sl(i) + xr(i)*zeta_l;
    
    Thur(i) = C(i)*(Su(i)-Sl(i))/2;    %calculate thickness distribution !TE thickness ignored!
    Thlr(i) = C(i)*(Sl(i)-Su(i))/2;    %calculate thickness distribution !TE thickness ignored!
    Cmr(i) = C(i)*(Su(i)+Sl(i))/2;    %calculate camber distribution !TE thickness ignored!
end


for i = 1:length(xt)
    
    %calculate Class Function for x(i):
    C(i) = (xt(i)^N1)*(1-xt(i))^N2;
    
    %calculate Shape Functions for upper and lower surface at x(i)
    Su(i) = 0;  %Shape function initially zero
    for j = 0:nut
        Krnu = factorial(nut)/(factorial(j)*factorial(nut-j));
        Su(i) = Su(i) + Aut(j+1)*Krnu*(1-xt(i))^(nut-j)*xt(i)^(j);
    end
    Sl(i) = 0;  %Shape function initially zero
    for k = 0:nlt        
        Krnl = factorial(nlt)/(factorial(k)*factorial(nlt-k));
        Sl(i) = Sl(i) + Alt(k+1)*Krnl*(1-xt(i))^(nut-k)*xt(i)^(k);
    end
    
    %calculate upper and lower surface ordinates at x(i)
    Yut(i) = C(i)*Su(i) + xt(i)*zeta_u;
    Ylt(i) = C(i)*Sl(i) + xt(i)*zeta_l;
    
    Thut(i) = C(i)*(Su(i)-Sl(i))/2;    %calculate thickness distribution !TE thickness ignored!
    Thlt(i) = C(i)*(Sl(i)-Su(i))/2;    %calculate thickness distribution !TE thickness ignored!
    Cmt(i) = C(i)*(Su(i)+Sl(i))/2;    %calculate camber distribution !TE thickness ignored!
end
Ytuk = (Ytur-Ytut)*(1-kinkloc) + Ytut;
Ytlk = (Ytlr-Ytlt)*(1-kinkloc) + Ytlt;
Ytu85 = (Ytur-Ytut)*(1-0.85) + Ytut;
Ytl85 = (Ytlr-Ytlt)*(1-0.85) + Ytlt;

Yustt = Yut';
Ylstt = Ylt';
Yusrt = Yur';
Ylsrt = Ylr';
Ytuks = Ytuk';
Ytlks = Ytlk';
Ytu85s = Ytu85';
Ytl85s = Ytl85';
%assemble airfoil coord matrix

xk = (xr-xt)*(1-kinkloc) + xt;
x85 = (xr-xt)*(1-0.85) + xt;
Xtur = [xr  Yusrt];
Xtlr = [xr  Ylsrt];
Xtut = [xt  Yustt];
Xtlt = [xt  Ylstt];
Xtuk = [xk  Ytuks];
Xtlk = [xk  Ytlks];
Xtu85 = [x85  Ytu85s];
Xtl85 = [x85  Ytl85s];

