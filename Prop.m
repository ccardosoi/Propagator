clear all
close all
% Parámetros de la ventana numérica
h = 1280;
v = 1024;
w0 = 1e-3;
xmax = 7.6e-3;              % Tamaño horizontal SLM Holoeye
ymax = 4.3e-3;              % tamañ vertical SLM

% h = 7.6e-3;              % spatial extent of the grid
% v = 4.3e-3;
% Nh = 1280;               % number of samples
% Nv = 1024;
% deltah = h/Nh;           % sample spacing
% deltav = v/Nv;

% xs = xmax*((1:1280)-(1280/2));
% ys = xmax*((1:1024)-(1024/2));
xs = xmax*(2/h)*(-h/2:h/2-1);
ys = ymax*(2/v)*(-v/2:v/2-1);
% xs = linspace(-Nh/2,Nh/2-1,Nh)*deltah;
% ys = linspace(-Nv/2,Nv/2-1,Nv)*deltav;
kxmax = 1/(1280*(xs(2)-xs(1)));
kymax = 1/(1024*(ys(2)-ys(1)));
% kx = kmax*((1:1280)-(1280/2))/(xs(2)-xs(1));
% ky = kmax*((1024/2)-(1:1024))/(ys(2)-ys(1));
kx = kxmax*(-h/2:h/2-1);
ky = kymax*(-v/2:v/2-1);
% kx = linspace(-Nh/2,Nh/2-1,Nh)/(Nh*deltah);
% ky = linspace(-Nv/2,Nv/2-1,Nv)/(Nv*deltav);

dx = xs(2)-xs(1); dy = ys(2)-ys(1);
dkx = kx(2)-kx(1); dky = ky(2)-ky(1);

% % Parámetros para SLM
% Oz = 0.3*pi/180;
% Oxy = 45*pi/180;
k = 2*pi/633e-9; 
% kx = k.*sin(Oz).*cos(Oxy);
% ky = k.*sin(Oz).*sin(Oxy);

% Generamos la malla
[Xs, Ys] = meshgrid(xs,ys);
[Kx, Ky] = meshgrid(kx,ky);

% Definicion del modo LG
U = @(x,y,w0,m) ((sqrt(x.^2 + y.^2)/w0).^abs(m)) .* exp(-(x.^2 + y.^2)/(w0^2)) .* exp(1i*m*atan2(y,x));

% Apertura circular
CA = zeros(1024,1280);
R = 0.1e-3;
A = (Xs.^2 + Ys.^2 <= R^2);             % circular aperture of radius R
CA(A) = 1;

% Círculos
C = 2;
S = 0.5e-3;
M = zeros(1024,1280);
for i=1:C
    A = ((Xs+S*cos((i-1)*pi/(C/2))).^2 + (Ys+S*sin((i-1)*pi/(C/2))).^2 <= R^2);
    M(A) = 1;
end

% Generación del modo LG
Up = U(Xs,Ys,w0,1);
c = sum(sum(Up.*conj(Up)));
LGp = Up./sqrt(c);

% Parámetros iniciales
    % Incrementos en z 
        dz = 6e-1;
    % Distancia de propagación
        dp = 100.2;
        n = dp/dz;
        
        dp = 300;
        n = 100;
        dz = dp/n;
    % Fase Global
        G = exp(1i*k*dz);
    % Fase de propagación
        kz = -(1i/2)*(((Kx.^2 + Ky.^2))./k);
        P = 1; %exp(1i*kz*dz);
    
    % Transformada del haz inciial
        Uin = LGp;
        U0 = fftshift(fft2(fftshift(Uin)));
    
    % Inicialización Haz de salida
%         Uout = CA;
%         Uout1 = fftshift(ifft2(fftshift(U0)));
%         Uout2 = fftshift(ifft2(fftshift(Uout)));

% Propagación
for j=1:n
    Uout(:,:,j) = fftshift(ifft2(fftshift(U0.*P)));
    P = P.*exp(kz*dz);
%     G = G.*G;
end
    
% Gráfica plano transversal
figure (1)
subplot(1,2,1)
imagesc(abs(Uin))
title('Gráfica plano transversal campo entrada')
xlabel('x')
ylabel('y')
subplot(1,2,2)
imagesc(abs(Uout(:,:,20)))
title('Gráfica plano transversal campo salida')
xlabel('x')
ylabel('y')

figure (2)
subplot(1,2,1)
imagesc(angle(Uout(:,:,1)))
title('Fase plano transversal campo entrada')
xlabel('x')
ylabel('y')
axis square
subplot(1,2,2)
imagesc(angle(Uout(:,:,n)))
title('Fase plano transversal campo salida')
xlabel('x')
ylabel('y')
axis square

% Uout1 = Uout(:,:,1);
% Uouta1 = Uout(:,640,1);
% Uout2 = Uout(:,:,2);
% Uouta2 = Uout(:,640,2);
% Uout3 = Uout(:,:,3);
% Uouta3 = Uout(:,640,3);
% Uout4 = Uout(:,:,4);
% Uouta4 = Uout(:,640,4);
% Uout5 = Uout(:,:,5);
% Uouta5 = Uout(:,640,5);
% Uout6 = Uout(:,:,6);
% Uouta6 = Uout(:,640,6);
% Uout7 = Uout(:,:,7);
% Uouta7 = Uout(:,640,7);
% Uout8 = Uout(:,:,8);
% Uouta8 = Uout(:,640,8);
% Uouta21 = Uout(:,640,21);
% Uouta20 = Uout(:,640,20);
% Uouta = Uout(512,:,50);



zp = zeros(1280,n);
for jj = 1:n
    zp(:,jj) = Uout(512,:,jj);
%     for ii = 1:1280
%         if zp(ii,jj)>1
%             zp(ii,jj) = 1;
%         elseif zp(ii,jj)<-1
%             zp(ii,jj) = -1;
%         end
%     end
end

% n = ones(1024,n);

% Gráfica plano de propagación
figure (3)
imagesc(abs(zp))
title('Gráfica plano de propagación')
xlabel('Prop z')
ylabel('x')


