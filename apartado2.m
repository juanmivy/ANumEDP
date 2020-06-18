function [x,uint] = apartado2
%
% Resuelve u_t + a u_x = 0  en [ax,bx] con condiciones de 
% contorno periodicas empleando el metodo de Bean Warming 
% con m nodos espaciales interiores.
%
% Representa la solucion en varios instantes de tiempo.
% Devuelve mallado y  aproximacion en la ultima iteracion.
% Modificado del código propuesto en clase para la simulación de 
% la ecuación de advección con el método de Lax-Wendroff


printf(" Comenzamos introduciendo los datos del problema, por motivos técnicos la
condicion inicial debe introducirse desde el código\n");
ax= input('Introduzca el inicio del intervalo espacial: ');
bx= input('Introduzca el final del intervalo espacial: ');
a= input('Introduzca el coeficiente de transporte: ');
tfinal= input('Introduzca el tiempo final: ');
eta= @(x) x; %Intoducir aqui la condicion inicual deseada
h= input('Introduzca el paso en espacio: ');
k= input('Introduzca el paso en tiempo: ');



% Datos de las discretizaciones
m=round((bx-ax)/h)-1;      % nodos en espacio
nu = a*k/h;                % numero de Courant
x = linspace(ax,bx,m+2)';  % Fijamos el mallado en el intervalo espacial
nsteps = round(tfinal / k);    % numero de pasos temp
nplot = 20;       % representamos cada nplot pasos

% Comprobamos si el ultimo paso llega a tfinal
if abs(k*nsteps - tfinal) > 1e-5
  disp(' ')
  disp(sprintf('OJO *** k no divide a tfinal, k = %9.5e',k))
  disp(' ')
end

% En el caso con condiciones periodicas se tienen:
% m+1 incognitas u(2:m+2)  y   u(1) = u(m+2)
I = 2:(m+2);   % indices de las incognitas

% Definicion de condiciones iniciales:
tn = 0;
u0 = eta(x);
u = u0;

% Condiciones de contorno periodicas:
u(1) = u(m+2); 
% OJO: nodo fantasma para cond periodicas necesitamos dos nodos "fantasma"  
u(m+3) = u(2);   
u(m+4)=u(3);
% representamos condicion inicial:
clf
plot(x,u0)
axis([ax bx -5 10]) %axis([0 1 -.2 1.2])
title('Condicion inicial a tiempo 0')
input('Pulse <return> para continuar ');

% Bucle temporal:
for n = 1:nsteps
  tnp = tn + k;   % = t_{n+1}

  % Implementamos el método de Beam-Warming:
  if a<0
    % Condiciones de contorno periodicas:
  u(1) = u(m+2); 
% OJO: nodo fantasma para cond periodicas necesitamos dos nodos "fantasma"  
  u(m+3) = u(2);   
  u(m+4)=u(3);
  u(I) = u(I) - 0.5*nu*(-3*u(I)+4*u(I+1)-u(I+2)) + 0.5*nu^2 * (u(I) - 2*u(I+1) + u(I+2));
endif
if a>0
  % Condiciones de contorno periodicas:
  u(m+2) = u(1); 
% OJO: nodo fantasma para cond periodicas necesitamos dos nodos "fantasma"  
  u(m+3) = u(2);   
  u(m+4)=u(3);
  u(I) = u(I) - 0.5*nu*(3*u(I)-4*u(I+1)+u(I+2)) + 0.5*nu^2 * (u(I) - 2*u(I+1) + u(I+2));
endif

  % Condiciones periodicas:
  u(1) = u(m+2);   % redundante para dato inicial bien preparado
  % OJO: los nodos fantasma para cond periodicas
  u(m+3) = u(2);   
  u(m+4)=u(3);
  % representamos resultados cada nplot pasos:
  if mod(n,nplot)==0 || n==nsteps
    uint = u(1:m+2);  % representamos nodos en el intervalo
    plot(x,uint)
    axis([ax bx -5 10]) % axis([0 1 -.2 1.2])
    title(sprintf('t = %9.5e  tras %4i pasos con %5i nodos',...
                       tnp,n,m+1))
    if n<nsteps, pause(.5); end;
  end

  tn = tnp;   % para el siguiente paso temporal
end % del for
uint = u(1:m+2);
endfunction 