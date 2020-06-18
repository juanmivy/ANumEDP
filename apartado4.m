function [x,uint] = apartado4
%
% Resuelve u_t + a u_x = 0  en [ax,bx] con condiciones de 
% contorno periodicas empleando el metodo de Bean Warming 
% con m nodos espaciales interiores.
%
% Representa la solucion en varios instantes de tiempo.
% Devuelve mallado y  aproximacion en la ultima iteracion.
% Modificado del código propuesto en clase para la simulación de 
% la ecuación de advección con el método de Lax-Wendroff

m=input('Introduzca el numero de nodos en espacio deseado: ')

% Comenzamos introduciendo los datos del problema
global a
a = -1;           % parametro de la ecuacion
ax = 0;           % inicio del intervalo esspacial
bx = 10;          % fin del intervalo espacial
tfinal = 10;                % Tiempo maximo
%eta= @(x) (0) .* (x < 3 ) +((4-(x-5)).^2 .* (x >= 3 & x <=7 ) + (0) .* (x >7 & x <=10 ); % Cond inicial
eta=@(x) (0).*(x <= 3)+(4-(x -5).^2).*(x >= 3 & x <= 7 ) +(0).*(x>7) 
% Datos de las discretizaciones
h = (bx-ax)/(m+1);         % h paso espacial
k = 0.2*h;                  % k paso temporal
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
axis([0 10 -1 5]) %axis([0 1 -.2 1.2])
title('Condicion inicial a tiempo 0')
 %%Procedemos a calcular la integral%%
   vector_integral=h.*u0;
   valor_integral=sum(vector_integral);
   disp("El valor de la integral a tiempo cero es="),disp(valor_integral);
   clear vector_integral,valor_integral;
input('Pulse <return> para continuar ');


% Bucle temporal:
for n = 1:nsteps
  tnp = tn + k;   % = t_{n+1}

  % Implementamos el método de Beam-Warming:
  u(I) = u(I) - 0.5*nu*(-3*u(I)+4*u(I+1)-u(I+2)) + 0.5*nu^2 * (u(I) - 2*u(I+1) + u(I+2));

  % Condiciones periodicas:
  u(1) = u(m+2);   % redundante para dato inicial bien preparado
  % OJO: los nodos fantasma para cond periodicas
  u(m+3) = u(2);   
  u(m+4)=u(3);
  % representamos resultados cada nplot pasos:
  if mod(n,nplot)==0 || n==nsteps
    uint = u(1:m+2);  % representamos nodos en el intervalo
    plot(x,uint)
    axis([0 10 -1 5]) % axis([0 1 -.2 1.2])
    title(sprintf('t = %9.5e  tras %4i pasos con %5i nodos',...
                       tnp,n,m+1))
    if n<nsteps, pause(.5); end;
    %%Procedemos a calcular la integral%%
   vector_integral=h.*uint;
   valor_integral=sum(vector_integral);
   disp("El valor de la integral a tiempo(primera salida) es (segunda salida)"),disp(tn),disp(valor_integral);
   clear vector_integral,valor_integral;
  end
 
  tn = tnp;   % para el siguiente paso temporal
end % del for
uint = u(1:m+2);
endfunction 
