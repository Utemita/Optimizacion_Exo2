% Analisis de posicion del mecanismo propuesto para el rehabilitador de
% mano

clear all;
close all;
clc;

% Por conveniencia se establece que las variables THETA escritas en
% mayusculas correspondan a angulos medidos en grados, y las variables
% "theta" en minuscula correspondan a angulos medidos en radianes (rad)

% Parametros de diseno propuestos por el disenador
%Dedo índice
Bancada1 = 18;
Bancada2 = 20;
Link1 = 35;
Link2 = 49;
Link3 = 25;
Link4 = 20;
Link5 = 25;
Link6 = 55;
Link7 = 35;
Link8 = 52;
Link9 = 25;
Link10 = 35;
c2 = 46.01;             %Distancia entre la art. IFP y el punto P3
TETHA2inicial = 0;
TETHA1inicial = 109;
THETA14B = 90;
hsp = 17; 
dsp = 18;
fp = 49;
fm = 26;
fd = 24;
THETAauxfm = 51.39; 
THETAauxfd = 38.78; 

% %Dedo Medio
% Link1 = 35;
% Link2 = 49;
% Link3 = 25;
% Link4 = 20;
% Link5 = 25;
% Link6 = 55;
% Link7 = 35;
% Link8 = 52;
% Link9 = 25;
% Link10 = 35;
% c2 = 46.81;             %Distancia entre la art. IFP y el punto P3
% TETHA2inicial = 0;
% TETHA1inicial = 109;
% THETA14B = 90;
% hsp = 17; 
% dsp = 18;
% fp = 47.91;
% fm = 27.52;
% fd = 24;
% THETAauxfm = 54.34795726; 
% THETAauxfd = 30; 
% 
% %Dedo Pulgar
% Link1 = 35;
% Link2 = 49;
% Link3 = 25;
% Link4 = 20;
% Link5 = 25;
% Link6 = 55;
% Link7 = 35;
% Link8 = 52;
% Link9 = 25;
% Link10 = 35;
% c2 = 46.81;             %Distancia entre la art. IFP y el punto P3
% TETHA2inicial = 0;
% TETHA1inicial = 109;
% THETA14B = 90;
% hsp = 17; 
% dsp = 18;
% fp = 47.91;
% fm = 27.52;
% fd = 24;
% THETAauxfm = 54.34795726; 
% THETAauxfd = 30;

%%%%%%%%%%%%%%%%%%%%
omega2 = 10;
reng = 2;
% Asignacion de los parametros de diseno a las variables utilizadas en las
% ecuaciones

% Primer mecanismo de 5 barras
r4 = Link1;
r5 = Link2;
r2 = Link3;
r1 = Link4;
r3 = Bancada1/2;

% Primer mecanismo de 4 barras
a = Link4;
b = Link5;
c = sqrt(hsp^2+dsp^2); % Hipotenusa de los soportes de la falange proximal
d = Bancada2;

% Segundo mecanismo de 5 barras
r1m2 = fp - 2*dsp;  %Distancia entres soportes S1 y S2
r2m2 = Link7;
r3m2 = Link5/2;     %Mitad del Link 5 del mecanismo completo
r4m2 = Link3; 
r5m2 = Link6;

% Segundo mecanismo de 4 barras
a2 = Link7;
b2 = Link8;
% c2 se ingresa como parametro de diseno directo debido a su posicion
d2 = sqrt( hsp^2 + dsp^2); %Hipotenusa de los soportes de la falange prox.
 
% Enseguida se presenta una breve descripcion de las variables utilizadas
% en el programa y los valores recomendados para pruebas

% Bancada1 = 18;
% Bancada2 = 20;
% Link1 = 35;
% Link2 = 49;
% Link3 = 25;
% Link4 = 20;
% Link5 = 25;
% Link6 = 55;
% Link7 = 35;
% Link8 = 52;
% Link9 = 25;
% Link10 = 35;
% c2 = 46.81;           %Distancia entre la art. IFP y el punto P3
% TETHA2inicial = 0;    %Angulo de la manivela de ENTRADA en la posicion 
%                       estacionaria en grados;
% TETHA1inicial = 109;  %Angulo de la manivela del mecanismo de 4 barras en 
%                       la posicion estacionaria en grados;
% THETA14B = 90;        %Angulo de la bancada del mecanismo de 4 barras
%                       inicial en grados con respecto a la horizontal
% hsp = 17;             %Altura de los soportes de la falanges
% dsp = 18;             %distancia del soporte falange prox. a la art.
% fp = 47.91;           %Largo de la falange proximal
% fm = 27.52;           %Largo de la falange medial
% fd = 24;              %Largo de la falange distal
% THETAauxfm = 54.34795726; 
% THETAauxfd = 30; 
% omega2=10;            %Velocidad de la manivela en RPM
% reng = 2;             %Relacion de engranaje propuesta

% THETA1 corresponde al angulo de entrada del mecanismo de 4 barras y
% angulo de la segunda junta actuada del primer mecanismo de 5 barras

% Se pasan los angulos a radianes y se calcula el tiempo total del recorrido
theta14B= THETA14B*pi/180; %Angulo de la Bancada en radianes
OMEGA2=omega2*2*pi/60; %Velocidad de la manivela en rad/seg.
%   Para la vuelta completa
tiempototal=(2*pi)/OMEGA2; %tiempo total del recorrido

%Division de 360 grados y del tiempo en una vuelta completa de la manivela
jnodos=360;
dtheta=360/(jnodos);
dtiempo=tiempototal/jnodos;

%for j=1:jnodos+1 %Para recorrido de 360 grados
for j=1:133     % SOlo recorre 132 grados por conveniencia (modificable)
    
tiempo(j)=(j-1)*dtiempo; %Iniciamos en tiempo cero
THETA2(j)=(j-1)*dtheta+TETHA2inicial; % Angulo de la manivela del M5Barras 
theta2(j)=THETA2(j)*pi/180; %angulo de la manivela del M5Barras en radianes

THETA1(j)=THETA2(j)/reng + TETHA1inicial; %De la rel. de engranaje maniv.4B
theta1(j)=THETA1(j)*pi/180;  %Angulo de manivela 4B en radianes




%CALCULOS PARA EL MECANISMO INICIAL DE 5 BARRAS

%Relaciones obtenidas de las ecuaciones
e = ( r1*sin(theta1(j))- r4*sin(theta2(j)) )/( r4*cos(theta2(j)) - r1*cos(theta1(j)) + 2*r3);
temp1 = ( 2*(r1*r3*cos(theta1(j)) + r3*r4*cos(theta2(j))) - r1^2 + r2^2 + r4^2 - r5^2 );
f = temp1/( 2*(r4*cos(theta2(j)) - r1*cos(theta1(j)) + 2*r3) );
%Coeficientes de la ecuacion cuadratica
daux = e^2 + 1;
g = 2*( e*f - e*r1*cos(theta1(j)) + e*r3 - r1*sin(theta1(j)) );
h = f^2 - 2*f*( r1*cos(theta1(j)) - r3 ) - 2*r1*r3*cos(theta1(j)) + r1^2 + r3^2 - r2^2;

pyP(j) = ( - g + sqrt( g.^2 - 4*daux*h) )/ (2*daux); 
pxP(j) = e*pyP(j) + f;

% pyN(j) = ( - g - sqrt( g.^2 - 4*daux*h) )/ (2*daux);
% pxN(j) = e*pyN(j) + f;



% CALCULOS PARA EL MECANISMO INICIAL DE 4 BARRAS
k1=d*cos(theta14B)+a*cos(theta1(j));
k2=d*sin(theta14B)+a*sin(theta1(j));
k3=k1^2+k2^2+c^2-b^2;
k4=-c^2+b^2+k1^2+k2^2;
A1=-k3-2*k1*c;
B1=4*k2*c;
C1=2*k1*c-k3;
A2=k4-2*k1*b;
B2=4*k2*b;
C2=2*k1*b+k4;

%Calculo de posiciones angulares
tantheta31=(-B2+sqrt(B2^2-4*A2*C2))/(2*A2);
tantheta32=(-B2-sqrt(B2^2-4*A2*C2))/(2*A2);
tantheta41=(-B1+sqrt(B1^2-4*A1*C1))/(2*A1);
tantheta42=(-B1-sqrt(B1^2-4*A1*C1))/(2*A1);
theta3c(j)=2*atan(tantheta31);
theta3a(j)=2*atan(tantheta32);
theta4c(j)=2*atan(tantheta41);
theta4a(j)=2*atan(tantheta42);
THETA3c(j)=2*atan(tantheta31)*180/pi; % theta3 cruzado en grados
THETA3a(j)=2*atan(tantheta32)*180/pi; % theta3 abierto en grados
if THETA3a(j)<0
    THETA3a(j) = 360 + THETA3a(j);
end
THETA4c(j)=2*atan(tantheta41)*180/pi; %theta4 cruzado en grados
THETA4a(j)=2*atan(tantheta42)*180/pi; %theta4 abierto en grados
if THETA4a(j)<0
    THETA4a(j) = 360 + THETA4a(j);
end

%Calculo de la posicion de la falange proximal
THETAfp(j) = THETA4a(j)+ atand(hsp/dsp);
thetafp(j)= deg2rad(THETAfp(j)); %Angulo de la falange proximal en radianes
pxIFP(j) = fp*cos(thetafp(j)) - d*cos(theta14B) - r3;
pyIFP(j) = fp*sin(thetafp(j)) - d*sin(theta14B);

%Calculo del sujetador de la falange proximal
pxS1(j) = c*cos(theta4a(j)) - d*cos(theta14B) - r3;
pyS1(j) = c*sin(theta4a(j)) - d*sin(theta14B);
THETAps2(j) = THETAfp(j) - atand(hsp/(fp-dsp));
thetaps2(j)= deg2rad(THETAps2(j)); %Angulo de la falange proximal en radianes
rs2 = sqrt(hsp^2 + (fp - dsp)^2);
pxS2(j) = rs2*cos(thetaps2(j)) - d*cos(theta14B) - r3;
pyS2(j) = rs2*sin(thetaps2(j)) - d*sin(theta14B);

%Calculo del punto Pm4
pxM4(j) = a*cos(theta1(j)) - r3;
pyM4(j) = a*sin(theta1(j));



% CALCULOS PARA EL SEGUNDO MECANISMO DE 5 BARRAS CON RESPECTO AL 
% SISTEMA SECUNDARIO
THETAroll(j) = atan2d(pyM4(j) - pyS1(j), pxM4(j) - pxS1(j));
thetaroll(j) = deg2rad(THETAroll(j)); %Angulo de cambio de los ejes en radianes

THETA1m2so(j) = atan2d(pyS2(j) - pyS1(j), pxS2(j) - pxS1(j));
if THETA1m2so(j)<0
    THETA1m2so(j) = 360 + THETA1m2so(j);
end
theta1m2so(j) = deg2rad(THETA1m2so(j)); %Angulo de entrada 1 con respecto 
                                        %al so en radianes
THETA1m2(j) = THETA1m2so(j) - THETAroll(j);
theta1m2(j) = deg2rad(THETA1m2(j)); %Angulo de entrada 1 en radianes

THETA2m2so(j) = atan2d(pyP(j) - pyM4(j), pxP(j) - pxM4(j));
if THETA2m2so(j)<0
    THETA2m2so(j) = 360 + THETA2m2so(j);
end
theta2m2so(j) = deg2rad(THETA2m2so(j)); %Angulo de entrada 2 con respecto 
                                        %al so en radianes
THETA2m2(j) = THETA2m2so(j) - THETAroll(j);
theta2m2(j) = deg2rad(THETA2m2(j)); %Angulo de entrada 2 en radianes

%Relaciones obtenidas de las ecuaciones
temp2 = ( r4m2*cos(theta2m2(j)) - r1m2*cos(theta1m2(j)) + 2*r3m2);
em2 = ( r1m2*sin(theta1m2(j))- r4m2*sin(theta2m2(j)) )/temp2;
temp3=(2*(r1m2*r3m2*cos(theta1m2(j))+r3m2*r4m2*cos(theta2m2(j)))-r1m2^2+r2m2^2+r4m2^2-r5m2^2);
fm2 = temp3/( 2*(r4m2*cos(theta2m2(j)) - r1m2*cos(theta1m2(j)) + 2*r3m2) );
%Coeficientes de la ecuacion cuadratica
dauxm2 = em2^2 + 1;
gm2 = 2*( em2*fm2 - em2*r1m2*cos(theta1m2(j)) + em2*r3m2 - r1m2*sin(theta1m2(j)) );
hm2=fm2^2- 2*fm2*(r1m2*cos(theta1m2(j))-r3m2)-2*r1m2*r3m2*cos(theta1m2(j))+r1m2^2+r3m2^2-r2m2^2;

pyP2(j) = ( - gm2 + sqrt( gm2.^2 - 4*dauxm2*hm2) )/ (2*dauxm2); 
pxP2(j) = em2*pyP2(j) + fm2;

%TRASLADO AL SISTEMA ORIGINAL DEL SEGUNDO MECANISMO DE 5 BARRAS
p2(j) = sqrt(pxP2(j)^2 + pyP2(j)^2);
THETA2p2(j) = atan2d(pyP2(j),pxP2(j));
theta2p2(j) = deg2rad(THETA2p2(j)); %Ang. del punto P2 respecto al sistema sec. en rad.

%Coordenadas del origen del sistema secundario con respecto al principal
pxAUX(j) = (pxS1(j) + pxM4(j))/2;
pyAUX(j) = (pyS1(j) + pyM4(j))/2;

%Coordenadas finales del punto P2 con respecto al sistema principal
pxP2so(j) = p2(j)*cos(theta2p2(j) + thetaroll(j)) + pxAUX(j);
pyP2so(j) = p2(j)*sin(theta2p2(j) + thetaroll(j)) + pyAUX(j);


%CALCULOS PARA EL SEGUNDO MECANISMO DE 4 BARRAS

THETA14B2(j) = atan2d( pyS2(j) - pyIFP(j) , pxS2(j) - pxIFP(j) );
theta14B2(j) = deg2rad(THETA14B2(j)); %Angulo de bancada mecanismo 4B 2 en radianes
THETA24B2(j) = atan2d( pyP2so(j) - pyS2(j) , pxP2so(j) - pxS2(j) );
if THETA24B2(j)<0
    THETA24B2(j) = 360 + THETA24B2(j);
end
theta24B2(j) = deg2rad(THETA24B2(j)); %Angulo de manivela mecanismo 4B 2 en radianes

k1m2=d2*cos(theta14B2(j))+a2*cos(theta24B2(j));
k2m2=d2*sin(theta14B2(j))+a2*sin(theta24B2(j));
k3m2=k1m2^2+k2m2^2+c2^2-b2^2;
k4m2=-c2^2+b2^2+k1m2^2+k2m2^2;
A1m2=-k3m2-2*k1m2*c2;
B1m2=4*k2m2*c2;
C1m2=2*k1m2*c2-k3m2;
A2m2=k4m2-2*k1m2*b2;
B2m2=4*k2m2*b2;
C2m2=2*k1m2*b2+k4m2;

%Calculo de posiciones angulares
tantheta31m2=(-B2m2+sqrt(B2m2^2-4*A2m2*C2m2))/(2*A2m2);
tantheta32m2=(-B2m2-sqrt(B2m2^2-4*A2m2*C2m2))/(2*A2m2);
tantheta41m2=(-B1m2+sqrt(B1m2^2-4*A1m2*C1m2))/(2*A1m2);
tantheta42m2=(-B1m2-sqrt(B1m2^2-4*A1m2*C1m2))/(2*A1m2);
theta3cm2(j)=2*atan(tantheta31m2);
theta3am2(j)=2*atan(tantheta32m2);
theta4cm2(j)=2*atan(tantheta41m2);
theta4am2(j)=2*atan(tantheta42m2);
THETA3cm2(j)=2*atan(tantheta31m2)*180/pi; % theta3 cruzado en grados
THETA3am2(j)=2*atan(tantheta32m2)*180/pi; % theta3 abierto en grados
if THETA3am2(j)<0
    THETA3am2(j) = 360 + THETA3am2(j);
end
THETA4cm2(j)=2*atan(tantheta41m2)*180/pi; %theta4 cruzado en grados
THETA4am2(j)=2*atan(tantheta42m2)*180/pi; %theta4 abierto en grados
if THETA4am2(j)<0
    THETA4am2(j) = 360 + THETA4am2(j);
end

% Obteniendo las coordenadas del punto P3
pxP3(j) = pxIFP(j) + c2*cos(deg2rad(THETA4am2(j)));
pyP3(j) = pyIFP(j) + c2*sin(deg2rad(THETA4am2(j)));

% Calculando la posicion de la falange medial
THETAfm(j) = THETA4am2(j)+ THETAauxfm;
thetafm(j)= deg2rad(THETAfm(j)); %Angulo de la falange medial en radianes

pxIFD(j) = fm*cos(thetafm(j)) + pxIFP(j);
pyIFD(j) = fm*sin(thetafm(j)) + pyIFP(j);

% Calculando la posicion de la falange distal
THETAfd(j) = THETAfm(j) + THETAauxfd;
thetafd(j)= deg2rad(THETAfd(j)); %Angulo de la falange distal en radianes
pxPF(j) = fd*cos(thetafd(j)) + pxIFD(j);
pyPF(j) = fd*sin(thetafd(j)) + pyIFD(j);

end


%A CONTINUACION GRAFICAMOS LOS RESULTADOS OBTENIDOS

% Para observar la grafica requerida solo basta con suprimir los
% comentarios. Algunas graficas poseen doble comentario en secciones
% especificas con la finalidad de alternar la forma de graficacion. La gran
% mayoria alterna entre graficar con respecto al tiempo o respecto a la
% manivela de entrada del mecanismo completo

% % Posicion en X del punto P
% figure(1)
% plot(tiempo,pxP)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxP)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada X (mm)')
% title('Coordenada X del punto P (Conf. +)')
% grid on


% % Posicion en Y del punto P
% figure(2)
% plot(tiempo,pyP)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyP)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada Y (mm)')
% title('Coordenada Y del punto P (Conf. +) ')
% grid on


% % Grafica manivela de entrada 5 barras con respecto al tiempo 
% figure(3)
% plot(tiempo,THETA2)
% xlabel('Tiempo (seg)')
% ylabel('\theta_{2} Manivela 5 Barras (grados)')
% title('Posicion de Manivela 5 barras respecto al tiempo')
% grid on



% % Grafica manivela del mecanismo inicial de 4 barras
% figure(5)
% plot(tiempo,THETA1)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA1)
% % xlabel('\theta_{2} Manivela 5 Barras (grados)')
% ylabel('\theta_{2} Manivela 4 Barras (grados)')
% title('Posicion de ambas manivelas iniciales')
% grid on

% % Angulo Theta3 cruzado del primer mecanismo de 4 barras
% figure(6)
% plot(THETA2,THETA3c)
% xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{3} cruzado (grados)')
% title('Posicion acoplador M4B (cruzado)')
% grid on

% % Angulo Theta3 abierto del primer mecanismo de 4 barras
% figure(7)
% plot(tiempo,THETA3a)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA3a)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{3} abierto (grados)')
% title('Posicion acoplador M4B (abierto)')
% grid on


% % Angulo Theta4 cruzado del primer mecanismo de 4 barras
% figure(8)
% plot(THETA2,THETA4c)
% xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{4} cruzado (grados)')
% title('Posicion balancin M4B (cruzado)')
% grid on

% % Angulo Theta4 abierto del primer mecanismo de 4 barras
% figure(9)
% plot(tiempo,THETA4a)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA4a)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{4} abierto (grados)')
% title('Posicion balancin M4B (abierto)')
% grid on

% % Angulo de la falange proximal
% figure(10)
% plot(tiempo,THETAfp)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETAfp)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{FP} (grados)')
% title('Posicion Falange Proximal')
% grid on

% % Posicion en X de la articulacion IFP
% figure(11)
% plot(tiempo,pxIFP)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxIFP)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada x_{IFP} (mm)')
% title('Coordenada X de la articulacion IFP')
% grid on

% % Posicion en Y de la articulacion IFP
% figure(12)
% plot(tiempo,pyIFP)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxIFP)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada y_{IFP} (mm)')
% title('Coordenada Y de la articulaciOn IFP')
% grid on

% % Posicion en X del punto S1
% figure(13)
% plot(tiempo,pxS1)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxS1)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada x_{S1} (mm)')
% title('Coordenada X del soporte S1')
% grid on

% % Posicion en Y del punto S1
% figure(14)
% plot(tiempo,pyS1)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyS1)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada y_{S1} (mm)')
% title('Coordenada Y del soporte S1')
% grid on

% % Posicion en X del punto S2
% figure(16)
% plot(tiempo,pxS2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxS2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada x_{S2} (mm)')
% title('Coordenada X del soporte S2')
% grid on

% % Posicion en Y del punto S2
% figure(17)
% plot(tiempo,pyS2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyS2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada y_{S2} (mm)')
% title('Coordenada Y del soporte S2')
% grid on

% % Posicion en X del punto M4
% figure(18)
% plot(tiempo,pxM4)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxM4)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada x_{M4} (mm)')
% title('Coordenada X del punto P_{M4}')
% grid on

% % Posicion en Y del punto M4
% figure(19)
% plot(tiempo,pyM4)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyM4)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada y_{M4} (mm)')
% title('Coordenada Y del punto P_{M4}')
% grid on

% % Angulo de giro respecto al sistema original del segundo mec. 5B
% figure(20)
% plot(tiempo,THETAroll)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETAroll)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{roll} giro respecto a ejes principales (grados)')
% title('Giro respecto a ejes principales')
% grid on

% % Angulo entrada 1 Mecanismo 2 5B
% figure(21)
% plot(tiempo,THETA1m2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA1m2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{1M2} Entrada 1 Mec. 5 Barras 2 (grados)')
% title('Entrada 1 Mecanismo 2 5 barras')
% grid on

% % Angulo entrada 2 Mecanismo 2 5B
% figure(22)
% plot(tiempo,THETA2m2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA2m2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{2M2} Entrada 2 Mec. 5 Barras 2 (grados)')
% title('Entrada 2 Mecanismo 2 5 barras')
% grid on

% % Posicion en X del punto P2 respecto al sistema original
% figure(23)
% plot(tiempo,pxP2so)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxP2so)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada X (mm) del punto P_{2}')
% title('Coordenada X del punto P_{2} (Conf. +)')
% grid on

% % Posicion en Y del punto P2 respecto al sistema original
% figure(24)
% plot(tiempo,pyP2so)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyP2so)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada Y (mm) del punto P_{2}')
% title('Coordenada Y del punto P_{2} (Conf. +)')
% grid on

% % Coordenadas del punto P
% figure(25)
% subplot(2,1,1);
% plot(tiempo,pxP)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxP)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada X (mm)')
% title('Coordenada X del punto P (Conf. +)')
% grid on
% subplot(2,1,2);
% plot(tiempo,pyP)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyP)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada Y (mm)')
% title('Coordenada Y del punto P (Conf. +) ')
% grid on

% % Coordenadas del punto P2
% figure(26)
% subplot(2,1,1);
% plot(tiempo,pxP2so)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxP2so)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada X (mm) del punto P_{2}')
% title('Coordenada X del punto P_{2} (Conf. +)')
% grid on
% subplot(2,1,2);
% plot(tiempo,pyP2so)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyP2so)
% % xlabel('\theta_{2} (grados)')
% ylabel('Coordenada Y (mm) del punto P_{2}')
% title('Coordenada Y del punto P_{2} (Conf. +)')
% grid on

% % Ruta de trazo de la articulacion IFP
% figure(27)
% plot(pxIFP,pyIFP,'DisplayName','Ruta Art. IFP')
% title('Ruta de la falange proximal')
% xlabel('Coordenadas en X (mm)')
% ylabel('Coordenadas en Y (mm)')
% grid on

% % Ruta de trazo de la articulacion IFP, punto P y punto P2
% figure(28)
% plot(pxIFP,pyIFP,'DisplayName','Ruta Art. IFP')
% hold on
% plot(pxP,pyP,'DisplayName','Ruta punto P')
% plot(pxP2so,pyP2so,'DisplayName','Ruta punto P2')
% scatter(-r3,-d,'filled','DisplayName','Art. MCF')
% scatter(0,0,'filled','DisplayName','Origen')
% title('Ruta de la falange proximal')
% xlabel('Coordenadas en X (mm)')
% ylabel('Coordenadas en Y (mm)')
% grid on
% legend({'Ruta Art. IFP','Ruta punto P','Ruta punto P2','Art. MCF', 'Origen'}
% ,'Location','southeast','Orientation','vertical');
% %drawnow
% hold off

% % PosiciOn en X del punto P3
% figure(29)
% plot(tiempo,pxP3)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pxP3)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada x_{P3} (mm)')
% title('Coordenada X del punto P3')
% grid on

% % Posicion en Y del punto P3
% figure(30)
% plot(tiempo,pyP3)
% xlabel('Tiempo (seg)')
% % plot(THETA2,pyS3)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('Coordenada y_{P3} (mm)')
% title('Coordenada Y del punto P3')
% grid on

% % Angulo de la bancada del segundo mecanismo 4B
% figure(31)
% plot(tiempo,THETA14B2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA14B2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{14B2} (grados)')
% title('Angulo de la bancada del segundo mecanismo 4B')
% grid on

% % Angulo de la manivela del segundo mecanismo 4B
% figure(32)
% plot(tiempo,THETA24B2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA24B2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{24B2} (grados)')
% title('Angulo de la manivela del segundo mecanismo 4B')
% grid on

% % Angulo del balancin del segundo mecanismo 4B
% figure(33)
% plot(tiempo,THETA4am2)
% xlabel('Tiempo (seg)')
% % plot(THETA2,THETA4am2)
% % xlabel('\theta_{2} Manivela 5B (grados)')
% ylabel('\theta_{4m2} (grados)')
% title('Angulo del balancin del segundo mecanismo 4B')
% grid on

% Ruta de trazo de la falange Proximal
figure(34)
plot(pxIFP,pyIFP,'DisplayName','Ruta Art. IFP')
title('Ruta de la falange Proximal')
xlabel('Coordenadas en X (mm)')
ylabel('Coordenadas en Y (mm)')
grid on
writematrix(pxIFP , 'proximal_I_mx.xlsx', 'Sheet', 1, 'Range', 'A1');
writematrix(pyIFP, 'proximal_I_my.xlsx', 'Sheet', 1, 'Range', 'A1');

% Ruta de trazo de la falange medial
figure(34)
plot(pxIFD,pyIFD,'DisplayName','Ruta Art. IFD')
title('Ruta de la falange medial')
xlabel('Coordenadas en X (mm)')
ylabel('Coordenadas en Y (mm)')
grid on
writematrix(pxIFD , 'medial_I_mx.xlsx', 'Sheet', 1, 'Range', 'A1');
writematrix(pyIFD, 'medial_I_my.xlsx', 'Sheet', 1, 'Range', 'A1');

% Ruta de trazo de la falange distal
figure(35)
plot(pxPF,pyPF,'DisplayName','Ruta Final')
title('Ruta de la falange distal')
xlabel('Coordenadas en X (mm)')
ylabel('Coordenadas en Y (mm)')
grid on
writematrix(pxPF , 'proximal_I_mx.xlsx', 'Sheet', 1, 'Range', 'A1');
writematrix(pyPF, 'proximal_I_my.xlsx', 'Sheet', 1, 'Range', 'A1');

% Rutas de trazo completas
figure(36)
plot(pxIFP,pyIFP,'DisplayName','Ruta Art. IFP')
hold on
plot(pxIFD,pyIFD,'DisplayName','Ruta Art. IFD')
plot(pxPF,pyPF,'DisplayName','Ruta punto Final')
scatter(-r3,-d,'filled','DisplayName','Art. MCF')
scatter(0,0,'filled','DisplayName','Origen')
title('trajectory of each phalanx')
xlabel('Coordenadas en X (mm)')
ylabel('Coordenadas en Y (mm)')
grid on
% legend({'Ruta Art. IFP','Ruta Art. IFD','Ruta punto Final','Art. MCF', 'Origen'}
% ,'Location','southeast','Orientation','vertical');
legend({'Ruta Art. IFP','Ruta Art. IFD','Ruta punto Final','Art. MCF', 'Origen'});
hold off
