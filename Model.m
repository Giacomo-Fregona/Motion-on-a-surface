%-------------------------------------------------------------------------
%--------------------------MOTION ON A SURFACE----------------------------
%-------------------------------------------------------------------------

%           Giacomo Fregona - Corso Modelli Matematici - UniTS



%CHOOSING PHYSICAL CONSTANTS
m=1; %la massa del punto
g=9.81; %costante di gravitazione universale


%CHOOSING STEP AND DIMENSION OF THE PLOTTED SURFACE
h=input('passo h = '); %passo della griglia
area_stampa=[-10,10]; %definiamo l'area della superficie rappresentata


%CHOOSING THE SURFACE
%scegliamo la funzione f che descrive la superficie, l'intervallo di tempo tspan = [0 t] in cui seguiremo il moto del corpo, il vettore q_iniziale (posizione e velocità iniziali), 
syms x y real
%f = 0.1*(y+x); tspan=[0 5]; q_iniziale=[0;9;2;-1]; %piano inclinato h=0.3
%f = -(sqrt(100-x^2-y^2)); tspan=[0 3]; q_iniziale=[0;0;0;0];area_stampa=[-8 8];%emisfera di raggio 5 con punto fermo al centro
%f = ((y+x)^2)/10; tspan=[0 7]; q_iniziale=[-4;8;1;-1];area_stampa=[-9 9];%paraboloide h=0.03
%f = 0.9*sin(x)+cos(y); tspan=[0 8]; q_iniziale=[-6;-3;3;2]; area_stampa=[-7 1 -12 1];% h=0.1
%f = sin(x)+cos(y); tspan=[0 6]; q_iniziale=[-6;-3;3;2]; area_stampa=[-7 2 -6 1];% h=0.1
f = -(sqrt(100-x^2-y^2)); tspan=[0 10]; q_iniziale=[8;0;0;6];%emisfera di %raggio 10 h=0.1


%SYMBOLIC CALCULATIONS
%calcoliamo il gradiente in ogni punto
gradiente_f=[diff(f,x),diff(f,y)];
%calcoliamo i due vettori u,v che generano il piano tangente
u=[1;0;0]+[0;0;1]*gradiente_f(1);
u=u/norm(u);
v=[0;1;0]+[0;0;1]*gradiente_f(2);
v=v/norm(v);
%calcoliamo il versore normale alla superficie
N=cross(u,v);
N=N/norm(N);
%costruiamo il vettore vel che rappresenta la velocità del punto
syms xp yp real
vel=[xp;yp;xp*gradiente_f(1)+yp*gradiente_f(2)];
%definiamo la curva gamma
syms t real;
gamma=[x+xp*t;y+yp*t;subs(f,[x,y],[x+xp*t,y+yp*t])];
%calcoliamo la curvatura (con segno)
k=subs( norm(cross(diff(gamma,t,1),diff(gamma,t,2))/(norm(diff(gamma,t,1) )^3+eps(0))) ,t,0);%*sign(subs(diff(gamma(3),t,2),t,0));
k=k*(cross(subs(cross(diff(gamma,t),diff(diff(gamma,t),t))/(norm(cross(diff(gamma,t),diff(diff(gamma,t),t)))+eps(0)),t,0),vel /(norm(vel)+eps(0)))'*N);
%calcoliamo il vettore forza normale risultante
Fn=N*(norm(vel)^2)*k;
%calcoliamo il vettore P (peso parallelo)
P=g*m*[0;0;-1]-(g*m*[0;0;-1]'*N)*N;


%PLOTTING THE SURFACE AND INITIAL STATE
ezmesh(f,area_stampa) %stampa il grafico di f nel quadrato  
pause(5)
hold on
plot3(q_iniziale(1),q_iniziale(2),subs(f,[x,y],[q_iniziale(1),q_iniziale(2)]),'o','Color','b','MarkerSize',5,'MarkerFaceColor','b') %stampa il punto di partenza della curva
quiver3(q_iniziale(1),q_iniziale(2),subs(f,[x,y],[q_iniziale(1),q_iniziale(2)]),q_iniziale(3),q_iniziale(4),subs(xp*gradiente_f(1)+yp*gradiente_f(2),[x,y,xp,yp],q_iniziale'),'Color',[0 .7 0],'LineWidth',5)
pause(0.01)


%NUMERICAL CALCULATIONS
%Risoluzione con RK4
Fx=@(t,x_,y_,xp_,yp_) double(subs(Fn(1)+P(1),[x,y,xp,yp,],[x_,y_,xp_,yp_]))/m;
Fy=@(t,x_,y_,xp_,yp_) double(subs(Fn(2)+P(2),[x,y,xp,yp,],[x_,y_,xp_,yp_]))/m;
N=ceil((tspan(2)-tspan(1))/h);%numero di punti della griglia
t=zeros(1,N);
X=zeros(1,N);
XP=zeros(1,N);
Y=zeros(1,N);
YP=zeros(1,N);
Z=zeros(1,N);
ZP=zeros(1,N);
t(1)=tspan(1);
X(1)=q_iniziale(1);
XP(1)=q_iniziale(3);
Y(1)=q_iniziale(2);
YP(1)=q_iniziale(4);
Z(1)=subs(f,[x,y],[X(1),Y(1)]);
% calcolo della traiettoria con il metodo RK4
for n= 1:(N-1)
    k1x=Fx(t(n),X(n),Y(n),XP(n),YP(n));
    k1y=Fy(t(n),X(n),Y(n),XP(n),YP(n));
    k2x=Fx(t(n)+h/2,X(n)+XP(n)*h/2,Y(n)+YP(n)*h/2,XP(n)+h/2*k1x,YP(n)+h/2*k1y);
    k2y=Fy(t(n)+h/2,X(n)+XP(n)*h/2,Y(n)+YP(n)*h/2,XP(n)+h/2*k1x,YP(n)+h/2*k1y);
    k3x=Fx(t(n)+h/2,X(n)+XP(n)*h/2+(h^2)*k1x/4,Y(n)+YP(n)*h/2+(h^2)*k1y/4,XP(n)+h/2*k2x,YP(n)+h/2*k2y);
    k3y=Fy(t(n)+h/2,X(n)+XP(n)*h/2+(h^2)*k1x/4,Y(n)+YP(n)*h/2+(h^2)*k1y/4,XP(n)+h/2*k2x,YP(n)+h/2*k2y);
    k4x=Fx(t(n)+h,X(n)+XP(n)*h+(h^2)*k2x/2,Y(n)+YP(n)*h+(h^2)*k2y/2,XP(n)+h*k3x,YP(n)+h*k3y);
    k4y=Fy(t(n)+h,X(n)+XP(n)*h+(h^2)*k2x/2,Y(n)+YP(n)*h+(h^2)*k2y/2,XP(n)+h*k3x,YP(n)+h*k3y);
    XP(n+1)=XP(n)+h/6*(k1x+2*k2x+2*k3x+k4x);
    YP(n+1)=YP(n)+h/6*(k1y+2*k2y+2*k3y+k4y);
    X(n+1)=X(n)+h*XP(n)+h^2/6*(k1x+k2x+k3x);
    Y(n+1)=Y(n)+h*YP(n)+h^2/6*(k1y+k2y+k3y);
    t(n+1)=t(n)+h;
    Z(n+1)=subs(f,[x,y],[X(n+1),Y(n+1)]);
    plot3(X(n+1),Y(n+1),Z(n+1),'o','Color','b','MarkerSize',5,'MarkerFaceColor','b') %stampa il punto appena calcolato
    pause(h)
end


%FINAL VISUALIZATION
plot3(X(1:N),Y(1:N),Z(1:N),'Color',[0 .7 0],'LineWidth',5); %stampa il grafico della curva
hold off