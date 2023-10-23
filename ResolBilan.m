function dydt = ResolBilan(t,y,p,I,nexp,tps,ts,Cste_Reacteur)
% ResolBilan Résout le bilan pour des conditions données
%
% Entrées:
%   t               - Temps actuel en sec
%   y               - Valeur actuelle de y
%   p               - Paramètres des modèle 1 à 5
%   p1              - Paramètres du modèle N°1 [a, f]
%   p2              - Paramètres du modèle N°2 [a, f, b2]
%   p3              - Paramètres du modèle N°3 [a, f, n]
%   p4              - Paramètres du modèle N°4 [a, n, b1]
%   p5              - Paramètres du modèle N°5 [a, b1, b2]
%   I               - Densité de flux en W/m^2 
%   nexp            - Numéro de l'expérience
%   tps             - Temps associés à l'expérience en sec
%   ts              - Longueur du shoulder initial en sec
%   Cste_Reacteur   - Géométrie du réacteur [s, vr, v] 
%
% Sortie:
%   dydt  - Dérivées de y par rapport au temps pour chacun des modèles
%   testés

% Constantes
s = Cste_Reacteur(1); %surface irradiée en m^2 - 0.015 
v_ir = Cste_Reacteur(2); % volume irradié en l - 0.3
v = Cste_Reacteur(3); % volume totale en l - 0.35

% Correction du volume en fonction du temps
maxIndex = length(tps);
v_index = find(t >= tps, 1, 'last');

if isempty(v_index)
   v = v;
else
   v = v-(0.003*v_index);
end

% Calcul des dydt
switch nexp
    case 1 % calcul uniquement avec le modèle N°1
        p1=p(1,:);
        dydt1(:,1) = Model1(p1, t, ts, y(1), s, v, v_ir, I);
        dydt2(:,1) = dydt1;
        dydt3(:,1) = dydt1;
        dydt4(:,1) = dydt1;
        dydt5(:,1) = dydt1;

    case 2 % calcul uniquement avec le modèle N°2
        p2=p(1,:);
        dydt2(:,1) = Model2(p2, t, ts, y(2), s, v, v_ir, I);
        dydt1(:,1) = dydt2;
        dydt3(:,1) = dydt2;
        dydt4(:,1) = dydt2;
        dydt5(:,1) = dydt2;

    case 3 % calcul uniquement avec le modèle N°3
        p3=p(1,:);
        dydt3(:,1) = Model3(p3, t, ts, y(3), s, v, v_ir, I);
        dydt1(:,1) = dydt3;
        dydt2(:,1) = dydt3;
        dydt4(:,1) = dydt3;
        dydt5(:,1) = dydt3;

    case 4 % calcul uniquement avec le modèle N°4
        p4=p(1,:);
        dydt4(:,1) = Model4(p4, t, ts, y(4), s, v, v_ir, I);
        dydt1(:,1) = dydt4;
        dydt2(:,1) = dydt4;
        dydt3(:,1) = dydt4;
        dydt5(:,1) = dydt4;

    case 5 % calcul uniquement avec le modèle N°5
        p5=p(1,:);
        dydt5(:,1) = Model5(p5, t, ts, y(5), s, v, v_ir, I);
        dydt1(:,1) = dydt5;
        dydt2(:,1) = dydt5;
        dydt3(:,1) = dydt5;
        dydt4(:,1) = dydt5;

    case 6 % calcul avec tous les modèles
        p1=p(1,:);
        dydt1(:,1) = Model1(p1, t, ts, y(1), s, v, v_ir, I);
        p2=p(2,:);
        dydt2(:,1) = Model2(p2, t, ts, y(2), s, v, v_ir, I);
        p3=p(3,:);
        dydt3(:,1) = Model3(p3, t, ts, y(3), s, v, v_ir, I);
        p4=p(4,:);
        dydt4(:,1) = Model4(p4, t, ts, y(4), s, v, v_ir, I);
        p5=p(5,:);
        dydt5(:,1) = Model5(p5, t, ts, y(5), s, v, v_ir, I);
    
    case 7 % calcul avec le model N°0
        p0=p;
        dydt1(:,1) = Model0(p0, t, ts, y(1), s, v, v_ir, I);
        dydt2(:,1) = dydt1;
        dydt3(:,1) = dydt1;
        dydt4(:,1) = dydt1;
        dydt5(:,1) = dydt1;

end
dydt=[dydt1; dydt2; dydt3; dydt4; dydt5];
end


function dydt = Model1(p1, t, ts, y, s, v, v_ir, I)
% Modèle N°1
a_1 = p1(1);
f_1 = p1(2);
if t < ts
    dydt(1,1)= 0;
% elseif y(1)<10000
%     dydt(1,1)=0;
else
    %dydt(1,1) = -(vr/v)*a_1*y(1)*I^f_1;
    dydt(1,1) = -(1/v)*a_1*(v_ir*y(1))*((s/v_ir)*I)^f_1;
end
end


function dydt = Model0(p0, t, ts, y, s, v, v_ir, I)
% Modèle N°0
a_0=p0(1);
f_0=p0(2);
n_0=p0(3);
b1_0=p0(4);
b2_0=p0(5);

if t < ts
    dydt(1,1)= 0;
% elseif y(1)<10000
%     dydt(1,1)=0;
else
    %dydt(1,1) = -(vr/v)*a_4*y(1)^n_4*(I/(b1_4*I+1));
    dydt(1,1) = -(1/v)*a_0*(v_ir*y(1))^n_0/(v_ir*y(1)*b2_0+1)*((s*I)^f_0/(b1_0*I*s+1));
end
end

function dydt = Model2(p2, t, ts, y, s, v, v_ir, I)
% Modèle N°2
a_2=p2(1);
f_2=p2(2);
b2_2=p2(3);

if t < ts
    dydt(1,1)= 0;
% elseif y(1)<10000
%     dydt(1,1)=0;
else
    %dydt(1,1) = -(vr/v)*a_2*I^f_2*(y(1)/(b2_2*y(1)+1));
    dydt(1,1) = -(1/v)*a_2*((s)*I)^f_2*(v_ir*y(1)/(b2_2*v_ir*y(1)+1));
end
end

function dydt = Model3(p3, t, ts, y, s, v, v_ir, I)
% Modèle N°3
a_3=p3(1);
f_3=p3(2);
n_3=p3(3);


% H=heaviside(t,ts)
% dydt(1,1) = -H*(1/v)*a_3*((s)*I)^f_3*(v_ir*y(1)*10)^n_3;

if t < ts
    dydt(1,1)= 0;
% elseif y(1)<10000
%     dydt(1,1)=0;
else
    %dydt(1,1) = -(vr/v)*a_3*I^f_3*y(1)^n_3;
    dydt(1,1) = -(1/v)*a_3*((s)*I)^f_3*(v_ir*y(1))^n_3;
    %dydt(1,1) = -a_3*((s/v_ir)*I)^f_3*(y(1))^n_3;
end
end

function dydt = Model4(p4, t, ts, y, s, v, v_ir, I)
% Modèle N°4
a_4=p4(1);
n_4=p4(2);
b1_4=p4(3);

if t < ts
    dydt(1,1)= 0;
% elseif y(1)<10000
%     dydt(1,1)=0;
else
    %dydt(1,1) = -(vr/v)*a_4*y(1)^n_4*(I/(b1_4*I+1));
    dydt(1,1) = -(1/v)*a_4*(v_ir*y(1))^n_4*((s)*I/(b1_4*(s)*I+1));
end
end

function dydt = Model5(p5, t, ts, y, s, v, v_ir, I)
% Modèle N°5
a_5=p5(1);
b1_5=p5(2);
b2_5=p5(3);

if t < ts
    dydt(1,1)= 0;
% elseif y(1)<10000
%     dydt(1,1)=0;
else
    %dydt(1,1) = -(vr/v)*a_5*(I/(b1_5*I+1))*(y(1)/(b2_5*y(1)+1));
    dydt(1,1) = -(1/v)*a_5*((s)*I/(b1_5*(s)*I+1))*(v_ir*y(1)/(b2_5*v_ir*y(1)+1));
end
end
