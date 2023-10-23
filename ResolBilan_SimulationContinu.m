function dydt = ResolBilan_SimulationContinu(t,y,p,I_data,q_data,C0,Cste_Reacteur,n_model,tps_I, tps_q)
% ResolBilan Résout le bilan pour des conditions données
%
% Entrées:
%   t               - Temps actuel
%   y               - Valeur actuelle de y
%   p               - Paramètres du model
%   I_data          - Densité de flux en W/m^2 pour différents temps
%   q_data          - Débit en l/sec pour différents temps
%   C0              - Concentration de la solution d'alimentation en cfu/L
%   tps             - Temps associés à I et/ou q
%   Cste_Reacteur   - Géométrie du réacteur [s, vr, v] 
%   n_model         - le numéro du model à utilisé pour la simulation
%
% Sortie:
%   dydt            - Dérivée de y par rapport au temps

% Constantes
s = Cste_Reacteur(1); %surface irradiée en m^2 
v_ir = Cste_Reacteur(2); % volume irradié en l
v = Cste_Reacteur(3); % volume totale en l
C0 = C0; 

% Initialisation de la densité de flux en fonction du temps
if length(I_data) == 1
    I = I_data(1);
else
    maxIndex = length(tps_I);
    I_index = find(t >= tps_I, 1, 'last');
    
    if isempty(I_index)
       I = I_data(1);
    elseif I_index == maxIndex
       I = I_data(I_index);
    else
       I = I_data(I_index);
    end
end


% Initialisation de la valeur du ts
ts=0;
% if 10>=I_data(1)<20
%     ts=1300;
% elseif 20>=I_data(1)<35
%     ts=1100;
% elseif 35>=I_data(1)<45
%     ts=750;
% elseif I_data(1)>=45
%     ts=350;
% else
%     ts=0;
% end

% Initialisation du débit d'alimentation en fonction du temps
if length(q_data) == 1
    q=q_data;
else
    maxIndex = length(tps_q);
    q_index = find(t >= tps_q, 1, 'last');
    
    if isempty(q_index)
       q = q_data(1);
    elseif q_index == maxIndex
       q = q_data(q_index);
    else
       q = q_data(q_index);
    end
end


% Choix du model et calcul de dydt
switch n_model
    case 1 % calcul avec le model N°1
        p = p(1,:);
        dydt = Model1(p, t, y(1), s, v, v_ir, I, q, C0,ts);
    case 2 % calcul avec le model N°2
        p = p(2,:);
        dydt = Model2(p, t, y(1), s, v, v_ir, I, q, C0,ts);

    case 3 % calcul uniquement avec le modèle N°3
        p=p(3,:);
        dydt = Model3(p, t, y(1), s, v, v_ir, I, q, C0,ts);

    case 4 % calcul uniquement avec le modèle N°4
        p=p(4,:);
        dydt = Model4(p, t, y(1), s, v, v_ir, I, q, C0,ts);

    case 5 % calcul uniquement avec le modèle N°5
        p=p(1,:);
        dydt = Model5(p, t, y(1), s, v, v_ir, I, q, C0,ts);
end

end



function dydt = Model1(p, t, y, s, v, v_ir, I, q, C0,ts)
% Modèle N°1
a_1 = p(1);
f_1 = p(2);
if t < ts
    dydt(1,1)= 0;
else
    dydt(1,1) = (q/v) * (C0 - y(1))-(1/v)*a_1*(v_ir*y(1))*((s)*I)^f_1;
end
end



function dydt = Model2(p, t, y, s, v, v_ir, I, q, C0,ts)
% Modèle N°2
a_2=p(1);
f_2=p(2);
b2_2=p(3);
if t < ts
    dydt(1,1)= 0;
else
    dydt(1,1) = (q/v) * (C0 - y(1))-(1/v)*a_2*((s)*I)^f_2*(v_ir*y(1)/(b2_2*v_ir*y(1)+1));
end
end


function dydt = Model3(p, t, y, s, v, v_ir, I, q, C0,ts)
% Modèle N°3
a_3=p(1);
f_3=p(2);
n_3=p(3);
if t < ts
    dydt(1,1)= 0;
else
    dydt(1,1) = (q/v) * (C0 - y(1))-(1/v)*a_3*((s)*I)^f_3*(v_ir*y(1))^n_3;
end
end

function dydt = Model4(p, t, y, s, v, v_ir, I, q, C0,ts)
% Modèle N°4
a_4=p(1);
n_4=p(2);
b1_4=p(3);
if t < ts
    dydt(1,1)= 0;
else
    dydt(1,1) = (q/v) * (C0 - y(1))-(1/v)*a_4*(v_ir*y(1))^n_4*((s)*I/(b1_4*(s)*I+1));
end
end


function dydt = Model5(p, t, y, s, v, v_ir, I, q, C0,ts)
% Modèle N°5
a_5=p(1);
b1_5=p(2);
b2_5=p(3);
if t < ts
    dydt(1,1)= 0;
else
    dydt(1,1) = (q/v) * (C0 - y(1))-(1/v)*a_5*((s*I)/(b1_5*(s)*I+1))*(v_ir*y(1)/(b2_5*v_ir*y(1)+1));
end
end
