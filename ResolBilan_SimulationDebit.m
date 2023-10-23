function q = ResolBilan_SimulationDebit(t, y, p, I_data, n_model, C0, Cs, Cste_Reacteur, tps, q0)
% ResolBilan Résout le bilan pour des conditions données
%
% Entrées:
%   t               - Temps actuel
%   C0              - Concentration dans la solution d'alimentation en
%   cfu/100ml
%   Cs              - Concentration dans la solution de sortie en cfi/100ml
%   p               - Paramètres du modèle [a, f]
%   I               - Densité de flux en W/m^2 pour différents temps
%   n_model         - Numéro du model utilisé
%   tps             - Temps associés à I
%   Cste_Reacteur   - Géométrie du réacteur [s, vr, v] 
%
% Sortie:
%   q               - Débit en l/s

% Constantes
s = Cste_Reacteur(1); % surface irradiée en m^2 
v_ir = Cste_Reacteur(2); % volume irradié en l 
v = Cste_Reacteur(3); % volume totale en l  

i = 1;

% Choix du model et calcul de q
switch n_model
    case 1 % calcul avec le model N°1
        p = p(1,:);
        q = Model1(t, tps, I_data, i, q0, v_ir, s, Cs, v, p, C0);

    case 2 % calcul avec le model N°2
        p = p(2,:);
        q = Model2(t, tps, I_data, i, q0, v_ir, s, Cs, v, p, C0);

    case 3 % calcul uniquement avec le modèle N°3
        p=p(3,:);
        q = Model3(t, tps, I_data, i, q0, v_ir, s, Cs, v, p, C0);

    case 4 % calcul uniquement avec le modèle N°4
        p=p(4,:);
        q = Model4(t, tps, I_data, i, q0, v_ir, s, Cs, v, p, C0);

    case 5 % calcul uniquement avec le modèle N°5
        p = p(5,:);
        q = Model5(t, tps, I_data, i, q0, v_ir, s, Cs, v, p, C0);
end

end


function q = Model1(t, tps, I_data, i, q0, v_ir, s, y, v, p, C0)
% Modèle N°1
a_1 = p(1);
f_1 = p(2);
for time = t
    % Initialisation de la densité de flux en fonction du temps
    maxIndex = length(tps);
    I_index = find(time >= tps, 1, 'last');
    if isempty(I_index)
        I = I_data(1);
    elseif I_index == maxIndex
        I = I_data(I_index);
    else
        I = I_data(I_index);
    end
    
    if I == 0
        q(i,1) = q0; 
    else
        q(i,1) = (1/v_ir)*a_1*(v_ir*10*y(1))*((s)*I)^f_1*v_ir/(C0-y(1));
    end
    i=i+1;
end
end


function q = Model2(t, tps, I_data, i, q0, v_ir, s, y, v, p, C0)
% Modèle N°2
a_2=p(1);
f_2=p(2);
b2_2=p(3);
for time = t
    % Initialisation de la densité de flux en fonction du temps
    maxIndex = length(tps);
    I_index = find(time >= tps, 1, 'last');
    if isempty(I_index)
        I = I_data(1);
    elseif I_index == maxIndex
        I = I_data(I_index);
    else
        I = I_data(I_index);
    end
    
    if I == 0
        q(i,1) = q0; 
    else
        q(i,1) = (1/v_ir)*a_2*((s)*I)^f_2*(v_ir*y(1)*10/(b2_2*v_ir*10*y(1)+1))*v_ir/(C0-y(1));
    end
    i=i+1;
end
end


function q = Model3(t, tps, I_data, i, q0, v_ir, s, y, v, p, C0)
% Modèle N°3
a_3=p(1);
f_3=p(2);
n_3=p(3);
for time = t
    % Initialisation de la densité de flux en fonction du temps
    maxIndex = length(tps);
    I_index = find(time >= tps, 1, 'last');
    if isempty(I_index)
        I = I_data(1);
    elseif I_index == maxIndex
        I = I_data(I_index);
    else
        I = I_data(I_index);
    end
    
    if I == 0
        q(i,1) = q0; 
    else
        q(i,1) = (1/v_ir)*a_3*((s)*I)^f_3*(v_ir*y(1)*10)^n_3*v_ir/(C0-y(1));
    end
    i=i+1;
end
end


function q = Model4(t, tps, I_data, i, q0, v_ir, s, y, v, p, C0)
% Modèle N°4
a_4=p(1);
n_4=p(2);
b1_4=p(3);
for time = t
    % Initialisation de la densité de flux en fonction du temps
    maxIndex = length(tps);
    I_index = find(time >= tps, 1, 'last');
    if isempty(I_index)
        I = I_data(1);
    elseif I_index == maxIndex
        I = I_data(I_index);
    else
        I = I_data(I_index);
    end
    
    if I == 0
        q(i,1) = q0; 
    else
        q(i,1) = (1/v_ir)*a_4*(v_ir*10*y(1))^n_4*((s)*I/(b1_4*(s)*I+1))*v_ir/(C0-y(1));
    end
    i=i+1;
end
end


function q = Model5(t, tps, I_data, i, q0, v_ir, s, y, v, p, C0)
% Modèle N°5
a_5=p(1);
b1_5=p(2);
b2_5=p(3);
for time = t
    % Initialisation de la densité de flux en fonction du temps
    maxIndex = length(tps);
    I_index = find(time >= tps, 1, 'last');
    if isempty(I_index)
        I = I_data(1);
    elseif I_index == maxIndex
        I = I_data(I_index);
    else
        I = I_data(I_index);
    end
    
    if I == 0
        q(i,1) = q0; 
    else
        q(i,1) = (1/v_ir)*a_5*((s)*I/(b1_5*(s)*I+1))*(v_ir*y(1)*10/(b2_5*v_ir*y(1)*10+1))*v_ir/(C0-y(1));
    end
    i=i+1;
end
end
