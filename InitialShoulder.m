function ts = InitialShoulder(yexp,texp,lim,I)
%
% Entrées:
%   texp            - Temps de mesure expérimentale en sec 
%   yexp            - Concentration bactérienne dans le réacteur bacth à l'instant t en cfu/100ml
%   lim             - la limite de variation toléré pour la définition du shoulder (en %) 
%
% Sortie:
%   ts              - La durée du Shoulder initial en sec
texp = texp.';
switch I
    case 10
        lim = (2+abs(I-35)/100)*lim;
    case 20
        lim = (2+abs(I-35)/100)*lim;
    case 30
        lim = (1+abs(I-35)/100)*lim;
    case 35
        lim = (1+abs(I-35)/100)*lim;
    case 45
        lim = (1.1+abs(I-35)/100)*lim;
end

y0 = yexp(1); % la concentration initale en bactérie dans le réacteur
i = 1;
lim;
for t = texp
    if abs(yexp(i)-y0)*100/y0 == lim
        ts(1,1) = t;
        i;
        break
    elseif abs(yexp(i)-y0)*100/y0 > lim
        j = i-1;
        %correction = (floor(abs(yexp(j)-y0)*100/y0)-lim)/(floor(abs(yexp(j)-y0)*100/y0)-floor(abs(yexp(i)-y0)*100/y0));
        %correction = (yexp(j)/y0-lim)/(yexp(j)/y0-yexp(i)/y0);
        %ts(1,1) = texp(j)+((t+texp(j))*correction);
        ts=(t+texp(j))/2;
        break
    else
        i = i+1;
    end
end

end