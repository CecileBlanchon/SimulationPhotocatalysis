function r = VitesseInactivation_Moyenne(q,y0,y)
r_instantane = q.*(y0-y);
r = mean(r_instantane)*3600;
end

  