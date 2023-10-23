function E = Optimisation(p,tempsexp1,tempsexp2,tempsexp3,tempsexp4,tempsexp5,ydata1,ydata2,ydata3,ydata4,ydata5,Cin,I,ts,Cste_Reacteur,nexp)

p=p(1,:)
if  (p(1)<0)||(p(2)<0)||(p(3)<0);
    E = 10000;
else
    p(1);
    p(2);
    p(3);

   replicas = [1 2 3];

    for i = replicas
        y0 = [Cin(1,i) Cin(1,i) Cin(1,i) Cin(1,i) Cin(1,i)];
        tps=tempsexp1(:,i);
        [t,ysim] = ode45(@ResolBilan_SimulationContinu,tempsexp1(:,i),y0,[],p,I(1),nexp,tps,ts(1,i),Cste_Reacteur);%Résolution équa diff
        y1(:,i) = ysim(:,1);
    end

    for i = replicas
        y0 = [Cin(2,i) Cin(2,i) Cin(2,i) Cin(2,i) Cin(2,i)];
        tps=tempsexp2(:,i);
        [t,ysim] = ode45(@ResolBilan,tempsexp2(:,i),y0,[],p,I(2),nexp,tps,ts(2,i),Cste_Reacteur);%Résolution équa diff
        y2(:,i) = ysim(:,1);
    end

    for i = replicas
        y0 = [Cin(3,i) Cin(3,i) Cin(3,i) Cin(3,i) Cin(3,i)];
        tps=tempsexp3(:,i);
        [t,ysim] = ode45(@ResolBilan,tempsexp3(:,i),y0,[],p,I(3),nexp,tps,ts(3,i),Cste_Reacteur);%Résolution équa diff
        y3(:,i) = ysim(:,1);
    end

    for i = replicas
        y0 = [Cin(4,i) Cin(4,i) Cin(4,i) Cin(4,i) Cin(4,i)];
        tps=tempsexp4(:,i);
        [t,ysim] = ode45(@ResolBilan,tempsexp4(:,i),y0,[],p,I(4),nexp,tps,ts(4,i),Cste_Reacteur);%Résolution équa diff
        y4(:,i) = ysim(:,1);
    end

    for i = replicas
        y0 = [Cin(5,i) Cin(5,i) Cin(5,i) Cin(5,i) Cin(5,i)];
        tps=tempsexp5(:,i);
        [t,ysim] = ode45(@ResolBilan,tempsexp5(:,i),y0,[],p,I(5),nexp,tps,ts(5,i),Cste_Reacteur);%Résolution équa diff
        y5(:,i) = ysim(:,1);
    end


     Y1=log10(y1);
     Y2=log10(y2);
     Y3=log10(y3);
     Y4=log10(y4);
     Y5=log10(y5);

     ydata1=log10(ydata1);
     ydata2=log10(ydata2);
     ydata3=log10(ydata3);
     ydata4=log10(ydata4);
     ydata5=log10(ydata5);
    
     
    E1_1=mean((abs((ydata1(:,1)-Y1(:,1))./ydata1(:,1))*100));
    E2_1=mean((abs((ydata2(:,1)-Y2(:,1))./ydata2(:,1))*100));
    E3_1=mean((abs((ydata3(:,1)-Y3(:,1))./ydata3(:,1))*100));
    E4_1=mean((abs((ydata4(:,1)-Y4(:,1))./ydata4(:,1))*100));
    E5_1=mean((abs((ydata5(:,1)-Y5(:,1))./ydata5(:,1))*100));

    E1_2=mean((abs((ydata1(:,2)-Y1(:,2))./ydata1(:,2))*100));
    E2_2=mean((abs((ydata2(:,2)-Y2(:,2))./ydata2(:,2))*100));
    E3_2=mean((abs((ydata3(:,2)-Y3(:,2))./ydata3(:,2))*100));
    E4_2=mean((abs((ydata4(:,2)-Y4(:,2))./ydata4(:,2))*100));
    E5_2=mean((abs((ydata5(:,2)-Y5(:,2))./ydata5(:,2))*100));

    E1_3=mean((abs((ydata1(:,3)-Y1(:,3))./ydata1(:,3))*100));
    E2_3=mean((abs((ydata2(:,3)-Y2(:,3))./ydata2(:,3))*100));
    E3_3=mean((abs((ydata3(:,3)-Y3(:,3))./ydata3(:,3))*100));
    E4_3=mean((abs((ydata4(:,3)-Y4(:,3))./ydata4(:,3))*100));
    E5_3=mean((abs((ydata5(:,3)-Y5(:,3))./ydata5(:,3))*100));

    

    E=E1_1+E2_1+E3_1+E4_1+E5_1+E1_2+E2_2+E3_2+E4_2+E5_2+E1_3+E2_3+E3_3+E4_3+E5_3;
end

end
