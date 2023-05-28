 clear;

Total_Sp_eff=0;
for k=1 : 1
    macrocell(k) = Macro(5,10,1000,1000,10,3);
    macrocell(k)=PathGain(macrocell(k));
    macrocell(k)=Update(macrocell(k),1);
    for i=2 : 30
        macrocell(k)=Update(macrocell(k),0);
    end
    for j=1 : 10
     Total_Sp_eff = Total_Sp_eff+ macrocell(k).Fcel(j).SP_Eff(30);
    end
end
Total_Sp_eff=Total_Sp_eff/10;
disp(Total_Sp_eff);
 


