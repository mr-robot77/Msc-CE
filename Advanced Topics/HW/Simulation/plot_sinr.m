% Compute SINR in Each SC Macrocell
for n=1 : macrocell(k).N
  for t=1 : 20
     sinr(n,t)= 0;
  end 
end

for t=1 : 20
   for n=1 : macrocell(k).N
      for i=1 : macrocell(k).M
          sinr(n,t)= sinr(n,t) + macrocell(k).sinr(i,n,t);
      end
  end
  %sinr(n,t)=sinr(n,t)/macrocell(k).M;
end
linespec = {'-+b','-+r','-+g','-+m', '-sk','-sr','-sg','-sm','-ob','-or','-og','-om',};
figure;
hold;
for n=1 : 6
  plot(1:1:t, sinr(n,:),linespec{n});
end
hold off;
xlabel('Time(Iteration)')
ylabel('Average SINR')
legend('SC l','SC 2','SC 3','SC 4','SC 5','SC 6','SC 7','SC 8');
grid

figure;



%Compute SINR in Each SC Femtocell
for c=1 : macrocell(k).nScel
    for n=1 : macrocell(k).Scel(c).N
        for t=1 : 20
            sinr(n,t)= 0;
        end
    end
    
    for t=1 : 20
        for n=1 : macrocell(k).Scel(c).N
            for i=1 : macrocell(k).Scel(c).M -1
                sinr(n,t)= sinr(n,t) + macrocell(k).Scel(c).sinr(i,n,t);
            end
        end
        %sinr(n,t)=sinr(n,t)/macrocell(k).M; 
    end
    
    hold;
    for n=1 : 6
        plot(1:1:t, sinr(n,:),linespec{n}); 
    end
    hold off;
    xlabel('Time(Iteration)')
    ylabel('Average SINR')
    legend('SC l','SC 2','SC 3','SC 4','SC 5','SC 6','SC 7','SC 8');
    grid
    
    figure;
end