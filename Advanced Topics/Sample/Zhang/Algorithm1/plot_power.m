% Compute SINR in Each SC Macrocell
for n=1 : macrocell(k).N
  for t=1 : 20
     power(n,t)= 0;
  end 
end

for t=1 : 20
   for n=1 : macrocell(k).N
      for i=1 : macrocell(k).M
          power(n,t)= power(n,t) + macrocell(k).P(i,n,t);
      end
  end
  %sinr(n,t)=sinr(n,t)/macrocell(k).M;
end
linespec = {'-+b','-+r','-+g','-+m', '-sb','-sr','-sg','-sm','-ob','-or','-og','-om',};
figure;
hold;
for n=1 : 6
  plot(1:1:t, power(n,:),linespec{n});
end
hold off;
xlabel('Time(Iteration)')
ylabel('Average Power')
legend('SC l','SC 2','SC 3','SC 4','SC 5','SC 6','SC 7','SC 8');
grid
figure;



%Compute SINR in Each SC Femtocell
for c=1 : macrocell(k).nFcel
    for n=1 : macrocell(k).Fcel(c).N
        for t=1 : 20
            power(n,t)= 0;
        end
    end
    
    for t=1 : 20
        for n=1 : macrocell(k).Fcel(c).N
            for i=1 : macrocell(k).Fcel(c).M -1
                power(n,t)= power(n,t) + macrocell(k).Fcel(c).P(i,n,t);
            end
        end
        %sinr(n,t)=sinr(n,t)/macrocell(k).M;
    end
    
    hold;
    for n=1 : 6
        plot(1:1:t, power(n,:),linespec{n}); 
    end
    hold off;
    xlabel('Time(Iteration)')
    ylabel('Average Power')
    legend('SC l','SC 2','SC 3','SC 4','SC 5','SC 6','SC 7','SC 8');
    grid
    
    figure;
end