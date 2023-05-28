function mnodes=test(n,x,y)
rx = unifrnd(1:n, x); 
ry = unifrnd(1:n, y);
for j=1 : n
    mnodes(j)= complex(rx(j), ry(j));
end
bsx=x/2;
bsy=y/2;
%plot(bsx,bsy,'ro') ;
%hold;
%plot(rx, ry, '*');



