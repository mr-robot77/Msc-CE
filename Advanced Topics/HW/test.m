function mnodes=test(n,x,y)
rx = unifrnd(1:n, x); 
ry = unifrnd(1:n, y);
for j=1 : n
    mnodes(j)= complex(rx(j), ry(j));
end
%rx(6)=50;
%ry(6)=50;
bsx=250;
bsy=250;
plot(bsx,bsy,'ro') ;
hold;
plot(rx, ry, '*');

