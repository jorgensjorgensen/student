load CIDNP;
x=CIDNP(:,2);
y=CIDNP(:,3);
z=CIDNP(:,4);
plot(x,y,'b',x,abs(z),'r')
