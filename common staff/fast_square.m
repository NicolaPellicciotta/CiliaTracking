j=1:5e4;
square = zeros(length(j),1);
alg = square;
square = j.^2;
alg = (j-25).*100+(50-j).^2;
plot(j,square,'o');
hold on;
plot(j,alg,'*g');