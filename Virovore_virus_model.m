function dydt = Virovore_virus_model(~,y,a,e)

dydt = zeros(size(y));
R = y(1);
C = y(2);

dydt(1) = - a*R*C;
dydt(2) = e*a*R*C;