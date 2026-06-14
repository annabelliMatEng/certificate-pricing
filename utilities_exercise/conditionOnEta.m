function [c, c_eq]=conditionOnEta(P, alpha)

sigma=P(1);
kappa=P(2);
eta=P(3);

c=-eta-(1-alpha)/(kappa*(sigma^2));
c_eq=[];

end