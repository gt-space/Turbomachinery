function [d_2,H,eta_H] = EffeciencyTester(Q,H)
n=23000;
rho = 804; %kg/m^3
g = 9.81;
P(1) = rho*g*H*Q;
n_q(1) = ( n * ( sqrt ( Q ) )) / ( H(1) ^ ( 3 / 4 ) )
if Q <= 1
    a = 1;
elseif Q > 1 
    a = .5;
end
m = 0.08 * a * ((1/Q)^(.15)) * ((45/n_q)^(0.06));
eta_H =  1 - (0.055*((1/Q)^m)) - .2*((.26 - log10(n_q/25))^2)*((1/Q)^(.1));

num(1)=1;
for i = 2:100
    num(i)=i;
    H(i) = H(1)/eta_H;
    P(i) = rho*g*H(i)*Q;
    n_q(i) = ( n * ( sqrt ( Q ) )) / ( H(i) ^ ( 3 / 4 ) );
    m = 0.08 * a * ((1/Q)^(.15)) * ((45/n_q(i))^(0.06));
    eta_H =  1 - (0.055*((1/Q)^m)) - .2*((.26 - log10(n_q(i)/25))^2)*((1/Q)^(.1));
end

H = H(length(H));
d_2 = 84.6/n * sqrt(H/1.089);
eta_H = eta_H(length(eta_H));



end

