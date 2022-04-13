% hill function
function H = hill(X,X0,lambda,n)
% H+ 
H1 = (X^n)/(X0^n+X^n);
% H-
H2 = (X0^n)/(X0^n+X^n);
 
% Hill function =  H- plus (lambda)* H+
H= H2 + lambda*H1;