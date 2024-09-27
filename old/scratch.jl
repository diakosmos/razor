

(x,DM,DeM,DoM) =cheb.chebdif(25,2);
De = DeM[:,:,1];
Do = DoM[:,:,1];
xp = x[1:13]
up = 1.0./(xp.-1.1).^4
A = [Do  Do*0.0;De*0.0 De];
B = [Do*0.0  -Do; -De -diagm(ones(13))];
AA = A*diagm([up;up])
