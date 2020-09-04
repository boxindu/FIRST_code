function [fk, gk] = GoodnessF(x,A,Q,I,S,a,b)
fk = norm((x*A*x' - Q),'fro')^2-a*trace(S*x')+b*norm((x*x' - I),'fro')^2;
gk = (2*x*A*(x')*x*A'+2*x*A'*(x')*x*A-2*Q*x*A'-2*Q'*x*A)-a.*S+b.*(4*x*(x')*x-4*x);
end
