function [cStar,G,mu,V,r,q,s] =  MILP_is_M_game_equivalent(M) 

%MILP_is_M_game_equivalent  Check if functional 'M' can be transformed
%into a deterministic game using renormalization, shift and NS degrees of
%freedom 
%
%   Problem:
%
%      c* = max_{G, V, mu, s, r, q, c}  c
%
%       s.t.   V(a,b,x,y) in {0,1}
%              0 <= mu(x,y) <= 1
%              G(a,b,x,y) := c*M(a,b,x,y)  + r(a,x,y) + q(b,x,y) + s(x,y)
%              sum_y r(a,x,y) = 0,    sum_x q(b,x,y) = 0
%              0 <= G(a,b,x,y) <= V(a,b,x,y)
%              mu(x,y) - (1 - V(a,b,x,y)) <= G(a,b,x,y) <= mu(x,y) + (1 - V(a,b,x,y))
%
%   The last two pairs of inequalities are the standard McCormick / big-M
%   linearization of  G = V * mu  with M = 1 (valid because mu is in [0,1]).
%
%
%   Inputs:
%     M    oa x ob x nx x ny array, representing the coefficients of a Bell
%     functional (in scenario (nx,ny,oa,ob) )
%
%   Outputs:
%     cStar     Optimal value of c 
%		G		Game functional ( G(a,b,x,y) = mu(x,y)*V(a,b,x,y) )
%		mu		Prior of the game
%		V		Predicate of the game
%		r		B->A NS freedom 
%		q		A->B NS freedom 
%		s		block shifts
%
%   Example
%     M = randn(3,3,2,2);
%     [cStar, G] = MILP_is_M_game_equivalent(M);
%	  If cStar > 10^-8 , G should be a game representation of M, otherwise M is likely game-inequivalent 



%Extract scenario:

oa = size(M,1); ob = size(M,2) ;

nx = size(M,3); ny = size(M,4) ;


%%%%SDP


G=sdpvar(oa,ob,nx,ny,'full','real');

V=binvar(oa,ob,nx,ny,'full','real');

mu=sdpvar(nx,ny,'full','real');

r = sdpvar(oa,nx,ny,'full','real');
q = sdpvar(ob,nx,ny,'full','real');
s = sdpvar(nx,ny,'full','real');

c = sdpvar(1,1); 


R = (mu(:) >= 0) + (sum(mu(:)) == 1); %normalized prior

for x=1:nx
	for y=1:ny

		for a=1:oa
			for b=1:ob


				R = R + ( G(a,b,x,y) == c*M(a,b,x,y) + r(a,x,y) + q(b,x,y) + s(x,y) );   % G = c*M + r + q + s

			end
		end

	end
end


R = R + (G(:) >= 0) + (V(:) >= G(:)); % 0 <= G <= V

R = R + (sum(r,3) == 0) + (sum(q,2) == 0) ; %sum_y r = 0, sum_x q = 0 (NS-orthogonal variables)



mu_clone = repelem(mu(:), oa*ob);     % [mu1 mu1 ... mu2 mu2 ...]

R = R + ( G(:) >= mu_clone(:) - (1-V(:)) ) + ( mu_clone(:) + (1-V(:)) >= G(:) ) ; % mu - (1-V) <= G <= mu + (1-V)






ops=sdpsettings('verbose',1,'warning',1,'solver','mosek');

optimize(R,-c,ops);


cStar = double(c);

G = double(G); mu = double(mu); V=double(V);

s=double(s); q=double(q); r=double(r);


end
