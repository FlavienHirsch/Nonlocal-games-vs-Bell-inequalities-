

Inputs = [4 4]

Outputs = [2 2]

nx = Inputs(1); ny = Inputs(2);

oa = Outputs(1); ob = Outputs(2);



	k = 1

	InhereM = In4422Full(k,:);

	M = ConvFull2Array(InhereM,Inputs,Outputs); %convert to table M(a,b,x,y)



	% extract the game

	[cStar, sol, exitflag, output] = MILP_is_M_det_transformable(M);


cStar


mu = sol.mu;

	V = sol.V;

	G = Game2Inequality(mu,V) ;

	canG = CanonicalizeInequality(G,Inputs,Outputs);

	canIn = CanonicalizeInequality(InhereM,Inputs,Outputs);



	check = norm(canG-canIn)


%%

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




f = -c ; 


ops=sdpsettings('verbose',1,'warning',1,'solver','mosek');

solvesdp(R,-f,ops);


double(f)



	G = Game2Inequality(double(mu),double(V)) ;

	canG = CanonicalizeInequality(G,Inputs,Outputs);

	canIn = CanonicalizeInequality(InhereM,Inputs,Outputs);



	check = norm(canG-canIn)
