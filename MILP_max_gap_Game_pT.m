


function [cStar, G, mu, V , beta] = MILP_max_gap_Game_pT(pT,D)

%MILP_max_gap_Game_pT Optimizes the Q-C gap of behavior 'pT' over all
%nonlocal deterministic games, where C is computed over all (optional)
%strategies 'D'


%More precisely, solve the following optimization problem:
%
%   max_{G,μ,V,β}   sum_{abxy} pT(a,b,x,y) G(a,b,x,y) - β
%
%   s.t.   β ≥ sum_{abxy} G(a,b,x,y) D_det(a,b,x,y,k)         ∀k
%		   0 ≤ μ(x,y) ≤ 1 , sum_xy μ(x,y) = 1
%          V(a,b,x,y) ∈ {0,1}
%          0 ≤ G(a,b,x,y) ≤ V(a,b,x,y)
%          μ(x,y) - (1-V(a,b,x,y)) ≤ G(a,b,x,y) ≤ μ(x,y) + (1-V(a,b,x,y))
%
%   The last two inequality blocks linearize G = μ·V (valid since
%   μ ∈ [0,1])
%
%   Inputs
%     pT     oa x ob x nx x ny array, representing a behaviour
%     in bipartite Bell scenario (nx,ny,oa,ob)
%Optional:
%	  D		oa x ob x nx x ny x N array, representing 'N' local deterministic
%	  strageties, indexed by the last index (that is, D(a,b,x,y,k) is a behavior indexed like 'pT', for all values of k)
%	Default: D is the complete set of deterministic strategies for the
%	corresponding (oa,ob,nx,ny) Bell scenario
%
%   Outputs:
%     cStar     Optimal value of the quantum-to-classical of behaviour pT
%     (over all deterministic nonlocal games)
%		G		Game functional ( G(a,b,x,y) = mu(x,y)*V(a,b,x,y) )
%		mu		Prior of the game
%		V		Predicate of the game
%		beta	Local value of the game
%
%   Example
%     [cStar, G] = MILP_max_gap_Game_pT(pT,D);
%	  If cStar > 10^-8 , G should be a valid nonlocal game which pT violates, otherwise pT is likely game-classical




%Extract scenario:

oa = size(pT,1); ob = size(pT,2) ;

nx = size(pT,3); ny = size(pT,4) ;



%%%%SDP


G=sdpvar(oa,ob,nx,ny,'full','real');

V=binvar(oa,ob,nx,ny,'full','real');

mu=sdpvar(nx,ny,'full','real');

beta=sdpvar(1,1);


R = (mu(:) >= 0) + (sum(mu(:)) == 1); %normalized prior



R = R + (G(:) >= 0) + (V(:) >= G(:)); % 0 <= G <= V


if nargin < 2 || isempty(D)
	Dmat = local_det_strategies(oa,ob,nx,ny);
else
	Dmat = sparse(reshape(D, [], size(D,5)));
end
R = R + (beta >= Dmat' * G(:));    % beta >= G*D_k  for all k


mu_clone = repelem(mu(:), oa*ob);     % [mu1 mu1 ... mu2 mu2 ...]

R = R + ( G(:) >= mu_clone(:) - (1-V(:)) ) + ( mu_clone(:) + (1-V(:)) >= G(:) ) ; % mu - (1-V) <= G <= mu + (1-V)



f = G(:)' * pT(:) - beta ;


ops=sdpsettings('verbose',1,'warning',1,'solver','mosek');

optimize(R,-f,ops);


cStar = double(f);

G = double(G); mu = double(mu); V=double(V);

beta = double(beta);



%helper function

	function Dmat = local_det_strategies(oa,ob,nx,ny)

		% Builds all local deterministic behaviours of the
		%bipartite Bell scenario (oa,ob,nx,ny) as columns of a sparse 0/1 matrix
		%
		%Strategy k = (ia,ib) assigns a = lamA(ia,x), b = lamB(ib,y), so
		%D(a,b,x,y) = 1 iff a = lamA(x) and b = lamB(y).
		%
		%   Inputs
		%     oa,ob  number of outputs of Alice / Bob
		%     nx,ny  number of inputs of Alice / Bob
		%
		%   Output
		%     Dmat   sparse (oa*ob*nx*ny) x (oa^nx * ob^ny) matrix. Column k is
		%            D(:,:,:,:,k)(:), with the same linear ordering as pT(:) and G(:).
		%            Each column has exactly nx*ny ones.
		%
		%   Examples
		%     Dmat = local_det_strategies(2,2,2,2);   % CHSH: 16 columns
		%     D    = reshape(full(Dmat),2,2,2,2,[]);  % back to the 5-D format
		%     beta = max(Dmat' * G(:));               % local value of a game G


		NA = oa^nx; NB = ob^ny; N = NA*NB; d = oa*ob*nx*ny; %number of local and global strategies

		lamA = mod(floor((0:NA-1)' ./ oa.^(0:nx-1)), oa) + 1;   % NA x nx, outputs 1..oa
		lamB = mod(floor((0:NB-1)' ./ ob.^(0:ny-1)), ob) + 1;   % NB x ny

		[IA,IB] = ndgrid(1:NA,1:NB); IA = IA(:); IB = IB(:);    % k = IA + NA*(IB-1)
		rows = zeros(N, nx*ny); c = 0;
		for y = 1:ny
			for x = 1:nx
				c = c + 1;
				rows(:,c) = sub2ind([oa ob nx ny], lamA(IA,x), lamB(IB,y), x*ones(N,1), y*ones(N,1));
			end
		end

		cols = repmat((1:N)', 1, nx*ny);
		Dmat = sparse(rows(:), cols(:), 1, d, N);

	end



end



