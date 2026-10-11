%% Construct a simple violation of F_\star^1 and F_\star^2

%Scenario 

Inputs = [2 2]

Outputs = [4 4]



%Pauli matrices:

X = [0 1;1 0]; Y = [0 -1i;1i 0]; Z = [1 0;0 -1];


rho = 1/2 * [1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]; %Two-qubit maximally entangled state

th = [ pi/2, 0]; %A's angles

ph = [ pi/4 , 3*pi/4 ]; %B's angles

P_x = [sin(th(1)) 0  cos(th(1)) ; sin(th(2)) 0  cos(th(2)) ]; %A's Bloch vectors

P_y = [sin(ph(1)) 0  cos(ph(1)) ; sin(ph(2)) 0  cos(ph(2)) ]; %B's Bloch vectors

%POVMs:

%2-outcome projective measurements:
PA = zeros(2,2,4,2);
for x=1:2
PA(:,:,1,x) = .5*(eye(2) + P_x(x,1)*X + P_x(x,2)*Y + P_x(x,3)*Z ); % PA_{0|x} is the projector with Bloch vector P_x(x,:)
PA(:,:,2,x) = eye(2) - PA(:,:,1,x); % PA_{1|x} = 1 - PA_{0|x} 
end

PB= zeros(2,2,4,2);
for y=1:2
PB(:,:,1,y) = .5*(eye(2) + P_y(y,1)*X + P_y(y,2)*Y + P_y(y,3)*Z ); % PB_{0|y} is the projector with Bloch vector P_y(y,:)
PB(:,:,2,y) = eye(2) - PB(:,:,1,y); % PB_{1|x} = 1 - PB_{0|x} 
end

%set POVM elements to 0
A = zeros(2,2,4,2); B = zeros(2,2,4,2);

%Non-zero elements will be 1 and 4
A(:,:,1,:) = PA(:,:,1,:) ;
A(:,:,4,:) = PA(:,:,2,:) ; 

B(:,:,1,:) = PB(:,:,1,:) ;
B(:,:,4,:) = PB(:,:,2,:) ; 


pQ = zeros(4,4,2,2);

for x=1:2
	for y=1:2

		for a=1:4
			for b=1:4

				pQ(a,b,x,y) = trace( kron(A(:,:,a,x),B(:,:,b,y)) * rho) ;  %p(a,b|x,y) = Tr( [A_{a|x} \otimes B_{b|y}] * rho) 
			end
		end

	end
end





%% G visibility 

load('Vertices_F_star_1.mat')


pN = sum(Vertices_F_star1,5)/size(Vertices_F_star1,5);  %Equal mixture of all vertices of F_\star^1

vis = .041 %choose weight of 'pQ'

pT = vis*pQ+(1-vis)*pN; %pT = vis*pQ + (1-vis)*1/N*sum_k V_k , where V_k are vertices of F_star^1


[cStar,G] = MILP_max_gap_Game_pT(pT);  %compute the max violation of pT over all (deterministic) nonlocal games


gap_pT = cStar  %if >10^-8 one can check that pT violates a (deterministic) nonlocal games, otherwise, pT is likely in G (the set of game-classical behaviours) 



