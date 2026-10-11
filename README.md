# Bell nonlocality does not imply violation of nonlocal games

Code and data accompanying the article

> R. Faleiro, L. E. A. Porto, M. T. Quintino, and F. Hirsch,
> *Bell nonlocality does not imply violation of nonlocal games*,
> [arXiv:2610.06799](https://arxiv.org/abs/2610.06799) (2026).

**In short.** Nonlocal games with a deterministic predicate, $G_{ab|xy} = \mu(x,y)\,V(a,b|x,y)$ with $V \in \{0,1\}$, are not complete witnesses of Bell nonlocality: there exist nonlocal behaviours (including quantum ones) that violate no such game. This holds in the scenario $(4,4,2,2)$ and in every bipartite scenario $(X,Y,A,B)$ with $X,Y \ge 2$ and $A,B \ge 4$.

This repository provides:

1. a MILP that decides whether a Bell functional is equivalent to a nonlocal game (up to positive scaling, block shifts and no-signalling transformations) — Appendix B;
2. a MILP that computes the maximal quantum-to-classical gap of a behaviour over all nonlocal games, i.e. tests membership in the game-classical polytope $\mathcal{G}$ — Eq. (11) and Appendix G;
3. explicit game representations of all game-equivalent facets of the local polytope in $(4,4,2,2)$ (173 of 174 classes) and $(2,2,4,4)$ (31 of 33 known classes);
4. the three game-inequivalent facets $J_\star$ (in $(4,4,2,2)$), $F^1_\star$ and $F^2_\star$ (in $(2,2,4,4)$), and a nonlocal-game representative of $J_\star$ after input lifting to $(5,5,2,2)$;
5. the script reproducing the critical $\mathcal{G}$-visibility computation of Appendix G.

---

## Requirements

- MATLAB (tested on R2023b)
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (free academic licence available); any MILP solver supported by YALMIP should work after changing `'solver','mosek'` in `sdpsettings`.

## Repository contents

| File | Description |
|---|---|
| `MILP_is_M_game_equivalent.m` | Tests whether a Bell functional `M` is game-equivalent (Appendix B, Algorithm 1). |
| `MILP_max_gap_Game_pT.m` | Maximises $\langle G, p\rangle - \omega_L(G)$ over all deterministic-predicate games $G$ (Eq. (11)). $c^\star = 0$ ⇔ $p \in \mathcal{G}$. |
| `QuantumViolation_F_star_and_G_visibility.m` | Builds the quantum behaviour of Appendix F violating $F^1_\star$ and $F^2_\star$, mixes it with the barycentre of the face of $\mathcal{L}$ defined by $F^1_\star$, and evaluates the game gap (Appendix G). |
| `List_of_173_FacetInequalities_4422_pfull.mat` | The 173 game-equivalent facet classes of $(4,4,2,2)$, full-probability notation. |
| `Game_representation_173_FacetInequalities_4422.mat` | A nonlocal-game representative $(\mu, V)$ for each of them. |
| `List_of_31_FacetInequalities_2244_pfull.mat` | The 31 game-equivalent facet classes of $(2,2,4,4)$, full-probability notation. |
| `Game_representation_31_FacetInequalities_2244.mat` | A nonlocal-game representative $(\mu, V)$ for each of them. |
| `Jstar.mat` | The game-inequivalent facet $J_\star$ of $(4,4,2,2)$ (Eqs. (7), (C1)). |
| `F_star_1.mat`, `F_star_2.mat` | The game-inequivalent facets $F^1_\star$, $F^2_\star$ of $(2,2,4,4)$ (Eqs. (8)–(9), (D1)–(D2)). |
| `Lifted_J_star_game_5522.mat` | A nonlocal game $(\mu, V)$ in $(5,5,2,2)$ equivalent to the input-lifted $J_\star$. |
| `Vertices_F_star_1.mat` | The 56 local deterministic behaviours saturating $F^1_\star$ (used to build $p_\star$ in Appendix G). |

## Conventions and data format

All behaviours and functionals are stored as 4-D arrays `P(a,b,x,y)` of size `oa × ob × nx × ny` (MATLAB 1-based indices: output `a` ↔ `a-1` in the paper). A Bell scenario $(X,Y,A,B)$ corresponds to `nx = X, ny = Y, oa = A, ob = B`.

**Facet lists.** `In4422_list_173_game_facets_pfull` (173 × 64) and `In2244_list_31_game_facets_pfull` (31 × 64), `int16`. Each row is one Bell functional `M`, with the inequality read as $\langle M, p\rangle \le \beta_L(M)$ (local bounds are not stored). Row entries are ordered with `x` slowest, then `y`, `a`, and `b` fastest. To obtain the `oa × ob × nx × ny` array:

```matlab
S   = load('List_of_173_FacetInequalities_4422_pfull.mat');
oa = 2; ob = 2; nx = 4; ny = 4;
row = double(S.In4422_list_173_game_facets_pfull(i,:));
M   = permute(reshape(row, [ob oa ny nx]), [2 1 4 3]);   % M(a,b,x,y)
```

**Game representations.** Cell arrays `V_173_facets_4422`, `mu_173_facets_4422` (resp. `V_31_facet_2244`, `mu_31_facet_2244`), with entry `i` corresponding to row `i` of the facet list:

- `V{i}`: `oa × ob × nx × ny` 0/1 predicate;
- `mu{i}`: `nx × ny` prior. For $(2,2,4,4)$ it is normalised; for $(4,4,2,2)$ it is given up to normalisation (and some entries are stored as `uint8`), so use `mu = double(mu)/sum(double(mu(:)))`.

They satisfy $\mu(x,y)V(a,b|x,y) = c\,M_{ab|xy} + r_{a|x}(y) + s_{b|y}(x) + d_{xy}$ with $c>0$, $\sum_y r_{a|x}(y) = 0$, $\sum_x s_{b|y}(x) = 0$ (Definition 3 and Eq. (A9)).

**Game-inequivalent facets.** `Jstar` (`2 × 2 × 4 × 4`), `F_star_1` and `F_star_2` (`4 × 4 × 2 × 2`), `int16`, indexed `M(a,b,x,y)`. They differ from the matrices printed in Eqs. (C1), (D1), (D2) by a block shift of $+1$ on the $(x,y) = (0,0)$ block, so their local bounds are shifted accordingly:

| Variable | Stored inequality | Paper's form |
|---|---|---|
| `Jstar` | $\langle J_\star, p\rangle \le 2$ | $\le 1$ |
| `F_star_1`, `F_star_2` | $\langle F^i_\star, p\rangle \le 1$ | $\le 0$ |

**Lifted game in $(5,5,2,2)$.** `V` (`2 × 2 × 5 × 5`, 0/1) and `mu` (`5 × 5`, normalised). The game $\mu(x,y)V(a,b|x,y)$ is equivalent, in the sense of Definition 3, to $J_\star$ lifted to five inputs per party (zero coefficients for the fifth inputs `x = 5`, `y = 5`). Note that the variables are named `V` and `mu`: load them into a struct (`S = load(...)`) to avoid overwriting other variables.

**Vertices.** `Vertices_F_star1` is `4 × 4 × 2 × 2 × 56` (`uint8`); `Vertices_F_star1(:,:,:,:,k)` is the $k$-th deterministic behaviour saturating $F^1_\star$.

## Usage

**Is a Bell functional game-equivalent?**

```matlab
M = randn(3,3,2,2);                          % any functional M(a,b,x,y)
[cStar, G, mu, V] = MILP_is_M_game_equivalent(M);
% cStar > 1e-8 : G = mu.*V is a game equivalent to M
% cStar ~ 0    : M is (numerically) game-inequivalent -- not a proof, see Appendix B
```

**Does a behaviour violate some nonlocal game?**

```matlab
[cStar, G, mu, V, beta] = MILP_max_gap_Game_pT(p);   % p(a,b,x,y); all local deterministic strategies used by default
% cStar > 0 : p violates the game G (score > beta)
% cStar = 0 : p is game-classical, p in G (up to solver precision)
```

The number of local deterministic strategies is `oa^nx * ob^ny`, so this becomes expensive quickly.

**Reproduce the $\mathcal{G}$-visibility of Appendix G.**

```matlab
QuantumViolation_F_star_and_G_visibility
```

With `vis = 0.042` the behaviour violates a nonlocal game; with `vis = 0.041` the MILP returns $c^\star = 0$ (up to solver precision), giving a quantum, Bell-nonlocal, game-classical behaviour.

## Citation

```bibtex
@article{faleiro2026bell,
  title   = {Bell nonlocality does not imply violation of nonlocal games},
  author  = {Faleiro, Ricardo and Porto, Lucas E. A. and Quintino, Marco T{\'u}lio and Hirsch, Flavien},
  journal = {arXiv preprint arXiv:2610.06799},
  year    = {2026}
}
```

## Acknowledgements

The facet lists of $(4,4,2,2)$ and $(2,2,4,4)$ come from E. Zambrini Cruzeiro (see Refs. [23, 24] of the article).

## License

MIT — see [`LICENSE`](LICENSE).
