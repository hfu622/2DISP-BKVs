# 2DISP-BKVs

This project identifies candidate best known (utilization) values (BKVs) for 25 Nesting instances among the audited publications on non-exact algorithms for the **Two-Dimensional Irregular Strip Packing (2DISP)** problem. It compares published results with values obtained using author code or our implementations. Sources with unconfirmed inputs or configurations are marked separately.

## 1. Dataset Source

The data come from the [official ESICUP 2d_irregular dataset](https://github.com/ESICUP/datasets/tree/main/2d_irregular). The [25 official XML files](official_esicup/xml) preserve the contours, quantities, initial orientations, and allowed rotations. `shapes2.xml` is an alias for the upstream `blaz/blaz.xml`; its contents stay unchanged.

## 2. Algorithms and Assignments

| Selected algorithm | Instances | Code source |
|---|---|---|
| Sparrow [1] | Albano, Dagli, Fu, Jakobs1/2, Mao, Marques, Shapes0/1/2, Shirts, Swim, Trousers | Author's MIT-licensed Rust source, pinned version |
| GCS [2] | Dighe1/2 | Our dedicated Rust/Jagua implementation |
| F057 supplementary SVGnest fixed-order vector baseline [3] | Han, Poly3b/4b | Author's JavaScript wrapper and pinned SVGnest version; source retrieval tool with hash verification |
| SNCPH [4] | Poly2a/3a/4a/5a | Our dedicated Rust/Jagua implementation |
| 3PM [5] | Poly1a | Our Rust/Jagua implementation with an exact native DBM backend and a Z3 rational continuous backend |
| CHS2 [6] | Poly2b | Our Rust/Jagua implementation with CBC |
| SemiDiscreteBLF [7] | Poly5b | Our dedicated Rust/Jagua implementation |

1. Gardeyn, Vanden Berghe, and Wauters (2025). [An open-source heuristic to reboot 2D nesting research](https://arxiv.org/abs/2509.13329v4). [Author's Sparrow code](https://github.com/JeroenGar/sparrow).
2. Elkeran (2013). [A new approach for sheet nesting problem using guided cuckoo search and pairwise clustering](https://doi.org/10.1016/j.ejor.2013.06.020).
3. (2026). [Raster-Convolution Search for Two-Dimensional Irregular Packing with Curved and Holed Parts](https://doi.org/10.1109/AANN71145.2026.11709270). We use the [supplementary fixed-order vector baseline](https://github.com/Hermione-HKX/Raster-Convolution-Search-and-CTM-Bench/blob/c7aaa7e8c31db15f4f7df7673b2cecfc362b21cd/results/svgnest_public_xml.csv), rather than the proposed raster method. [Author's wrapper source](https://github.com/Hermione-HKX/Raster-Convolution-Search-and-CTM-Bench) and [SVGnest](https://github.com/Jack000/SVGnest).
4. Baldacci et al. (2014). [Algorithms for nesting with defects](https://doi.org/10.1016/j.dam.2012.03.026). We use its strip-packing heuristic branch.
5. (2016). [A model-based heuristic for the irregular strip packing problem](https://doi.org/10.1590/0101-7438.2016.036.03.0447).
6. Antonio Martinez Sykora (2013). [Nesting problems: exact and heuristic algorithms](https://www.uv.es/marsyan/docs/thesis.pdf). The CHS2 constructive heuristic is described in Chapter 6 of this doctoral thesis.
7. Chehrazad, Roose, and Wauters (2022). [A fast and scalable bottom-left-fill algorithm to solve nesting problems using a semi-discrete representation](https://doi.org/10.1016/j.ejor.2021.10.043).

## 3. Reimplementation Methods

Selected methods without available author code are implemented in Rust, using [Jagua](https://github.com/JeroenGar/jagua-rs) 1.1.0 for the geometry infrastructure. NFP/IFR geometry, contact candidates, price updates, insertion and local compaction follow the definitions in the respective papers, preserving concavities and degenerate point and line contacts. Compactness is assessed against each method's objective, neighborhood, or discrete envelope: SNCPH checks local phi improvements; semi-discrete BLF seeks the earliest column and lowest feasible position; 3PM/CHS2 retain the model's minimum-length objective and distinguish time-limited feasibility from proven subproblem optimality. Geometric feasibility alone does not establish algorithmic fidelity or global compactness.

3PM uses native DBM discrete search and [Z3](https://github.com/Z3Prover/z3) rational continuous optimization, with [CBC](https://github.com/coin-or/Cbc) as a fallback for models the specialized backend cannot recognize. CHS2 uses CBC in place of CPLEX. Exact optimization subproblems within an overall heuristic are distinct from solving the complete 2DISP problem with an exact algorithm. Layouts are independently checked against the original official contours using Shapely/GEOS; rational geometry checks are also available for degenerate fixed-0° cases. Shrinking pieces or substantially changing the strip width is not an acceptable way to obtain a utilization value.

The repository includes [reimplementation source](reproductions/src), [functional tests](reproductions/tests), [author source and licenses](vendor), and [actual configurations and decisions for unspecified parameters](config/reproduction_decisions.json). Licenses apply separately to each source: Sparrow is MIT-licensed and Jagua is MPL-2.0-licensed. The F057 wrapper's license is unconfirmed, so its body is not redistributed; a pinned source list with SHA256 checks is provided instead.

## 4. Published References and Measured Values

Measured utilization uses the official nominal strip width; higher is better. Lengths are expressed in the paper comparison units; lower is better. Raw XML scaling units must not be mixed with those units. The table shows the best measured value under the disclosed configuration. For the seven instances whose selected source reports only length, utilization is derived from the official XML and marked with a superscript footnote; the original length is retained. Differences are measured minus reference, in percentage points (pp) for utilization and in L for length. F057 supplementary values include both rounded utilization and length, so very small differences are largely attributable to display precision.

† The author's original input or some configuration details are unavailable; the published value is a conditional reference. ‡ An accepted small numerical variant; see the limitations below.

| Instance | Algorithm | Published / Derived Reference | Measured Utilization % | Measured Length L | Difference |
|---|---|---|---:|---:|---|
| albano | Sparrow | 89.82% | 89.785369 | 9695.863000 | -0.034631 pp |
| dagli | Sparrow | 90.17% | 89.825012 | 56.303917 | -0.344988 pp |
| dighe1† | GCS | 100.00% | 96.863871 | 103.237667 | -3.136129 pp |
| dighe2† | GCS | 100.00% | 99.737820 | 100.262870 | -0.262180 pp |
| fu‡ | Sparrow | 92.41% | 92.400787 | 30.843893 | -0.009213 pp |
| han | SVGNest-F057-fixed-order | 63.60%; L=53.131 | 63.603421 | 53.130952 | +0.003421 pp; -0.000048 L |
| jakobs1‡ | Sparrow | 89.26% | 89.090488 | 11.000052 | -0.169512 pp |
| jakobs2‡ | Sparrow | 87.73% | 86.071262 | 22.423280 | -1.658738 pp |
| mao | Sparrow | 86.87% | 86.829851 | 1697.535400 | -0.040149 pp |
| marques | Sparrow | 92.02% | 91.212726 | 75.837090 | -0.807274 pp |
| poly1a† | 3PM | 64.873418%[^1] (L=15.8) | 66.386209 | 15.439954 | +1.512791 pp; -0.360046 L |
| poly2a† | SNCPH | 73.792811%[^1] (L=27.75) | 67.368523 | 30.396251 | -6.424288 pp; +2.646251 L |
| poly2b† | CHS2 | 69.363497%[^1] (L=32.6) | 63.342667 | 35.698686 | -6.020830 pp; +3.098686 L |
| poly3a† | SNCPH | 73.748517%[^1] (L=41.65) | 72.701201 | 42.250000 | -1.047316 pp; +0.600000 L |
| poly3b | SVGNest-F057-fixed-order | 64.64%; L=47.182 | 64.643931 | 47.181536 | +0.003931 pp; -0.000464 L |
| poly4a† | SNCPH | 73.792811%[^1] (L=55.5) | 71.026150 | 57.661876 | -2.766661 pp; +2.161876 L |
| poly4b | SVGNest-F057-fixed-order | 66.49%; L=58.166 | 66.491027 | 58.165743 | +0.001027 pp; -0.000257 L |
| poly5a† | SNCPH | 73.926011%[^1] (L=69.25) | 70.063058 | 73.068124 | -3.862953 pp; +3.818124 L |
| poly5b† | SemiDiscreteBLF | 66.043165%[^1] (L=69.5) | 64.647887 | 71.000000 | -1.395278 pp; +1.500000 L |
| shapes0‡ | Sparrow | 69.98% | 69.366562 | 57.520510 | -0.613438 pp |
| shapes1‡ | Sparrow | 76.73% | 76.706854 | 52.016212 | -0.023146 pp |
| shapes2‡ | Sparrow | 86.23% | 86.009366 | 25.113544 | -0.220634 pp |
| shirts | Sparrow | 90.92% | 89.991151 | 60.005900 | -0.928849 pp |
| swim‡ | Sparrow | 79.83% | 78.915916 | 5604.745600 | -0.914084 pp |
| trousers | Sparrow | 92.62% | 92.333556 | 235.888020 | -0.286444 pp |

[^1]: These utilization values are calculated by this project from the official XML, rather than reported directly in the original publication: $\eta = 100 A / (W L)$, where $A$ is the total area of all requested pieces, $W$ is the fixed strip width, and $L$ is the published length. All three quantities are first expressed in the paper comparison units. For Poly2a/3a/4a/5a, XML lengths are divided by 100 and XML areas by 10,000; Poly1a, Poly2b, and Poly5b require no such scaling. Because the author's original input or some settings are unavailable, these remain conditional references and do not establish identical author inputs. Differences use unrounded values before being displayed to six decimal places.

The [machine-readable per-instance summary](data/results_summary.json) retains the original source metrics, measured means, worst values, sample standard deviations where applicable, and source links. It contains scientific result summaries without operational records.

## 5. Limitations, unavailable information, and adopted configurations

The assignments identify candidate BKVs among audited sources; they do not prove the best value across the complete literature. A measured improvement is not automatically a new record for an identical published instance. Sixteen instances meet this project's input comparison rules: nine exact matches and seven quantified small variants. For Dighe1/2, Poly1a, Poly2a/b, Poly3a, Poly4a, and Poly5a/b, the author's original inputs or some configurations are unavailable. These nine cases retain conditional assignments and use official XML. Strict equivalence to the author's experiment is not claimed when original inputs, pre-rotations, or settings cannot be confirmed. Allowed Poly rotations follow the official XML.

Official Shirts contains 99 pieces; the original local copy contained 85. The official version is used here. Official Jakobs2 contains complete endpoints; the small contour variants affecting Jakobs2/Swim retain comparability limitations. Han's official XML total area is 1960, whereas another XLS/PDF version has area 1970; these versions are not mixed. Sparrow uses a small 0.01% strip-width expansion for Fu, Jakobs1/2, and Shapes0/1/2. Utilization in this table still uses the official nominal width, so these values should not be treated as records under a strictly unchanged width and zero geometric error.

| Algorithm | Unspecified or Unavailable Information | Adopted Configuration |
|---|---|---|
| SNCPH | K, MaxIt, pattern counts, and per-instance raster/panel settings are unspecified; price increments and decrements are only required to be positive | K=20, MaxIt=2, 10 patterns per configuration, price=penalty=0.1, initial price 1; finite panel length 4H; raster grids 512×128, 1024×256, and 2048×512. SNH1/2 and SNH3 initialized by either of them are each combined with PH1/2/3, giving 12 configurations. These are disclosed project choices, not recovered unpublished author controls |
| GCS | The Lévy exponent, boundary handling, and replacement rounding are unspecified; the prose and pseudocode conflict on the expansion base | Dedicated zero-clustering Dighe branch; 15 nests, 30 generations, abandonment 0.25, stagnation 200; Mantegna β=1.5, clipping, and replacement of the worst four nests per generation. Following the prose, expand the current length by 1% and shrink by 4%, retaining the best feasible layout independently; 600 seconds per trial including NFP preparation |
| 3PM | Author reference vertices and s in Figure 6 are undefined; the batch-selection prose conflicts with Algorithm 4; CPLEX is unavailable | First actual XML vertex; s=gmin, initial neighborhood side 4gmin; follow Algorithm 4; g=2→1→0.5, md=3000; integer supporting planes from original vertices, exact DBM/Z3 backends, and CBC fallback. Subproblem optimality does not prove global packing optimality |
| CHS2 | The CPLEX environment, recommended CBC settings, and warm-start details are unavailable | n_ini=12, 100 seconds per MIP; serial CBC2.10.13 with integer preprocessing and LP presolve disabled and solver random seeds fixed at 1. A feasible MIP start is submitted to the same model; existing slices are inherited while coordinates remain free, and the complete model is still validated. Solver substitution can affect time-limited search trajectories |
| SemiDiscreteBLF | Tie-breaking and some geometry implementation details are unspecified | R=0.5; one-sided extension of zero-length vertex intervals with L/R labels preserved; interior bridging disabled. Equal bounding-box areas are ordered stably by input index. All intervals are rechecked after each upward move; rounding stalls advance to the next representable y. If original-contour validation rejects a candidate, contact search continues in the same column |
| Sparrow | No nightly date or historical transitive dependencies are fixed | Pinned author source; nightly-2026-01-15, SIMD, release, target-cpu=native, only_final_svg, and locked dependencies; three workers, 960 seconds of exploration and 240 seconds of compression. Hardware, operating-system, and dependency differences can still affect the amount of search within a fixed time |
| F057 supplementary vector baseline | The wrapper license is unconfirmed; the supplementary baseline must not be confused with the proposed method | Pinned author wrapper and SVGnest versions, fixed order, official angles, and one placement-worker call; no added genetic search or generic post-compaction |

Geometric checks preserve the original contours and coordinates. The linear tolerance is `1e-6 × max(strip width, 1)` and the area tolerance is `1e-8 × max(total piece area, 1)`. Original MIP constraints use a row-scaled residual check of `1e-7`. Acceptance of small numerical errors differs from strictly zero overlap and zero boundary excess; passing a tolerance-based check does not establish exact geometric feasibility. Substantial width changes, piece scaling, or geometry errors beyond tolerance do not qualify as valid BKVs.

Under the disclosed CHS2 configuration, measured utilization ranges from 25.158157% to 63.342667%, with a mean of 55.857504% and a sample standard deviation of 9.084133 percentage points. A gap against the conditional published length remains; feasibility does not establish equivalent reproduction performance. Unavailable information and the actual choices are recorded in the [configuration file](config/reproduction_decisions.json).
