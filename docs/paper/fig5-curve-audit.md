# Full-profile comparison of the mouse and human proximal tubule

## Three-row Fig. 5 draft

**a–c. Ortholog profiles and positional filtering.** The first heatmap pair shows all 3,912 screened
orthologs in order of their mouse positional bias. Each row is one gene, shown at the same row position
in both heatmaps. The panels display equal-weight section means of centered log₂(CP10k + 0.5) on the
shared 0.19–0.875 COAST interval. The filter retains a gene when its fitted range reaches 1 log₂ unit
in both sections of at least one group. Of the 278 retained genes, 195 meet this rule only in mouse,
30 only in human and 53 in both. Four k-means clusters group the full human-minus-mouse difference
curves of these 278 genes after scaling each gene by the largest within-section SD across its four fits.
The heatmap cluster sizes are 130, 72, 42 and 34. Their labels describe the contrast curves. For example,
an early-positive and late-negative contrast can arise when a late mouse gradient is shallower in human.

**d–e. Functional programs and shape agreement.** The left grid counts membership of each selected
local annotation in the four clusters. The right grid restricts the same annotated sets to directly
comparable genes with reproducible twofold variation in both groups, using the shape criteria below.
Of the 53 comparable genes, 37 have conserved shapes, nine are intermediate and seven diverge. Cell tints
represent the fraction of that row in the corresponding cluster or shape category. The annotations can
overlap. For programs with only one or a few comparable members, the panel shows their coverage and
individual outcomes. All 14 comparable members of the selected solute-carrier transport set are conserved.

**f–m. All seven divergent profiles and a conserved reference.** The fitted peaks are earlier in both
human sections for SLC17A3, UGT3A1, APOC3 and MME, and later for TMEM252, PCK1 and MSRB1. DPEP1 is
the conserved reference for MME. Lines represent section fits and points show fixed-width position-bin
counts normalized by non-mitochondrial library size. Species comparisons use within-section centered
rates. The MME–DPEP1 correlations use observed position-bin rates within each section. Throughout, the
human comparison uses two cortical sections from one donor.

Run `python analysis/paper/fig5_curve_audit.py` to prepare the fixed-width counts and then
`python analysis/paper/fig5_three_row_figure.py` to export this draft. The exports and the row order and
program counts are in `results/paper_figures/fig5/three_row/`.

## Global comparison

Of 3,912 orthologs examined, 175 had reproducible positional variation in both mouse kidneys and both human
cortical sections. Within these 175, 97 retained the same profile shape (55%), 48 diverged (27%), and 30 were
intermediate or sensitive to the comparison (17%). Among the remaining 3,737, 479 reached the amplitude cutoff
in both mice only, 95 in both human sections only, 3,161 in neither pair and two in both pairs with insufficient
within-group reproducibility. The broader Fig. 5 screen flagged 474 of 873 eligible genes by allowing positional
variation supported in either group. It includes all 48 of the stricter divergent profiles. In fixed-width
raw-count bins, 48 of 48 divergent and 94 of 97 conserved profiles also met their observed shape-correlation
cutoffs. The total and split depend on how much positional variation is required in all four sections.

| Peak-position pattern among the 48 divergent profiles | Genes |
|---|---:|
| Late in mouse, early in human | 20 |
| Early in mouse, late in human | 10 |
| Between the middle and a terminal band | 12 |
| Same broad peak band with a different profile shape | 2 |
| Peak band varied between sections of one group | 4 |

The peak bands describe fitted profiles after the full-curve comparison. The late-to-early count exceeded the
early-to-late count at all three amplitude cutoffs (32 vs 22, 20 vs 10, and 6 vs 2). Among the 97 conserved
profiles, 48 peaked early in all four sections and 37 peaked late in all four. This points to stable early and late
expression patterns alongside changes in both directions.

| Functional annotation | Shape-supported | Conserved | Divergent |
|---|---:|---:|---:|
| All orthologs | 175 | 97 | 48 |
| Reactome SLC-mediated transport | 29 | 22 | 2 |
| Reactome amino-acid transport | 10 | 10 | 0 |
| KEGG lysosome | 8 | 7 | 1 |
| Reactome amino-acid, carbohydrate or fatty-acid metabolism | 38 | 20 | 11 |
| KEGG peroxisome | 8 | 3 | 4 |

These memberships overlap, and the three metabolic sets are counted once per gene in their union. The clearest
functional theme is the conservation of several transport and lysosomal profiles. The divergence fraction in the
metabolic union is similar to the full supported set, and its divergent members span both directions of change.
Peroxisome annotations offer a smaller mixed example with a higher fraction of differences. The fractions are
descriptive, with two mouse kidneys and two sections from one human donor. A shared direction of expression
change establishes a recurring positional pattern, while a common regulatory mechanism remains unresolved.

## Biological result

Early and late amino-acid transporter profiles are preserved across the sampled kidneys. Selected brush-border
peptidases and carbohydrate-metabolism enzymes show reproducible changes in where their transcripts are high
relative to other genes in the same functional group. The human samples are two cortical sections from one donor.
The measurements describe RNA profiles along the reconstructed axis.

| Program | What the individual profiles and counts show |
|---|---|
| Amino-acid transport | All 10 screened members of the archived Reactome plasma-membrane amino-acid transport set meet the conserved-shape rule. Seven stay early-weighted, including SLC6A19, SLC3A2, SLC7A7, SLC7A8 and SLC7A9. SLC3A1, SLC6A18 and SLC6A6 stay late-weighted. SLC6A19 and SLC7A9 track one another across fixed-width bins in all four sections (r = 0.91, 0.93, 0.86 and 0.88). The association persists in the lower-mitochondrial-fraction human profiles (r = 0.69 and 0.78). SLC16A10 has a shallow profile and less agreement in the low-mito check. |
| Lysosomal genes | Seven of eight shape-supported genes in the archived KEGG lysosome set retain their profile, including early-weighted CTSB, CTSL and LGMN and late-weighted HEXB. CTSH is the divergent member. |
| Brush-border peptide processing | MME and DPP4 are late-weighted in the two mice and early-weighted in both human sections. Mouse DPP4 is comparatively shallow. DPEP1 and ENPEP stay early-weighted. Across the fixed-width count bins, MME tracks DPP4 in all four sections, while the MME–DPEP1 correlation goes from −0.91 and −0.92 to +0.94 and +0.97. In the lower-mito human subsets it is +0.91 and +0.94. Of 16 screened genes in the archived KEGG protein digestion and absorption annotation, 10 qualify for the full-profile comparison. Eight transporter genes retain their shape, and MME and DPP4 diverge. ANPEP is weakly graded in the displayed interval. |
| Carbohydrate metabolism | PCK1, DCXR and MIOX are earlier-weighted in mouse and later-weighted in human. The MIOX human gradient is shallower than the mouse gradient. FBP1 and ALDOB retain a decline toward the late PT, even though positions of their interior maxima vary. PCK1–FBP1 correlation across count bins goes from +0.94/+0.90 to −0.72/−0.64 and stays negative in the lower-mito human subsets (−0.69/−0.62). |
| Peroxisome-associated genes | The shapes vary by component. CAT is middle-to-late in mouse and early in human, ACOX2 has a pronounced late rise in mouse and a shallow human profile, and HAO2 and AMACR are more late-weighted in human. IDH1 stays early-weighted. Four of eight shape-supported KEGG peroxisome genes meet the divergence rule. The signs and strengths vary with binning and mitochondrial filtering. |
| Transporter exceptions | SLC16A9 goes from late toward early and SLC16A12 goes from early toward late. Both changes survive reaggregation and the lower-mito check. Of 29 eligible genes in the archived Reactome SLC-mediated transport set, 22 retain their shape and two meet the divergence rule. This annotation includes a subset of the SLC families. |

The clearest comparison for Fig. 5 is amino-acid transport beside brush-border peptide processing. It pairs
a well-preserved early/late organization with selected enzymes whose relative profiles change. The PCK1–FBP1
comparison is a second example of a change within a named functional group.

## Data and checks

- Inspected the 3,912 one-to-one orthologs from the current Fig. 5 universe: shared detection, proximal-tubule
  versus neighboring segment filter, and convergence in the four sections. The inputs were per-section COAST
  log₂ rates and observed counts aggregated by position. No differential expression statistics or existing
  change calls determined the new groups. Local gene-set membership was joined after shape comparison.
- Compared shapes on the common 0.19–0.875 mouse-reference interval. Added 0.5 CP10k before taking logs.
  Interpolated the 20 observed equal-count bin rates to the fitting grid and smoothed with a Gaussian
  kernel of 0.025 COAST units for the shape comparison. Every qualified gene has a fitted range of at
  least 0.5 log₂ units in each section and within-group fit and observed correlations of ≥0.65 and ≥0.55,
  respectively. The 175 genes that qualify have positional variation supported in all four sections.
- Conserved shapes require all four mouse–human correlations to reach 0.70 for fits and 0.50 for observed bins.
  Divergent shapes require all four to stay at or below 0.30 and 0.45, respectively. The alternative COAST
  coordinates, constructed without the profiled gene, must agree at 0.65/0.45 for conserved or at most
  0.40/0.50 for divergent. These descriptive cutoffs yield 97 conserved, 48 divergent and 30 intermediate or
  sensitive profiles. At fitted range cutoffs of 0.35 and 0.75, the respective conserved/divergent counts are
  133/94 and 56/19. The balance of conserved and divergent curves is a feature of the well-supported subset.
- In the 97 conserved profiles, 48 peak in the early interval (≤0.35) in all four sections and 37 peak in the
  late interval (≥0.67) in all four. Among the 48 divergent profiles, 20 have all four maxima in opposite
  bands in the late-mouse-to-early-human direction, and 10 in the reverse direction. Other examples change
  between a middle and terminal peak or differ in profile shape despite peaking in the same broad region.
  These peak labels describe the full-curve groups after classification. The heatmap shows the shapes
  standardized within each section, then averaged with equal weights. Shape clusters order its rows.
- At the 0.5 log₂ range cutoff, 479 genes meet the range criterion in both mice but fail to meet it in both
  human sections, and 95 show the reverse. Different probe panels, one human donor and technical variation
  limit biological interpretation of the widespread reduction in human amplitude. The broader 474-gene
  Fig. 5 call set captures shape differences even when one group has limited amplitude. All 48 of the divergent
  shapes in this stricter audit are in that call set. Twenty genes called by the broader screen meet this
  audit's conserved-shape rule, since the two rules test different properties of the same curves.
- Re-aggregated the original counts into 20 fixed-width bins spanning the four sections' common central
  support. Rates use the sum of non-mitochondrial library sizes as the denominator. Repeated within each
  human section's bottom mitochondrial tertile. For the numeric comparisons, kept positions in the comparison
  interval with ≥20 tubules in both bottom-tertile subsets, yielding 13 shared bins between 0.28 and 0.872.
  The pair correlations above describe overlap across position bins in the sampled specimens.

## Reproduce

Run `python analysis/paper/fig5_curve_audit.py` in the analysis environment.
Outputs under `results/paper_figures/fig5/curve_audit/` include the two-panel shape heatmap, a 175-gene
atlas with section-level curves and re-aggregated counts, a 25-gene program figure, membership and per-gene
tables, and the observed pairwise correlations.
