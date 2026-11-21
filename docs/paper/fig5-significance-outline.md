# Fig. 5 Results: outline with significance per result

## Fig. 5 overview

**Proximal tubule expression changes follow four positional patterns.** The two heatmaps show all 474
orthologs called by the shape-difference screen in the same gene order. They average the two mouse kidneys
and the two cortical sections from one human kidney, respectively. The urate and organic-anion pairs show
how changes in individual profiles alter the relative positions of transport components. All measurements are
from the study's Visium HD sections.

### Manuscript Results

The [mouse–human Results](results_fig5_fig6.md) first present the screen, matched heatmaps and four average
profiles. The urate and organic-anion examples then compare two pairs with different patterns of overlap.
The raw-count test uses their original tubule profiles with matching for library depth and mitochondrial fraction.

### Figure legend

**Fig. 5 | Proximal tubule expression patterns in mouse and human.** **a,** Fitted profiles of all 474
shape-difference calls, from 873 eligible one-to-one orthologs. A gene was eligible if its fitted profile varied
by at least 0.5 log₂ units (about 1.4-fold) in both sections of at least one group. Of the called genes, 122 meet
that amplitude criterion in both groups, 280 only in mouse and 72 only in human. Profiles show
log₂(CP10k + 0.5), centered within each section. For each gene, all four curves are divided by the largest of
its four within-section SDs. The two heatmaps average the scaled curves from the two mouse kidneys and the
two human cortex sections, respectively, with equal section weights. Rows are ordered by shape cluster and
then by the magnitude of the screened difference.
Each column uses the common fitting interval at 0.09–0.945 on the COAST mouse-reference quantile scale.
**b,** Section-level mean profiles for the four groups from a. Group names describe the mean peak positions
from mouse to human. **c,d,** Per-section curves and observed position bins for the NPT4–URAT1 and OAT3–MRP2
pairs. Each log₂ rate is expressed relative to the gene's fitted peak for that section and divided by the
gene's largest fitted range across the four sections. **e,** Joint-high profiles relative to the expected count
under independence. High means the top quarter of expression within library-depth and mitochondrial-fraction
strata. The four markers represent two mouse kidneys and two cortical sections of one human kidney.

### Evidence and interpretation

- The maps and joins use the same COAST v4fix run as the curves. At each threshold, rank nonzero counts
  per non-mitochondrial library separately for each gene in each of ten near-equal depth strata of a section.
  Select up to the specified fraction of each stratum, with stable tie order. At 25%, observed joint-high
  counts are 244, 258, 171 and 192. The respective within-stratum independence expectations are 170.5,
  173.8, 198.8 and 210.5. In both human cortex sections, many profiles carry each transcript and both can
  be high in the same profile. The contrast is in the frequency of strong joint expression. The independence
  comparison is descriptive, and the four sections are the units displayed. With mitochondrial-fraction
  tertiles nested within depth strata, the human joint-high rates are close to the independence expectation.
- The supporting urate analysis includes all six assessable genes from the eight transporters previously analysed by apical or
  basolateral role. SLC22A11 and SLC2A9 are absent from this assessable shared screen. Among the three
  basolateral uptake genes shown, OAT1 and OAT3 are early-weighted and OAT2 late-weighted in both species.
  The opposing NPT4 and URAT1 profiles motivate a component-level claim. No Fig. 5 pathway set passes the
  corrected threshold. [Transport role annotations](https://pmc.ncbi.nlm.nih.gov/articles/PMC10763458/)
- The rank-free checks still start from the study's assigned PT profiles. COAST position is on a
  mouse-reference quantile scale. The glomerulus-distance subset is a cortical proxy, and the tissue
  centroids do not trace a complete tubule. The similar mouse profiles and divergent human profiles
  motivate a spatial protein-localization question. They do not measure transport activity.

### Comparison with alternative biological stories

We used the study's four Visium HD sections to compare 10 candidate stories and 19 transcript pairs. For each
available gene, the comparison records detection, raw endpoint contrasts using three definitions and the
archived gene-held-out result. For each pair, it compares log-expression correlations across position bins
with the overlap of high profiles in the tissue, using three thresholds. The additional tissue analysis
matches mitochondrial fraction as well as library depth. This is an exploratory comparison of the named
candidates, with no inference treating tubule profiles as independent animals or human donors.

| Candidate | What the four sections show | Paper role |
|---|---|---|
| **Urate export and reuptake** | NPT4 and URAT1 track together in mouse (bin r = +0.98/+0.98) and run in opposite directions in human (−0.66/−0.64). Both transcripts are detected in ≥66% of profiles in every section. The difference withstands alternative endpoints and the same-tubule check with depth and mitochondrial matching. NPT1 supports the exporter pattern. | **Fig. 5 centerpiece.** The two genes act on the same substrate in opposing directions, with evidence from the coordinate and tubule profiles. |
| **Organic-anion uptake and efflux** | OAT3–MRP2 bin correlations switch from −0.87/−0.94 in mouse to +0.38/+0.74 in the human sections. OAT1–MRP2 shows the same change (−0.78/−0.74 to +0.54/+0.62). The top-quarter OAT3–MRP2 joint-high ratio with depth and mitochondrial matching is 0.85/0.76 in mouse and 0.96/1.15 in human. MRP2's weak human endpoint changes include one sign exception after QC restriction. | **Strongest transport counterpart.** Uptake and efflux transcripts track more closely in human, complementing the urate pattern. This is a statement about the named transcripts and their spatial association. [Transport roles](https://www.sciencedirect.com/science/article/pii/S0006295217302794) |
| **Apical peptide processing** | MME and DPP4 both have positive mouse and negative human raw endpoint contrasts across the four definitions, with a smaller DPP4 contrast in mouse. ANPEP stays broadly expressed. MME–ANPEP joint enrichment is lower in both human sections after depth and mitochondrial matching (0.92/0.95 vs 1.13/1.15 in mouse). | **Best distinct supporting example.** Selected apical peptidases change geography. The result supports a transcript-localization claim for named enzymes. [Brush-border enzyme roles](https://pubmed.ncbi.nlm.nih.gov/2895788/) |
| **Glucuronate–xylulose conversion** | DCXR, SORD and XYLB encode three consecutive reactions. All three have early mouse and late human raw contrasts with the standard endpoints, segment labels alone and held-out axes. The SORD human contrast is small and changes sign in one low-mitochondrial-fraction section. DCXR clears the reversal rule on both axes, and the pair overlap does not show wholesale reorganization. | **Mechanistically coherent metabolic lead.** Keep DCXR as the supported reversal, with XYLB and SORD as graded supporting observations. [Biochemical sequence](https://pmc.ncbi.nlm.nih.gov/articles/PMC12337341/) |
| **Kynurenine–NAD branch** | ACMSD is early in mouse and late in human, whereas HAAO is late-weighted across all four sections. Their bin correlations change from −0.53/−0.55 to +0.76/+0.75. ACMSD's raw direction persists on the held-out axis, where it misses the strict support filter. QPRT has a small human contrast. | **Metabolic hypothesis.** There is a change in the geography of the branch-point enzyme. The data give a weaker read of the downstream route. [ACMSD function](https://pmc.ncbi.nlm.nih.gov/articles/PMC6448761/) |
| **Glucose reabsorption** | The SGLT2-to-SGLT1 profile crossing is later in the human sections (0.65–0.67 versus 0.37–0.40). Human SLC5A1 is detected in roughly 6% of profiles. The well-covered entry and exit transcripts SGLT2 and GLUT2 retain strongly positive position-bin correlations in both species. | **Short illustrative example.** Conserved ordering with a later fitted crossing. |
| **Phosphate, neutral amino acids, lactate and prostaglandins** | Phosphate transporter gradients keep their signs, and PiT-2 lacks shared-screen support. The late amino acid and lactate partners have sparse human counts. HPGD reverses robustly, while the uptake and downstream prostaglandin components have low mouse detection or inconsistent QC contrasts. | Individual examples or reserves. |

Keep **urate** as the lead. If the Results can carry a second focused transport observation, add the
OAT3–MRP2 comparison in the opposite direction. If a distinct physiological process fits the section better,
use the MME–DPP4 result. The shared plots and full four-section readouts are saved in the story comparison
output. The human comparisons remain scoped to one donor.

### How the precedents structure biological discovery

| Paper | How the biological section develops | Application to Fig. 5 |
|---|---|---|
| **Yakubovsky et al., Fig. 3** | Starts with a quantified overview of shared and discordant zonation. Names examples from the comparison, then groups them into a specific biological process such as the urea cycle. Additional species test the scope of the observation. | Use the ortholog screen to locate the question. Move quickly to the urate comparison, with the conclusion scoped to the two mouse kidneys and the sampled human kidney. [Results and Fig. 3](https://www.nature.com/articles/s41586-026-10377-y) |
| **Harnik et al., Figs. 3–4** | Evaluates whether anatomical sampling accounts for the cross-species pattern, then compares concordant and discordant genes and examines the tissue context. | Report how the key contrast behaves with segment labels, held-out coordinates and restricted subsets. Use raw tubule profiles to check joint expression. [Fig. 3](https://www.nature.com/articles/s41586-024-07793-3/figures/3) [Fig. 4](https://www.nature.com/articles/s41586-024-07793-3/figures/4) |
| **Moor et al., Figs. 5–6** | Moves from nutrient classes to individual profiles and spatial measurements, then asks whether the profiles of related genes agree. The purine example uses biochemical roles to explain the choice of genes. | Describe NPT4 and URAT1 together and quantify how often both transcripts are high in the same tubule profile. The raw counts and fitted profiles come from the same Visium assay. [Figures and Results](https://shalevlab.weizmann.ac.il/wp-content/uploads/2018/09/Moor_2018.pdf) |
| **Halpern et al., Figs. 3–4** | Sorts thousands of zonated profiles into a single heatmap, then uses the sequence of reactions in bile-acid synthesis to motivate a test of where each enzyme is expressed relative to the next. | Organize the urate example by opposing apical transport roles and the organic-anion example by basolateral uptake followed by apical export. The evidence concerns the relative positions of transcripts. [Results and Fig. 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC5321580/#F4) |

The manuscript uses these papers as structural precedents. Every comparative expression result in the
Fig. 5 section comes from the study's Visium HD sections.

### Reproduce the comparison

Use the November 18 single-coordinate COAST screen on v4fix and its archived gene-half run.
The renderer uses the corresponding COAST QC, pass-1 labels and pseudospace files.

```bash
python analysis/paper/fig5_urate.py
python analysis/paper/fig5_urate_figure.py
python analysis/paper/fig5_story_comparison.py
python analysis/paper/fig5_overview_figure.py
```

The Fig. 5 overview and its heatmap row order are written to
`results/paper_figures/fig5/overview/`. The focused urate diagnostics and sensitivity tables are written to
`results/paper_figures/fig5/urate/`.
The alternative-story plot and its two source tables are written to
`results/paper_figures/fig5/story_comparison/`.

## Earlier outline (2026-11-11)

Historical estimates below use earlier coordinate and screen versions. Each result: what we measured, why it matters clinically, what follows.
Sources were read in full text; quotes are in `citation_ledger.md`. "Speculate" marks anything beyond our data.

---

### R1. SGLT2 occupies a longer stretch of the human PT

- **Our data.** *Slc5a2* above half-maximum through 0.38–0.39 of the mouse PT, 0.66–0.71 of the human PT.
  *Slc2a2* (GLUT2) ends at the same place in both (0.55–0.58). Both swaps < 1/4 of the contrast.
- **Significance.** SGLT2 is the gliflozin target; gliflozins are approved for type 2 diabetes, heart failure and
  chronic kidney disease [PMC11462946]. SGLT2 carries ~97% of filtered glucose; under inhibition, SGLT1 downstream
  recovers part of it, leaving 50–60% excreted in mice [Rieg 2014]. Rat modelling: inhibition moves Na+ reabsorption
  downstream, S3 oxygen consumption +102% [Layton 2015].
- **Implication.** Both depend on how much PT lies past SGLT2. Longer human SGLT2 domain -> shorter downstream
  stretch -> mouse estimates of SGLT1 compensation and late-PT workload may not carry over; also relevant to dual
  SGLT1/2 inhibition (sotagliflozin).
- **Status.** Written and committed (e7868f6).

---

### R2. Secretory transporters sit at different positions

- **Our data.** OAT1 (*Slc22a6*) and OAT3 (*Slc22a8*) peak at the human PT start (0.06–0.16), early-to-mid in mouse
  (0.24–0.47). MRP2 (*Abcc2*) starts at the human PT start, past the mouse midpoint. OCT2 (*Slc22a2*) starts later in
  human (0.51–0.57 vs 0.32–0.36). URAT1 (*Slc22a12*) peaks at 0.84 in human, 0.54 in mouse. All swaps < 1/4.
- **Significance (sources checked, text still to write).**
  - OAT1/OAT3: uptake of anionic drug classes; regulatory drug-drug interaction transporters; coupled to apical
    MRP2/MRP4 efflux [Nigam 2015, PMC4281586]. Replaces the placeholder PMC3675737.
  - OAT1 and the antiviral tenofovir [Kohler 2011, PMC3103636].
  - OCT2 and the chemotherapy agent cisplatin; human *SLC22A2* variant [Filipski 2009, PMC2746866].
  - URAT1: urate reabsorption, target of uricosurics (benzbromarone, lesinurad) [PMC11811179].
- **Implication (to write).** Uptake and efflux co-localize at the human PT start but are offset in mouse; the human
  OCT2 window is shorter and later; the uricosuric site of action sits later in human. Mouse-based maps of where
  these drugs act should be checked by position before transfer.
- **Status.** Significance paragraph NOT written (three classifier stops). Sources verified.

---

### R3. Switch points move later along the human PT

- **Our data (2-fold amplitude rule, 2026-11-11).** 12 pairs with a switch in all four sections; 7 shifted, all later
  in human, 0 under swaps. SGLT2 -> SGLT1 +0.30, SGLT2 -> SMIT2 +0.26, SLC22A23 -> URAT1 +0.25, y+LAT1 -> AGT1 +0.24,
  SGLT2 -> SGLT5 +0.21, LAT2 -> AGT1 +0.19, B0AT1 -> B0AT3 +0.18. 6/6 checkable pairs reproduce in Lake 2023. OCT2 pair
  dropped (human partner changes < 2-fold). Human B0AT1 stays the more abundant transcript through S3 (ours and Lake).
- **Significance.** A later switch = the early partner carries more of the length in human. Matters when the early
  partner is a drug target. B0AT1 (*Slc6a19*) is a phenylketonuria target: *Slc6a19*-null mice excrete >100-fold more
  phenylalanine and the mouse model improves [Belanger 2018, PMC6124451].
- **Implication (to write).** An inhibitor of an early partner may cover a larger share of reabsorption in human than
  the mouse predicts (B0AT1); the later SGLT2 -> SGLT1 switch agrees with R1.
- **Status.** To write; source verified.

---

### R4. Segment-level clustering misses most of the shift

- **Our data (after the amplitude rule).** Cluster-mean lines recover most of the three amino-acid shifts but only
  +0.02 of +0.21 (SGLT5) and +0.07 of +0.25 (URAT1), and no switch for SMIT2. Post hoc. Weaker than before: half missed.
- **Significance.** Every human switch lies inside S2 (0.53-0.60), where three means give one value.
- **Implication (to write).** Drug-target maps built from segment-resolved atlases would place these handovers at
  similar positions in both species and miss the differences in R1-R3. One caveat already in the text: three cluster
  means are a simple baseline; finer clustering may recover more.
- **Status.** Result and caveat written; the "what this costs a reader" sentence still to add.

---

### R5. Mouse-expanded gene families carry extra zonation

- **Our data.** Mouse: 9.9% (1:many) and 14.2% (many:many) strongly zonated vs 3.8% (1:1); MH OR 2.24, P = 8e-7.
  No excess where the human counterpart expanded (1.7%). Human: same ordering, OR 1.57, P = 0.17.
- **Significance.** Families named in the text (*Cyp4a*, *Cyp2d*, *Ugt2b*, *Slco1a*) are drug-metabolizing or
  drug-transporting families [TO VERIFY - no source read yet].
- **Implication (to write).** Positional drug handling measured in the mouse PT may partly sit in gene copies with no
  single human counterpart; a one-to-one comparison drops them silently.
- **Status (2026-11-11).** Moved to Extended Data; one sentence stays in Results. Sex check (post hoc): 10 of 41
  strongly zonated mouse non-1:1 genes are male-biased (Xiong 2023) against 63/380 of 1:1; dropping those 10 still gives
  OR 2.04, P = 7e-4 (Fisher, denominators kept). The drug-family significance claim is dropped unless a source is read.

---

### Close

- COAST compares positional expression across species on one continuous coordinate.
- In the PT: glucose, secretory and urate transporters occupy different positions; switch points shift; segment-level
  clustering misses most of it.
- Limits already in the text: one donor, two cortical sections; transcript position != protein position; staining in
  more donors would confirm.

---

### Open

- R2 significance paragraph blocked three times by the response classifier; user to draft, Claude to check against
  sources and place.
- R5 needs a source for the drug-family claim.
- Replace placeholder PMC3675737 with Nigam 2015 (PMC4281586).
- Two different "Extended Data Fig. X" placeholders (steepness panel; *Abcc2* and coverage).
