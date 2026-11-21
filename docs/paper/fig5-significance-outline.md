# Fig. 5 Results: outline with significance per result

## Recommended Fig. 5 story (2026-11-21)

**Urate export and reuptake transporters occupy different positions relative to one another in mouse and human proximal tubules.**
The apical urate exporter NPT4 and the apical reabsorptive exchanger URAT1 have closely tracking transcript
profiles in mouse. In the human sections, NPT4 is weighted toward the early PT and URAT1 toward the late PT.
This gives the gene-level screen a physiological question to answer. Does a transporter with the same role sit
in the same part of the nephron relative to a transporter that moves the same substrate in the opposite direction?
Transport roles are based on functional studies. [Jutabha et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2966126/)
[Enomoto et al.](https://www.nature.com/articles/nature742) [Hosoyamada et al.](https://pubmed.ncbi.nlm.nih.gov/14747372/)

### Suggested Results text

#### Urate export and reuptake transporters show different axial expression relationships in mouse and human

We used the shared COAST coordinate to compare gene expression along the proximal tubule in two control mouse
kidneys and two healthy cortical sections from one human kidney. Among 3,908 assessable one-to-one orthologs,
60 had opposite early-to-late expression contrasts in mouse and human, with matching directions in each species'
two sections and agreement between fitted curves and raw counts. We examined this pattern in renal urate handling,
where NPT4 (SLC17A3) exports urate into the tubular lumen and URAT1 (SLC22A12) mediates apical urate reuptake.

In mouse, NPT4 and URAT1 transcript profiles rose together through the early proximal tubule and reached their
fitted maxima near the middle of the coordinate. In both human sections, NPT4 was highest early and URAT1 was
highest late. Raw expression from 20 position bins agreed with the fitted profiles. NPT4 and URAT1 were positively
correlated in both mice (Pearson \(r = 0.98\) in each) and negatively correlated in both human sections
(\(r = -0.66\) and \(-0.64\)). In independent microdissected mouse tubules, both transcripts were higher in
S2 and S3 than in S1 in the male and female profiles. In the Lake human kidney atlas, NPT4 was higher in
PT-S1 than PT-S3 in 11 of 12 eligible healthy cortical donors, while URAT1 was higher in PT-S3 than PT-S1
in all 12. NPT1 (SLC17A1), a second apical urate exporter, followed the direction of NPT4 in both datasets.
These profiles suggest that two opposing components of apical urate transport are more separated along the
human proximal tubule at the transcript level. [Mouse microdissection](https://doi.org/10.1681/asn.0000000000000089)
[Lake human atlas](https://doi.org/10.1038/s41586-023-05769-3)
[NPT1 transport assay](https://pubmed.ncbi.nlm.nih.gov/20566650/)

**Discussion implication.** Mouse models both transport processes, with enrichment of the two transcripts in
a similar part of the tubule. If the human transcript separation is reflected in apical protein
localization, experiments based on mouse tubules may assign the wrong relative positions to the sites of NPT4
export and URAT1 reuptake. Measuring both proteins in anatomically identified proximal tubules would test that
prediction.

### Evidence for the urate comparison

| Readout | Mouse | Human |
|---|---|---|
| NPT4 raw late-minus-early log2 rate | +1.15, +1.43 | −1.57, −1.93 |
| URAT1 raw late-minus-early log2 rate | +0.94, +1.24 | +1.58, +1.04 |
| NPT1 raw late-minus-early log2 rate | +0.83, +0.86 | −0.47, −0.75 |
| URAT1 minus NPT4 fitted peak position, mouse-reference coordinate | 0.10, 0.03 | 0.58, 0.59 |
| NPT4–URAT1 Pearson correlation, raw position bins | +0.98, +0.98 | −0.66, −0.64 |

- NPT4 was detected in 92% of PT tubule profiles in each mouse and in 74% and 66% of profiles in the human
  sections. Its direction persists after restricting late mouse tubules to the cortical proxy and restricting
  each human section to its lowest mitochondrial-fraction tertile. Both genes retain their endpoint directions
  on the archived gene-half axes. In that run, each transcript was assessed on a coordinate built without it.
  Raw endpoints use early and late cores on COAST with concordant PT labels. The independent atlas endpoint
  check uses author-labelled PT-S1 and PT-S3, without COAST or our probe panels.
- The atlas check uses paired log2(CPM + 1) differences from cortex-only normal donors with at least 50 nuclei
  in each endpoint. Median differences were −1.24 for NPT4, +1.24 for URAT1 and −0.46 for NPT1. At a 25-nucleus
  floor, NPT4 and NPT1 were earlier in 12 of 13 donors each, and URAT1 was later in all 13. NPT1 has the same
  endpoint directions internally, with one human fitted contrast below the screen's 0.5 log2 cutoff. It is
  supporting evidence from a specified transport role.
- Sakaguchi et al. already reported early weighting of NPT4 and NPT1 and late weighting of URAT1 in human
  nuclei. The contribution here is the cross-species contrast between overlapping mouse profiles and separated
  human profiles. The external sources support the endpoint directions. The peak gaps describe positions on
  the mouse-reference quantile scale, and the bin correlations summarize profiles within each section.
  [Sakaguchi et al., Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10763458/)
- The focused comparison is post hoc and supports axial patterning of the named transcripts. Quantifying
  protein abundance and net urate flux would require separate measurements. In Lake, the median human OAT1
  endpoint contrast is slightly positive, opposite to both Visium sections. Other components of the secretory
  and reabsorptive pathways have incomplete shared-panel coverage or fail expression checks.

### Figure and supporting examples

1. Keep the tissue overview and put the global screen in one compact panel. The current single-coordinate
   run has 873 shape-eligible genes, 474 descriptive shape calls and 60 stricter, raw-count-supported reversals
   among 3,908 assessable genes. These are two different denominators. No pathway passes the corrected threshold
   in either the 895-set endpoint screen or the 129-set shape screen. Move descriptive pathway means out of the
   main argument.
2. Give NPT4 and URAT1 the central panel. Put mouse and human side by side, with both genes, both section curves
   and raw bin points. Add the paired donor endpoint contrasts and sex-stratified mouse microdissection profiles
   alongside it. Use gene-wise scaling for the display and place any mouse boundaries below the mouse panel.
   The continuous profiles show overlap around the mouse mid-PT and separate expression maxima in human. A
   discrete endpoint comparison independently captures the direction of the difference.
3. Keep glucose as a shorter example of conserved ordering with a different expression transition. SGLT2
   precedes SGLT1 in both species. The crossing of their separately min–max-scaled curves is at 0.37–0.40 in
   mouse and 0.65–0.67 in the human sections on the mouse-reference coordinate. Human SLC5A1 is detected in
   about 6% of tubule profiles, and published human brush-border staining supports the ordering at segment
   resolution. Keep the crossing as an illustrative transcript transition in this cohort.
   [Vrhovac et al.](https://pubmed.ncbi.nlm.nih.gov/25304002/)
4. Keep MME and DCXR as short examples or in Extended Data. Both reverse in the current screen, with the same
   directions in mouse microdissection and paired human atlas donors. Retain the B0AT1/B0AT3 curves as descriptive
   if shown. Human SLC6A18 transcripts cannot establish an active counterpart to mouse B0AT3. Phosphate has a
   strong independent physiological literature, but PiT-2 (SLC20A2) is missing from the current shared screen,
   limiting a three-transporter comparison here. [B0AT3 function](https://pmc.ncbi.nlm.nih.gov/articles/PMC4591816/)
   [Phosphate transport study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13396815/)

### Reproduce the comparison

Use `results/paper_figures/fig5/screen` from the November 18 single-coordinate COAST run on v4fix and the
normal-donor Lake cortex/medulla pseudobulk cache from notebook 38. The mouse TPMs come from GSE212213.
The focused comparison reads the existing screen and public references without refitting or changing the
frozen screens.

```bash
python analysis/paper/fig5_urate.py \
  --lake-pseudobulk results/pt_literature_deep/stage_cache/lake_pt_pseudobulk__80a7dbd8b164f678.npz
```

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
