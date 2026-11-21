# Results: Figs. 5 and 6

Fig. 5 uses the overview produced by `analysis/paper/fig5_overview_figure.py`. Panel a shows the matched
heatmaps, b the class means, c the urate pair, d the organic-anion pair and e the raw-count overlap. The
figure legend and supporting analyses are in [the Fig. 5 outline](fig5-significance-outline.md).

The Fig. 6 draft uses the 2026-11-12 v4fix COAST results. Its figure legend appears below the Results.

### Proximal tubule gene expression follows distinct positional patterns in mouse and human

To compare how the proximal tubule (PT) is organized along its length, we aligned expression profiles from
two control mouse kidneys and two cortical sections from one human kidney on a shared COAST coordinate.
We assessed profile shape in 873 one-to-one orthologs with at least 1.4-fold variation along the PT in
both sections of at least one group. Of these, 474 passed the shape-difference screen, which compared
the observed group contrast with balanced section swaps (Fig. 5a). Among the 474, 122 reached the
amplitude threshold in both groups, 280 in mouse only and 72 in human only. The heatmaps show the same
genes in all four sections, on a scale that preserves differences in the amplitude of their profiles.

Clustering the difference profiles yielded four groups. In two groups, early or late enrichment exchanged
positions between the species. The early-to-late group contained 114 genes, including the L-xylulose
reductase DCXR. The late-to-early group contained 161 genes, including the transporters NPT4 and MRP2.
In the other two groups, the mean profiles peaked near the middle in mouse and late in human (105 genes),
or early in mouse and near the middle in human (94 genes) (Fig. 5a,b). The group means summarize the
expression patterns seen across the two mouse kidneys and the two sections from the sampled human kidney.
[DCXR biochemical role](https://pmc.ncbi.nlm.nih.gov/articles/PMC12337341/)

We asked how these differences change the relative positions of transport components. NPT4 (*SLC17A3*)
exports urate into the tubular lumen, and URAT1 (*SLC22A12*) takes it back up. Both transcripts rose together
and peaked near the middle of the coordinate in mouse. In both human sections, NPT4 peaked early and
URAT1 peaked late. Across the observed position bins, their expression profiles were positively correlated
in mouse (Pearson r = 0.98 in each kidney) and negatively correlated in the human sections (r = −0.66
and −0.64) (Fig. 5c). NPT1, a second apical urate exporter, showed the same direction of endpoint change
as NPT4. The relative URAT1 and NPT4 endpoint contrast also held with segment labels alone, gene-held-out
coordinates and restricted tissue subsets (Source Data). [NPT4 transport](https://pmc.ncbi.nlm.nih.gov/articles/PMC2966126/)
[URAT1 transport](https://www.nature.com/articles/nature742)
[NPT1 transport](https://pubmed.ncbi.nlm.nih.gov/20566650/)

A second pair showed a change in the opposite direction. The position-bin profiles of the basolateral
organic-anion uptake transporter OAT3 and the apical exporter MRP2 were negatively correlated in both
mice (r = −0.87 and −0.94) and positively correlated in both human sections (r = +0.38 and +0.74)
(Fig. 5d). OAT1 showed a similar change in its relationship with MRP2 (Source Data). In a separate
readout from the same assay, we counted PT profiles in which both transcripts ranked in the top quarter
within strata matched for library depth and mitochondrial fraction. Relative to the expectation under
independence, NPT4–URAT1 joint-high expression was lower in the two human sections than in either
mouse, and OAT3–MRP2 joint-high expression was higher (Fig. 5e). These comparisons describe the relative
placement and joint enrichment of transporter transcripts in the sampled kidneys, with both human
sections drawn from one donor. [Organic-anion transporter roles](https://www.sciencedirect.com/science/article/pii/S0006295217302794)

### COAST resolves where acute injury reshapes the proximal tubule

We applied COAST to two control mouse kidneys and two kidneys after ischaemia–reperfusion injury (IRI), placing all
four on one shared pseudostructure coordinate (Fig. 6a). We rebuilt the axis without the 3,635 injury-responsive genes;
injured tubules kept their order as well as control tubules did (Spearman ρ 0.91 and 0.92 against 0.90 and 0.95;
Methods). Because both conditions are measured on one probe panel, levels can be compared directly here.

*Havcr1*, which encodes KIM-1, rose by 10.3, 12.6 and 9.3 log2 in S1, S2 and S3 after injury (Fig. 6b,c).

### Clustering of shape differences identifies four classes of positional change after injury

We applied a screen for shape differences between the control and injured kidneys. Of 5,496 PT-expressed genes, 1,788 changed
at least 1.4-fold along the PT in both sections of at least one condition. The screen called 668 genes, and none under
either balanced swap of one control and one injured section.

Clustering the 668 difference curves gave four classes (Fig. 6f). In the first and smallest (37 genes), the mean control
shape peaked mid-PT (0.49) and the injured shape earlier (0.24); it included *Khk* (ketohexokinase), which peaked at
0.48 in controls and 0.26–0.29 after injury, and *Slc26a1* (SAT1), which kept its early peak but no longer fell in S3
after injury (Fig. 6d). In the second (180 genes), expression rose toward the end of the injured PT (peak 0.88) but not
of the control PT (peak 0.34); it included *Hmox1* (heme oxygenase 1), flat in controls and peaking at 0.84–0.88 after
injury, and the DNA synthesis gene *Tk1*. In the third (206 genes), control expression rose into S3 (peak 0.83),
whereas after injury it peaked mid-PT (0.46) and fell in S3; it included *Slc6a13* (GAT2), which peaked at 0.86–0.89 in
controls and 0.42–0.44 after injury, *Slc22a2* (OCT2; 0.80–0.82 against 0.50–0.58) and *Cyp4a14*. The fourth and
largest (245 genes) followed the same direction (control peak 0.88, injured peak 0.50); it included *Acadm*, which
peaked at 0.84–0.85 in controls and 0.51 after injury, *Havcr1* and *Cyp24a1*. *Slc6a8* (CRT1), also in this class,
kept its early peak but held its level through S2 after injury before falling in S3. In the third and fourth classes,
which together hold 451 of the 668 genes, the position of peak expression moved from S3 toward the middle of the PT
after injury.

Several genes graded within a segment in controls were nearly flat there after injury (Fig. 6e). *Cyp4a14*, *Agt* and
*Cyp24a1* each changed by more than 2-fold from the start to the end of S2 in both controls (+1.1, +1.35 and −1.5 to
−1.9 log2), but by −0.5 to −0.4, +0.1 to +0.3 and +0.1 log2 after injury. To ask whether this is general, we chose
gradients in one control mouse and measured how much of each the other mice kept. The second control kept 0.78–1.33 of
each S2 gradient, whereas each injured mouse kept 0.39–0.56, and the injured mice kept less in all 12 comparisons across
the three segments. Only 4 and 7 S2 gradients belonged to genes whose overall level barely changed after injury, too
few to separate the loss due to injury from noisier placement of injured tubules (Methods).

### Pathway-level shape differences after injury exceed the genome-wide background

Of 447 sets tested, 70 changed shape by the gene-level criterion and none under either swap; 24 of these also changed
more than random sets of the same size (q ≤ 0.05). Clustering the 70 set-level difference curves gave four groups, and
the set with the strongest competitive evidence in each is shown in Fig. 6g. E2F targets (38 genes; q = 0.01) and
Myc targets (56; q = 0.01) stayed at their lowest through the first half of the injured PT and peaked in S3 (0.81
and 0.89); in controls they were also highest at the end of the PT but higher in S1 than in mid-PT. The citric acid
cycle and respiratory electron transport (43; q = 0.01) and lysosome genes (38; q = 0.01) both peaked at the start of
the control PT (0.12 and 0.13) and mid-PT after injury (0.42 for both), and fell steeply into S3 after injury.

Taken together, these results demonstrate that COAST resolves injury along the PT. After IRI, the peak of 451 genes
moves from S3 toward the middle of the PT, respiration and lysosome programmes are lost from S3, proliferation
programmes rise in S3, and gradients within segments are about half erased. With two mice per condition these
comparisons are descriptive, and they describe one early time point (about 12 h–2 d after reperfusion).

**Fig. 6 | COAST resolves where acute injury reshapes the proximal tubule.** **a** One control and one injured (IRI)
mouse kidney section assayed by Visium HD. PT tubules are coloured by their position on the shared pseudostructure
coordinate; other tubules are grey. **b** *Havcr1* (KIM-1) expression in the same PT tubules. **c** *Havcr1*
expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each control (blue) and injured (orange)
mouse over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Dashed lines mark
the S1/S2 and S2/S3 boundaries of the control segment labels. **d** As in c for one gene from each shape class in f:
*Slc26a1* (SAT1), *Hmox1*, *Slc6a13* (GAT2) and *Slc6a8* (CRT1). **e** Three genes whose within-segment gradient is
lost after injury: *Cyp4a14*, *Agt* and *Cyp24a1*. **f** The four classes of shape difference. Each mouse's curve is
centred and scaled to unit variance; lines show the mean over the class's genes for each mouse (thin) and each
condition (thick). Background shading marks S1, S2 and S3. Gene counts above each panel; two member genes are named.
**g** As in f for one pathway from each group of pathway-level shape differences.

> Methods sentences still to write: shape screen (eligibility 1.4-fold in both sections of one condition; RMS ≥ 0.5;
> both swaps ≤ 0.5 × contrast; Ward k = 4), pathway step (≥ 10 shaped members; same rule on the mean difference
> curve; competitive p from 10,000 size-matched random sets, BH), and the 4-cluster choice (fixed, not tuned).
