# Visualisation and Statistical Analysis
## GPCR kinse knockout cells reveal the impact of individual GRKs on arrestin-binding and GPCR regulation

Julia Drube*[1] , Raphael Silvanus Haider*[1] , Edda Sofie Fabienne Matthees[1], Mona Reichel[1], 
Julian Zeiner[2], Sebastian Fritzwanker[3], Clara Ziegler[1], Saskia Barz[1], Laura Klement[1],
Jenny Filor[1], Verena Weitzel[1], Andrea Kliewer[3], Elke Miess[3], Evi Kostenis[3], 
Stefan Schulz[3] and Carsten Hoffmann[1]

[![DOI: 10.1038/s41467-022-28152-8](https://zenodo.org/badge/DOI/10.1038/s41467-022-28152-8.svg)](https://doi.org/10.1038/s41467-022-28152-8)

[1]: Institut für Molekulare Zellbiologie, CMB – Center for Molecular Biomedicine, Universitätsklinikum Jena, Friedrich-Schiller-Universität Jena, Hans-Knöll Straße 2, D-07745 Jena, Germany

[2]: Molecular, Cellular and Pharmacobiology Section, Institute for Pharmaceutical Biology, University of Bonn, Nussallee 6, 53115 Bonn, Germany 

[3]: Institut für Pharmakologie und Toxikologie, Universitätsklinikum Jena, Friedrich-Schiller-Universität Jena, Drackendorfer Straße 1, D-07747 Jena, Germany

[*] contributed equally

---


This code was written and used for statistical analysis and visualisation of data included in 
[Drube et al. 2021](https://doi.org/10.1038/s41467-022-28152-8). Cite this code [![DOI: 10.5281/zenodo.5764248](https://zenodo.org/badge/DOI/10.5281/zenodo.5764248.svg)]

---

### Statistical Analysis of all twelve tested GPCRs
[**ST1_vehicle_vs_stim.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/ST1_vehicle_vs_stim.R)
: All concentration response curves presented, were statistically analysed to determine functional recruitment. 
Results are listed in Supplementary Table 1.

Beta-arrestin recruitment data from all twelve tested GPCRs (Supplementary Figure 5 and 6) was preprocessed using
[**curve_formatting.py**](https://github.com/mo-yoda/Drube_2021/blob/main/Preprocessing/curve_formatting.py). 
Concentration response curves (Supplementary Figure 5) were analysed utilizing
[**F3i_ST3_curve_analysis.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F3i_ST3_curve_analysis.R). 
Results (Supplementary Table 3) were ultimately plotted as heatmap (Figure 3i, 
[**F3i_heatmap.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F3i_heatmap.R)). 
To identify increased baselines 
[**SF6,8_baseline_analysis.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/SF6,8_baseline_analysis.R)
was employed (Supplementary Figure 6, 8).

### Phosphorylation Pattern Analysis
Potential phosphorylation sites (P), clusters (PPP, PXPP 
[ref](https://pubmed.ncbi.nlm.nih.gov/10542263/)
) and patterns (PXPXXP, PXXPXXP
[ref](https://pubmed.ncbi.nlm.nih.gov/28753425/)
) were identified in 
the intracellular loop 3 (IL3) and the C-terminus (C-term) of the investigated GPCRs using 
[**moving_frame.py**](https://github.com/mo-yoda/Drube_2021/blob/main/Phosphorylation_pattern/moving_frame.py).
Data of following format was used as input in which each GPCR - beta-arrestin pair is listed with their GRK preference, the 
class of the GPCR (according to 
[Oakley et al. 1999](https://pubmed.ncbi.nlm.nih.gov/10542263/)
) and the IL3 or C-term amino acid sequence. Aminio acid sequences were obtained from [GPCRdb](https://gpcrdb.org/).

| GPCR | barr | GRK | arr_class | IL3_or_Cterm | seq |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| b2AR | barr1 | GRK2356 | A | IL3 | GRFHVQNLSQVEQDGRTGHGLRRS |

The number of identified phosphorylation sites, clusters and patterns are listed as part of Supplementary Table 2.
The count and relative position of these phosphorylation sites were visualized grouped by GRK preference or class 
[**F7_pattern_analysis.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Phosphorylation_pattern/F7_pattern_analysis.R). 
Generated plots are presented in Figure 7. A possible association between the relative position of PXPP clusters and 
GRK preference or class was investigated in 
[**F7_Fisher_details.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F7_Fisher_details.R).

### Statistical Analysis of Miscellaneous Datasets
Data processed with the following scripts was imported in following format in which each signal recorded 
by the respective method ("Output") can be linked to a certain factor which was variied between the samples 
("Condition"). A within factor variable ("State") was included for data measured before and after GPCR stimulation.

| Condition | State | Output |
| ----------- | ----------- | ----------- |
| dQ+EV | baseline | 1.03923 |
| dQ+EV | stimulated | 0.81442 |
| dQ+GRK2 | baseline | 1.93874|
| dQ+GRK2 | stimulated | 1.51442 |


[**F1d,f_SF3c_GRK_expression.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F1d%2Cf_SF3c_GRK_expression.R):
Statistical analysis of GRK expression data presented in Figure 1d, f and Supplementary Figure 3c was statistically analysed.

[**F2b,c,d_EC50.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F2b,c,d_EC50.R):
Comparison of several EC50 as presented in Figure 2b-d.

[**F3g,h_SF7c,e_confocal.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F3g%2Ch_SF7c%2Ce_confocal.R):
Comparison of colocalisation quantified from confocal microscopy before and after stimulation (Figure 3g ,h and 
Supplementary Figure 7c, e).

[**F5_SF11_ST4_AT1R.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/F5_SF11_ST4_AT1R.R):
 Statistical analysis of beta-arrestin recruitment to AT1R under various conditions as presented in Figure 5 and 
Supplementary Figure 11. Results are listed in Supplementary Table 4. Data was preprocessed before statistical 
analysis using
[**compare_format.py**](https://github.com/mo-yoda/Drube_2021/blob/main/Preprocessing/compare_format.py).

[**SF13_Losartan_Tolvaptan_analysis.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/SF13_Losartan_Tolvaptan_analysis.R):
Statistical analysis of beta-arrestin recruitment under application of inverse agonists as displayed in Supplementary Figure 13.

[**SF1,2_stats.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/SF1%2C2_stats.R):
Statistical analysis of datasets presented in Supplementary Figure 1 and 2.

[**SF9,10d,f_V2R_endo_KD_analysis.R**](https://github.com/mo-yoda/Drube_2021/blob/main/Statistical_Analysis/SF9%2C10d%2Cf_V2R_endo_KD_analysis.R):
Statstical analysis of beta-arrestin recruitment to V2R (Supplementary Figure 9) and AT1R 
(Supplementary Figure 10 d and f) in presence of catalitically inactive GRKs or endogenous expression 
of one specific GRK isoform.
