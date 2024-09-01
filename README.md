[![CC BY 4.0][cc-by-shield]][cc-by]

# PDBpepDS

For this comprehensive database of peptides, the PDB files of all the proteins were downloaded from [RCSB](https://www.rcsb.org/) on November 1st, 2023. 

## File information

This folder contains the following datasets:

* `output_full.csv`, which contains all the peptides, including peptides with a length of less than 15 AA as well as those with a length of precisely 15 AA.

* `peptides_length_15aa.csv`, this file contains the peptides strictly with a length of 15 AA only.

* `peptides_less_than_15aa.csv`, this file contains the PDB Code, chain ID and start/end points, and sequence for peptides that have a length less than 15 AA only.

* `peptides_with_unknown_residues.csv`, this file contains the PDB Code, chain ID and start/end points for peptides that have unknown residues in their PDB file (these are not included in either of the above files.)

## Descriptors

* `MOLECULAR_WEIGHT`:  This calculates the molecular weight of a protein sequence. It is calculated as the sum of the mass of each amino acid using expasy scale. 

* `HYDROPHOBIC_MOMENT`: Compute the maximal hydrophobic moment of a protein sequence. This computes the hydrophobic moment based on Eisenberg et al (1984). Hydrophobic moment is a quantitative measure of the amphiphilicity perpendicular to the axis of any periodic peptide structure, such as the α-helix or β-sheet. The calculation is done with angle=100, and window=11.

* `ALIPHATIC_INDEX`: The aliphatic index of a protein was proposed in Ikai (1980). It is defined as the relative volume occupied by aliphatic side chains (Alanine, Valine, Isoleucine, and Leucine). It may be regarded as a positive factor for the increase of thermostability of globular proteins.

* `CHARGE_PH7`: This function computes the theoretical net charge of a peptide sequence, based on the Henderson-Hasselbach equation described by Dexter S. Moore (1985). The net charge is computed at pH 7 using Lehninger pKa scale.

* `HYDROPHOBICITY_INDEX`: This function calculates the hydrophobicity index of an amino acid sequence by averaging the hydrophobicity values of each residue using the KyteDoolittle scale.

* `INSTABILITY_INDEX`: This function calculates the instability index proposed by Guruprasad et al (1990). This index predicts the stability of a protein based on its dipeptide composition. A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.

* `ISOELECTRIC_POINT`: Compute the isoelectric point of a protein sequence using EMBOSS pKscale. The isoelectric point (pI), is the pH at which a particular molecule or surface carries no net electrical charge.

* `MASS_OVER_CHARGE`: This calculates the (monoisotopic) mass over charge ratio (m/z) for peptides, as measured in mass spectrometry.

* `POSITIVELY_CHARGED_RESIDUES`: The number of positively charged amino acids in the given sequence. Arginine, Histidine, and Lysine.

* `NEGATIVELY_CHARGED_RESIDUES`: The number of positively charged amino acids in the given sequence. Glutamate and Aspartate are negatively charged amino acids

* `HYDROPHOBIC_RESIDUES`: The number of hydrophobic residues in the given sequence. Alanine, Isoleucine, Leucine, Methionine, Valine, Phenylalanine, Tryptophan, and Tyrosine.

* `HYDROPHILIC_RESIDUES`: The number of hydrophilic residues in the given sequence. Aspartic Acid, Glutamic Acid, Asparagine, Glutamine, Arginine, and Lysine.

* `SOLVENT_ACCESSIBLE_SURFACE_AREAS`: Calculation of solvent accessible surface areas. Uses the “rolling ball” algorithm developed by Shrake & Rupley algorithm, which uses a sphere (of equal radius to a solvent molecule) to probe the surface of the molecule. Radius of the probe is 1.40, roughly the radius of a water molecule. Resolution of the surface of each atom is 100. A higher number of points results in more precise measurements, but slows down the calculation. 

* `AROMATICITY`: The aromaticity value of a protein according to Lobry, 1994. It is simply the relative frequency of Phe+Trp+Tyr.

* `GRAVY`: Grand Average of Hydropathy value for protein sequences, calculated according to Kyte and Doolittle.

* `BLOSUM1-10`: BLOSUM (BLOcks SUbstitution Matrix) matrix is a substitution matrix used for sequence alignment of proteins. BLOSUM indices were derived of physicochemical properties that have been subjected to a VARIMAX analysis and an alignment matrix of the 20 natural AAs using the BLOSUM62 matrix.

* `PP1-3`: The Cruciani properties are a collection of scaled principal component scores that summarize a broad set of descriptors calculated based on the interaction of each amino acid residue with several chemical groups (or “probes”), such as charged ions, methyl, hydroxyl groups, and so forth.
  - PP1: Polarity
  - PP2: Hydrophobicity
  - PP3: H-bonding

* `F1-6`: The FASGAI vectors (Factor Analysis Scales of Generalized Amino Acid Information) are a set of amino acid descriptors, that reflect hydrophobicity, alpha and turn propensities, bulky properties, compositional characteristics, local flexibility, and electronic properties, that can be utilized to represent the sequence structural features of peptides or protein motifs.
  - F1: Hydrophobicity index,
  - F2: Alpha and turn propensities,
  - F3: Bulky properties,
  - F4: Compositional characteristic index,
  - F5: Local flexibility,
  - F6: Electronic properties

* `KF1-10`: The Kidera Factors were originally derived by applying multivariate analysis to 188 physical properties of the 20 amino acids and using dimension reduction techniques.
  - KF1: Helix/bend preference,
  - KF2: Side-chain size,
  - KF3: Extended structure preference,
  - KF4: Hydrophobicity,
  - KF5: Double-bend preference,
  - KF6: Partial specific volume,
  - KF7: Flat extended preference,
  - KF8: Occurrence in alpha region,
  - KF9: pK-C,
  - KF10: Surrounding hydrophobicity

* `MSWHIM1-3`: MS-WHIM scores were derived from 36 electrostatic potential properties derived from the three-dimensional structure of the 20 natural amino acids.

* `E1-5`: The Physical-Chemical Properties descriptors of a peptide. The PCP descriptors were constructed by performing multidimensional scaling of 237 physical-chemical properties.

* `PD1-2`: The Physical Descriptors of a peptide. The PP descriptors were constructed by improving on existing PCA-derived descriptors (Z-scales, MS-WHIM and T-scales) after correcting for the hydrophilicity of Methionine, Asparagine and Tryptophan based on Feng et al. 
  - PD1 - A descriptor related to residue volume. 
  - PD2 - A descriptor related to hydrophilicity.

* PROTFP1-8: The ProtFP set was constructed from a large initial selection of indices obtained from the AAindex database for all 20 naturally occurring amino acids.

* `SV1-4`: These vectors were obtained in Sneath (1996) by running PCA on the ϕ coefficient to explain the dissimilarity between the 20 natural amino acids based on binary state encoding of 134 physical and chemical properties (such as presence/absence of a —CH₃ group, step-wise optical rotation, etc.).
  - SV1 - A descriptor representing mainly aliphatic properties of each residue (AAindex:SNEP660101).
  - SV2 - A descriptor putatively modeling the number of reactive groups (AAindex:SNEP660102).
  - SV4 - A descriptor representing the aromatic properties of each residue (AAindex:SNEP660103).
  - SV4 - A descriptor with uncertain interpretation (AAindex:SNEP660104).

* `ST1-8`: The ST-scales were proposed in Yang et al (2010), taking 827 properties into account which are mainly constitutional, topological, geometrical, hydrophobic, electronic, and steric properties of a total set of 167 amino acids.

* `SVGER1-11`: SVGER descriptors were constructed by Principal Component Analysis of 74 geometrical descriptors (SVGER1 to SVGER6), 44 eigenvalue descriptors (SVGER7, SVGER8 and SVGER9), and 41 Randić descriptors (SVGER10 and SVGER11) computed for the 20 proteinogenic amino acids.

* `T1-5`: The T-scales are based on 67 common topological descriptors of 135 amino acids. These topological descriptors are based on the connectivity table of amino acids alone, and to not explicitly consider 3D properties of each structure.

* `VHSE1-8`: The VHSE-scales (principal components score Vectors of Hydrophobic, Steric, and Electronic properties), are derived from principal components analysis (PCA) on independent families of 18 hydrophobic properties, 17 steric properties, and 15 electronic properties, respectively, which are included in total 50 physicochemical variables of 20 coded amino acids.
  - VHSE1 - A descriptor representing hydrophobic properties. 
  - VHSE2 - Another descriptor representing hydrophobic properties. 
  - VHSE3 - A descriptor representing steric properties. 
  - VHSE4 - Another descriptor representing steric properties. 
  - VHSE5 - A descriptor representing electronic properties. 
  - VHSE6 - A second descriptor representing electronic properties. 
  - VHSE7 - A third descriptor representing electronic properties. 
  - VHSE8 - A fourth descriptor representing electronic properties.

* `Z1-5`: The Z-scales were proposed in Sandberg et al (1998) based on physicochemical properties of proteogenic and non-proteogenic amino acids, including NMR data and thin-layer chromatography (TLC) data.
  - Z1 - A descriptor quantifying lipophilicity.
  - Z2 - A descriptor modeling steric properties like steric bulk and polarizability.
  - Z3 - A descriptor quantifying electronic properties like polarity and charge.
  - Z4 - A descriptor relating to electronegativity, heat of formation, electrophilicity and hardness.
  - Z5 - Another descriptor relating to electronegativity, heat of formation, electrophilicity and hardness.

## Citing the dataset

If you use this dataset for academic work, please cite it using
```
PDBpepDS, Luckman Qasim and Laleh Alisaraie, 2024
```

This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
