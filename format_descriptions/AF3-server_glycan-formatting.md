To describe the glycan chains, we are using 3-letter CCD codes (Chemical Components in the PDB) of the corresponding glycans. Please note that stereoisomers are described by different CCD codes, e.g. mannose (C6H12O6) could be described as MAN for alpha-D-mannose and BMA for beta-D-mannose.

The Server supports the following glycan residues to be attached to a protein residue
- N (Asparagine): BGC, BMA, GLC, MAN, NAG
- T (Threonine): BGC, BMA, FUC, GLC, MAN, NAG
- S (Serine): BGC, BMA, FUC, GLC, MAN, NAG

Branched glycans can be constructed in the form of a tree with either one or two downstream connections per glycan, attached to a protein residue. Up to 8 glycan residues in total are supported. Here are some examples that demonstrate how to input branching glycans:

- NAG: NAG is a single glycan residue.  
![NAG](./images/nag.png)
  
- NAG(BMA): NAG has a single child which is BMA.  
![NAG(BMA)](./images/nag-bma.png)  
  
- NAG(BMA(BGC)): NAG has 1 child which is BMA; BMA has one child which is BGC.  
![NAG(BMA(BGC))](./images/nag-bma-bgc.png)

- NAG(FUC)(NAG): NAG has 2 children which are FUC and NAG.  
![NAG(FUC)(NAG)](./images/nag-fuc-nag.png)

- NAG(NAG(MAN(MAN(MAN)))): linear glycan chain.  
![NAG(NAG(MAN(MAN(MAN))](./images/man3.png)

- NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN))))): branched ligand chain.  
![NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN))))](./images/man6_branched.png)


Glycan - glycan connections should also be chemically valid. For example, GLC(NAG)(MAN) is not a valid branched glycan because NAG and MAN can’t form glycosidic bonds to GLC.

The Server assumes that glycosidic bonds are formed between atoms at positions that have the highest frequency of occurrence in similar bonds from the PDB - this might lead to different bond positions in the modeled structure than expected. Specifying exact atoms for the glycosidic bond is not currently supported.