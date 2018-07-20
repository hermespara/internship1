# Plan Stage

## Introduction

* Introduction de la transcription
* Epissage alternatif
* Intron retention (spliceosome, U2, U12)
* NMD
* Analyse par RNAseq
* Données GTEx
* L'étude sur les données : Comparaison de la rétention des introns U2 et U12

## Matériels et Méthodes

### Matériels

* Données GTEx
* Travail dèjà fait par l'équipe

### Méthodes

* Séléction des individus
* Pipeline bio-informatique (Snakemake)
    * fastqc + multiqc
    * sra toolkit
    * STAR
    * IRFinder
    * HTSeq count
    * Supprimer les fichiers brut
    * anonymization
    * junctionCover2IRF
    * indexation samtools
* Analyse statistique:
    * DESeq2 : expression des gènes
    * kissDE

## Résultats

* Automatisation du pipeline d’analyse
* Analyse de 40 échantillons (20 individus 2 tissus)
* Résultats DESeq2  
    * Graph
* Résultats kissDE

## Discussion


## Conclusion
