# Intro
U ovom dokumentu navodimo što smo diskutirali na sastancima, zapisujemo ideje i TODO.
Ideja je da Katarina, Filip i Dalibor ovo popunjavaju nakon svakog sastanka. 
Ne treba biti pre detaljno, ali će nam pomoći da ostanemo u fokusu i da ne zaboravimo na dobre ideje :)

## sastanak 22.02.2024.
_Katarina_

Napraviti poseban branch za sve prezentacije, postere i predavanja. U naslov bi trebao biti timestamp i za koga je bila prezentacija. Timestamp u obliku godina-mjesec (npr. 2024_02_ZEL_seminar). 

Consider a visit to Davor Lessel lab to Salzburg. 

**Dalibor**

Preparing powerpoint for JC. 

**Katarina**
Add genome shuffling data for metabolic syndromes and cardiac.

**Filip**

Prepared tensor board tool for visualization dimensionality reduction methods. 

## sastanak 26.01.2024.
_Katarina_


**Anja** 
* Submitted the annual project report (each group member could read the report to get a feeling of how te project stands)
*  One task that included adding on-the-fly option for targPred dataset was abandoned, and insted Filip will focus on redefining enhancers with genomic LLMs
*  Another task that should be done is integrating other data sources for defining enhancer-promoter interactions i.e. Roadmap, Hi-C (Radina did this), GTEx, STAR-seq. Also there are resources that already incorported this such as OpenTargets and GeneHancer. 


**Filip** 
* Compared epiEnh with non-epiEnh SNPs and their respective TFBS. The results indicate show developmental TFs in the epiEnh group.
* Next task is to make an enrichment analysis on GO ontology of those TFs, the background used should be constrained only on TFs.
* Will update the About page on targPred platform. 

**Katarina**
* Needs to test if trait associated SNPs are more likely to appear in GRBs with genome shuffling.
* Additonally, consider revisiting ClinVar coding variants and their ontology to expand the ontology connection. 

**Dalibor**
* Incorporated eniEnh data set into SOM results.
* Consider preparing a swarm (box, violin) plot of fractions of the epiEnh dataset per superclass. 

**Miha** 
* Is preparing a review paper about NSAI for a conference.
* Evaluate the about page from targPred. 

**General**
* Anja, Filip and Katarina should map out the thesis topics
* Plan out tasks for writing the targPred paper
* Book conferences, educations
* Half-project meeting
* Look up books

## sastanak 05.01.2024.
_Katarina_

**Phenotpes with LRR**

Using LD blocks for SNPs from the whole GWAS catalog (2021), that overlap enhancers (fantom5 enhancers) were separated into two categories: inGRB and outGRB. For each locus there are multiple EP interactions (limited to interactons with significance levels empP<0.05). For each trait and their loci EP distances (outGRB versus inGRB) were tested using Wilcoxon rank sum test, and corrected for multiple testing with bonferroni. The resulting table should indicate traits with LRR in GRBs. 

The list of traits from the GWAS catalog should be annotated and grouped into diseases and traits that are a type of measurement. The latter are traits that contain words such as: measurement, count, volume, index, width, age at, etc. This should be performed using EFO identifiers and searching parent terms that may unambiguously separate the two groups. It's worth checking out if a good grouping can be obtained by searching parent terms of identifeiers material property (these should be disease) or information entitiy (these should be measurements). A possible way is to use [rols library for R]([url](https://bioconductor.org/packages/release/bioc/vignettes/rols/inst/doc/rols.html)).



## sastanak 27.11.2023.
_Dalibor_

- napravljena je analiza pojačivača u GRB-u. za svaku ćeliju iz SOM-a 12x12 i 24 PAM superklase je izračunan postotak pojačivača u GRB-u. Pridruženi su pojedini GRB-i svakom pojačivaču, kao i povezani geni

- rezultati su u /storage/PROJECTS/targPred/dalibor/targPred/GRB_vs_SOM/ 

- sljedeće treba reducirati fantom5, tj. odbaciti sve pojačivače s razinom ekspresije višom od 100, te epidermis jezika kao tkivo, te ponovno istrenirati SOM na takvom skupu
- takodjer treba razdvojiti povezane gene iz GRB-a na *bystander* i *target* gene, te povezati pojačivače s podatcima iz GO

## sastanak 13.11.2023.
_Katarina_

Dalibor SOMs - We decided  that we will proceded with 12x12 PAM and 24 superclasses. Dalibor also prepared bar plots with medain expression numbers per superclass. We noticed that some of the superclasses in the left upper corner are mostly enhancers expressed in immunological cell types and lines. For the next meeting Dalibro will add number of enhancer in the median bar plots. 


## sastanak 02.11.2023.
_Katarina_

Dalibor SOMs - pasted tissue names on the unique IDs. (We should check Katarinas list for duplicate names). Send the clustering script to everyone, prepare the averege expression plots by supercalsses vertically for better visualization. Prepare powerpoint slides wiht same SOM clustering results with different grids and superclustering methods.  
For 8x8 and 12x12 grids 
- more and less superclasses
- different PAM and heierarchical
Also by considering SOM evaluation metrics and the okay results are 8x8 but much better at 10x10.

Katarina - preparing twotailed tests. The results do not show bimodal p-values anymore. Added wilcoxon and t.test but slows down things significantly so find ways to speed this up.

Miha - preparing slides for javni razgovor

## sastanak 19.10.2023.
_Katarina_

We discussed about writing the methods for the paper. Filip has a task to try and write out the workflow from the src/core and any missing things he should ping group members to fill out the unknowns. 

Katarina prepared a EP interaction html with an analysis of EP number of interaction, distance and distribution of p-values. Also, Katarina has a task to make the empP calculation into a two-tailed test. 

## sastanak 10.10.2023.
_Katarina_

Dalibor performed hierarchical clustering of tissue samples from FANTOM5 enhacner datasets using cluster package from R.
He prepared SOM outputs and cells with number of enhancers per cells and also labeld with superclasses. There are a couple of superclasses with a okay number of enhancers. Some superclasses cointain one cell with only one enhancer. Additionally, sorted samples (colnames) by hierarchical clustering and for each superclass we have bar plots that contain median of all enhancers across tissues from superclass. For next time, increase fontsize of all axis labels, y axis to log scale, additionally prepare three different SOM grids 6x6, 8x8 and 12x12 with two different superclustering methods (PAM and hierarchical). 

Filip is making entries to the manuscript. And his todo is to write the methodology for targPred.

Katarina has to fix byGRB list, check why there are samples that appear in multiple pooled samples and then finish her parts from the paper. 

## sastanak 25.09.2023.
_Katarina_

Filip showcased his recent work with transformer models train on genomic sequences. He used (or fine-tuned) [The Nucleotide Transfomrer](https://www.biorxiv.org/content/10.1101/2023.01.11.523679v1) model for classificaton of enhancers and promoters in order to find sequence features that distinct them. He selected layers form the classification model and plotted them with their labels in an interactive scatterplot [Altair](https://altair-viz.github.io/). He further showed a paper [HyenaDNA](https://arxiv.org/abs/2306.15794) which is a ("foundational") transfomrer model that has a large receptive field (1Mb) and a one nucleotide resolution. They have a different approach than The Nucleotide Transformer for adapting the model for different tasks which required fine-tuning, while HyenaDNA does prompt tuning. 

Dalibor is preparing SOMs with enhancers across all tissues from FANTOM5. 

## sastanak 11.09.2023.
_Katarina_

Nastavno na Daliborov JC s paperom:
“Integrative analysis of metabolite GWAS illuminates the molecular basis of pleiotropy and genetic correlation” https://elifesciences.org/articles/79348#content 
Svaki član grupe dobio je dio metoda za obraditi i prezentirati pred ostatkom grupe. Poanta je da bolje razumijemo aktualne metode koje se koriste u za analizu pleiotropije i da razmislimo možemo li nešto iskoristiti za naš projekt. 

@Anja - BOLT-LMM -metoda za GWAS calling,  a mixed-model association method. Inače modeli pretpostavljaju infinitesimal genetic architecture gdje su effect size-ovi normalno distribuirani. BOLT-LMM korsiti Bayesian mixture prior na marker effect size-ove (non-infinitesimal). 
Each locus makes an infinitesimal contribution to the variance
Polygenic = "the random sampling of alleles at each gene produces a continuous, normally distributed phenotype in the population„
Trait isn’t normally distributed, but it’s genetic component is (around parents’ mean)

=> not true for loci under the influence of few genes

loci of large effect conditioned out as fixed effects
Mixed model applied to a selection of markers
Apply Bayesian perspective to mixed model: modelling SNP effects with non-Gaussian prior distributions that better accommodate loci of both small and large effect. -> BOLT-LMM

@Filip Coloc-based colocalization, LDSC genetic correlation i stratified LD score regression

@Dalibor primjena paketa Hess na procjenu nasljednosti svojstava unutar određenih puteva. Lokalna varijanta (ro-Hess) za procjenu lokalne kovarijance LD blokova koje se preklapaju s regijama puteva. Program HESS služi za računanje korelacije između svojstava za određivanje genomskih regija koje pridonose više svojstava (traits).

@Katarina  Mendelian Randomization (MR) za causal inference. MR koristi nasumičnu raspodijelu gena u populaciji kako bi odredila kauzalnost između exposure i outcome. Za MR određuje se set SNP-ova koji su asocirani s exposure varijablom i zadovoljavaju neke uvjete (instrumental variables). Zatim se računa Beta vrijednost (causal effect) kao omjer asocijacije exposure-a i asocijacije s outcome-om. Zaključili smo da trenutno nemamo primjenu ovoga za naš projekt, ali svakako je koristan alat za neka buduća istraživanja.


## sastanak 02.08.2023.
_Katarina_  
@Nina - Raspravljali smo o Nininom zadatku. Nina je pogledala na primjeru jednog VCF-a i pronašla neke overlappove s GWAS datasetovima (oko 90). Isfiltirala je jedan trait, ali na žalost nije bilo overlapping enhancera. Svi vcf file-ovi nalaze se na orthusu i riješen je access. Prije nego što skida fileove kod sebe lokalno će procesirati podatke i sve pacijente staviti u jedan file te će tako izgubiti informacija o pacijentima. Naknadno, kada dobije fenotipkse informacije lako će se ponovno to povezati s pacijentima.
Nakon što pripremi SNP-ove će te varijante provući kroz čitav targPred što uključuje računje LD blokova. To će zahtijevati prilagodbu nekih dijelova koda… 

@Katarina pokazala neke talk-ove koji su potencijalno zanimljivi, poput modela koji predviđa CAGE signal. Jedan poster je pokazao da enformer ne gleda enhancer regije za predviđanje gene expression. Spominjao se [Gasperini et al (2019) CRISPRi dataset]([url](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6690346/)) koji je kandidat za JC. Slajdove i bilješke će biti prebačene u shared drive.   
- [x] Napraviti folder u shareanom drive-u sa slajdovima i bilješkama relevantnih prezentacija s ISCB/ECCB

@svi Početi popunjavati dijelove paper-a. Dodali smo neke stvari, malo izmjenili raspored i izdelegirali paragrafe. Odredili smo 16. Ili 17.08. Za sljedeći sastanak na kojem ćemo raspraviti o update-u. 
- [x]  Do 16. ili 17.08. Napraviti neki entry za paper. [Link na overleaf]([url](https://www.overleaf.com/project/6409fa3589bb0064d5cdfea0))

@Miha Ideja korištenje LLM-a i prompt engineering za izvlačenje dijagnoza i drugih informacija iz free format teksta u neki standardizirani oblik. Isprobati to u kontekstu uzoraka pacijenata i exomskih sekvenci sa MEF-a. Generalizirati za sve nalaze iz BIS-ova. 



## sastanak 10.07.2023.
_Katarina_  
Raspravljali smo kako napraviti som na enhancer-promotors. 
Prvi korak je napraviti SOM na svim enh iz fantom5. Analiza rezultata obuhvaća usporedbu s analizom iz FANTOM5-a. Očekivamo dobiti tissue-specific enhancere i ubiquitous enhancere. Zatim na to overlayati informacije koje imamo o tim enhancerima, njihove karakteristike primjerice empP, spearman. Nadalje, same klastere povezati s promoterima, na način da obilježimo u kojim se klasterima pojavljuju ti enhanceri. 

- [ ] **Dalibor** napraviti SOM na svim enhancer ekspresijskim vrijednostima iz fantom5

S obzirom da su ovo prve bilješke, dodajem stari TODO za koje znam, slobodno popraviti
- [ ] **Filip** Provjeriti je li ova tablica ono što mu je trebalo /storage/PROJECTS/targPred/data/outputs/2021/loci_in_GWAS.df.csv
- [ ] **Katarina** - kmeans na enhancerima i promoterima (bilješke sa sastanka 15.06.)
  Što trebam promjeniti u trenutnoj analizi (manje promjene)
  - [ ] Za E_active i combined_score_mean napraviti i medijan
  - [ ] Oni geni kojima neki transkripti ispadaju iz GRB-a promaknuti u bystandere
  - [x] Provjeriti što je točno is_elite i što znači za naše podatke (is_elite za EP parove znači da je ta interakcija potvrđena barem sa dvije baze, za sumiranje po promotoru ili enhanceru zato ima više smisla gledati mean jer ćemo dobiti točan udio is_elite veza)
  - [ ] Dodati i emP vrijednosti (medijan) i plotat
  - [ ] Slijedeći veći korak je dodati GWAS podatke
  - [ ] Pogledati koji enhancer overlapaju disease
  - [ ] Usporediti genome wide i GRB only
