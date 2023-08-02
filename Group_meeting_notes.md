# Intro
U ovom dokumentu navodimo što smo diskutirali na sastancima, zapisujemo ideje i TODO.
Ideja je da Katarina, Filip i Dalibor ovo popunjavaju nakon svakog sastanka. 
Ne treba biti pre detaljno, ali će nam pomoći da ostanemo u fokusu i da ne zaboravimo na dobre ideje :)

## sastanak 02.08.2023.
_Katarina_  
@Nina - Raspravljali smo o Nininom zadatku. Nina je pogledala na primjeru jednog VCF-a i pronašla neke overlappove s GWAS datasetovima (oko 90). Isfiltirala je jedan trait, ali na žalost nije bilo overlapping enhancera. Svi vcf file-ovi nalaze se na orthusu i riješen je access. Prije nego što skida fileove kod sebe lokalno će procesirati podatke i sve pacijente staviti u jedan file te će tako izgubiti informacija o pacijentima. Naknadno, kada dobije fenotipkse informacije lako će se ponovno to povezati s pacijentima.
Nakon što pripremi SNP-ove će te varijante provući kroz čitav targPred što uključuje računje LD blokova. To će zahtijevati prilagodbu nekih dijelova koda… 
@Katarina pokazala neke talk-ove koji su potencijalno zanimljivi, poput modela koji predviđa CAGE signal. Jedan poster je pokazao da enformer ne gleda enhancer regije za predviđanje gene expression. Spominjao se [Gasperini et al (2019) CRISPRi dataset]([url](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6690346/)) koji je kandidat za JC. Slajdove i bilješke će biti prebačene u shared drive.   
- [ ] Napraviti folder u shareanom drive-u sa slajdovima i bilješkama relevantnih prezentacija s ISCB/ECCB

@svi Početi popunjavati dijelove paper-a. Dodali smo neke stvari, malo izmjenili raspored i izdelegirali paragrafe. Odredili smo 16. Ili 17.08. Za sljedeći sastanak na kojem ćemo raspraviti o update-u. 
- [ ]  Do 16. ili 17.08. Napraviti neki entry za paper. [Link na overleaf]([url](https://www.overleaf.com/project/6409fa3589bb0064d5cdfea0))

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
