# Intro
U ovom dokumentu navodimo što smo diskutirali na sastancima, zapisujemo ideje i TODO.
Ideja je da Katarina, Filip i Dalibor ovo popunjavaju nakon svakog sastanka. 
Ne treba biti pre detaljno, ali će nam pomoći da ostanemo u fokusu i da ne zaboravimo na dobre ideje :)


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
