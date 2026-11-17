# Regjistri i rezultateve

## Gjenerimi i kandidatëve, phiX174 simuluar (seed=42, 30x, 1% gabim)

Parametrat: k=15, min_shared=3, band=5, min_overlap=40, max_occ=100

| Metrika | Vlera |
|---|---|
| Reads | 1292 |
| Çifte brute force | 833,986 |
| Çifte kandidate | 26,555 |
| Zvogëlimi | 31.4x |
| Precision | 0.9997 |
| Recall | 0.9997 |
| Saktësia e orientimit | 1.0 |
| Koha e indeksit | 0.13 s |
| Koha e kandidatëve | 0.56 s |

## Alinjimi i overlap-eve, phiX174 simuluar (seed=42, 30x, 1% gabim)

Parametrat: max_err=0.05, align_band=8, min_overlap=40

| Metrika | Vlera |
|---|---|
| Kandidatë të kontrolluar | 26,555 |
| Overlap-e të pranuara | 24,998 (94.1%) |
| Alinjime Hamming (rruga e shpejtë) | 16,672 |
| Alinjime me DP të brezuar | 9,883 |
| Gjatësia mesatare e overlap-it | 82.83 |
| Error rate mesatar | 0.0181 |
| Quality mesatar | 0.9819 |
| Precision i overlap-eve | 0.9998 |
| Koha | 4.07 s |

Llojet: a_to_b = 11,140; b_to_a = 11,144; b_in_a = 1,406; a_in_b = 1,308

Shënim: DP përdoret kur Hamming > 2.5% të gjatësisë, pra për overlap-e me indel ose shumë gabime.

## Grafi i overlap-it, phiX174 simuluar (seed=42, 30x, 1% gabim)

Parametrat: fuzz=10

| Metrika | Vlera |
|---|---|
| Overlap-e hyrëse | 24,998 |
| Reads të përmbajtura | 856 (66.3%) |
| Reads në graf | 436 |
| Nyje | 872 |
| Skaje para reduktimit | 6,378 |
| Skaje pas reduktimit | 902 (-85.9%) |
| Saktësia e skajeve para / pas | 0.9991 / 0.9933 |
| Koha e reduktimit | 0.007 s |
| Out-degree mesatar / maksimal | 1.034 / 3 |
| Nyje me degëzim | 31 |
| Nyje pa dalje | 5 |
| Komponente të lidhura | 1 |

## Rindërtimi OLC kundrejt greedy, phiX174 simuluar (seed=42, 30x, 1% gabim)

Parametrat: tip_reads=4, bubble_len=6, simplify_rounds=5, fuzz=10

### Thjeshtimi i grafit

| Metrika | Para | Pas |
|---|---|---|
| Skaje | 902 | 832 |
| Nyje me degëzim | 31 | 0 |
| Nyje pa dalje | 5 | 40 |
| Komponente | 1 | 20 |

Tips të hequra: 3 | Bubbles të hequra: 16

Shënim: komponentet rriten sepse reads të hequra nga tips/bubbles mbeten si unitigs të vetme (~140 baza). Ato janë kopje të tepërta dhe përdoren te konsensusi.

### Krahasimi i metodave

| Metrika | OLC | Greedy |
|---|---|---|
| Koha | 0.008 s | 0.005 s |
| Contigs | 20 | 19 |
| Gjatësia totale | 7,794 | 7,739 |
| Contig-u më i gjatë | 5,357 | 1,436 |
| N50 | 5,357 | 757 |
| Mbulimi i gjenomit nga contig-u kryesor | 99.4% | 26.7% |
| Edit distance (contig-u kryesor) | 125 | 33 |
| Identiteti (contig-u kryesor) | 0.9767 | 0.9770 |

Vërejtje: edit distance matet mbi contigs me gjatësi të ndryshme (5357 kundrejt 1436), ndaj identiteti është krahasimi i drejtë.
Përfundim: OLC prodhon një contig të vetëm 3.7 herë më të gjatë se i pari i greedy, me N50 7 herë më të madh. Greedy e fragmenton gjenomin sepse pranon skaje vetëm sipas peshës, pa e parë strukturën e grafit.
Identiteti 97.7% para konsensusit është gabimi i grumbulluar i reads të veçanta; java 6 e korrigjon.

## Konsensusi dhe variantet, phiX174 simuluar (seed=42, 30x, 1% gabim)

Parametrat: polish_rounds=2, polish_band=10, min_alt=3, min_alt_frac=0.25, seq_error=0.01, pvalue_max=1e-6

| Metrika | Vlera |
|---|---|
| Edit distance para konsensusit | 98 |
| Edit distance pas konsensusit | 0 |
| Identiteti pas konsensusit | 100% |
| Gjatësia e rindërtuar | 5,384 |
| Reads të alinjuara | 1,292 / 1,292 (100%) |
| Coverage mesatar | 29.9x |
| Variante të raportuara | 17 (para bashkimit të indels) |
| TP / FP / FN (para bashkimit) | 15 / 2 / 1 |
| TP / FP / FN (pas bashkimit) | 16 / 0 / 0 |
| Koha totale | 8.89 s |
| Memoria maksimale | 80 MB |

Shënim: FP/FN para bashkimit vijnë nga një insertion dy-bazësh (GC, pozicioni 1881) e raportuar si dy insertions një-bazëshe. Sekuenca e rindërtuar është 100% e saktë.

## Eksperimentet (Java 7), phiX174 simuluar (seed=42)

Të gjitha me reads të simuluara nga phiX174 real me 16 mutacione të njohura (10 SNP, 3 INS, 3 DEL), përveç `scale` që përdor gjenome të rastësishëm. Parametrat bazë: k=15, min_shared=3, band=5, min_overlap=40, max_err = max(0.05, 4 × gabimi i simuluar).

### Ablation (30x, 1% gabim)

| Metoda | Contig më i gjatë | Genome fraction | Reads të përfshira | TP / FP / FN | Koha |
|---|---|---|---|---|---|
| OLC i plotë | 5,357 | 0.9994 | 100% | 16 / 0 / 0 | 8.4 s |
| Pa transitive reduction | 151 | 0.0278 | 3.79% | 0 / 0 / 16 | 4.2 s |
| Pa quality scores | 5,357 | 0.9994 | 100% | 16 / 0 / 0 | 7.9 s |
| Greedy (baseline) | 1,436 | 0.2679 | 28.56% | 6 / 0 / 10 | 5.1 s |

Interpretim: transitive reduction është thelbësore. Pa të, grafi ka shumë degëzime dhe contig-u ndahet në çdo degëzim. Greedy gjen vetëm 6 nga 16 mutacionet, konfirmon se nuk është i përshtatshëm. Quality scores nuk ndryshojnë rezultatin në 30x, sepse vota e shumicës i korrigjon gabimet edhe pa peshim.

### Zgjedhja e k (30x, 1% gabim)

| k | Kandidatë | Overlap-e | Reads të përfshira | TP / FP / FN | Koha |
|---|---|---|---|---|---|
| 11 | 26,993 | 25,006 | 100% | 16 / 0 / 0 | 8.3 s |
| 15 | 26,555 | 24,998 | 100% | 16 / 0 / 0 | 7.8 s |
| 21 | 26,008 | 24,745 | 100% | 16 / 0 / 0 | 7.6 s |
| 31 | 23,242 | 22,615 | 99.85% | 16 / 0 / 0 | 6.9 s |

Interpretim: në një gjenom të vogël pa repeats k ka pak ndikim në saktësi. Numri i kandidatëve bie me k, dhe k=31 humbet disa overlap-e (një k-mer i gjatë ka më shumë gjasa të përmbajë gabim).

### Gabimi i sekuencimit (30x)

| Gabimi | Contig më i gjatë | Genome fraction | Reads të përfshira | TP / FP / FN |
|---|---|---|---|---|
| 0.1% | 5,378 | 0.9989 | 100% | 16 / 0 / 0 |
| 1% | 5,357 | 0.9994 | 100% | 16 / 0 / 0 |
| 2% | 5,335 | 0.9974 | 100% | 16 / 0 / 0 |
| 5% | 3,679 | 0.6933 | 69.2% | 11 / 0 / 5 |

Interpretim: me pragun fillestar `max_err = 2.5 × gabimi`, rezultati te 2% ishte i keq (contig 470, TP 1). Shkaku: dy reads me 2% gabim secili kanë rreth 4% mospërputhje mes tyre, afër pragut 5%, ndaj overlap-e të vërteta refuzoheshin. Me `4 × gabimi` problemi u zgjidh (contig 5,335, TP 16). Përfundim: pragu i overlap-it duhet të ndjekë shkallën e gabimit. Te 5% rindërtimi mbulon vetëm 69% të gjenomit, ky është kufiri aktual i algoritmit.

### Coverage (1% gabim)

| Coverage | Reads | Contig më i gjatë | Genome fraction | TP / FP / FN | Koha |
|---|---|---|---|---|---|
| 10x | 430 | 4,149 | 0.7717 | 13 / 0 / 3 | 1.8 s |
| 20x | 861 | 4,157 | 0.7741 | 13 / 0 / 3 | 4.0 s |
| 30x | 1,292 | 5,357 | 0.9994 | 16 / 0 / 0 | 8.6 s |
| 50x | 2,154 | 5,110 | 0.9545 | 16 / 0 / 0 | 18.7 s |

Interpretim: 10x e 20x japin rezultat të ngjashëm dhe jo të plotë. Rezultati i 50x është më i keq se i 30x, që nuk pritet, sepse më shumë reads duhet të ndihmojë. Arsyen nuk e kemi verifikuar ende (do ta hetojmë me diagnostikim).

### Skalueshmëria (gjenom i rastësishëm, coverage 30x, 1% gabim)

| Reads | Gjenomi | Brute force | Kandidatë | Zvogëlimi | Koha | Memoria | Genome fraction |
|---|---|---|---|---|---|---|---|
| 999 | 4.2 kb | 498,501 | 20,722 | 24x | 5.7 s | 115 MB | 0.877 |
| 2,000 | 8.3 kb | 1,999,000 | 41,066 | 49x | 11.5 s | 164 MB | 0.853 |
| 4,999 | 20.8 kb | 12,492,501 | 102,703 | 122x | 24.6 s | 309 MB | 0.572 |
| 10,001 | 41.7 kb | 50,005,000 | 204,915 | 244x | 42.5 s | 553 MB | 0.336 |
| 20,002 | 83.3 kb | 200,030,001 | 409,820 | 488x | 83.1 s | 1,034 MB | 0.313 |

Interpretim: kandidatët rriten afërsisht lineare (rreth 20 për read), ndërsa brute force rritet me N². Koha dhe memoria rriten afërsisht lineare. Kjo është kontributi kryesor i indeksimit me k-mers.

### Kufizime të njohura

- **Plotësia e rindërtimit nuk shkallëzohet:** në gjenomet e mëdha genome fraction bie nga 0.88 në 0.31. Pipeline-i raporton vetëm contig-un më të gjatë, dhe gjenomi ndahet në disa contigs. Shkaku nuk është verifikuar ende.
- **Rezultat jo-monoton te coverage:** 50x jep rezultat më të keq se 30x.
- **Rezultatet janë nga një seed i vetëm (42);** për përfundime më të forta duhen disa seed-a.
- **Metrika `identity`** matet vetëm mbi pjesën e rindërtuar, ndaj duhet lexuar bashkë me `genome_fraction`.
- **Memoria (`peak_memory_mb`)** ka kuptim vetëm te eksperimenti `scale`, sepse te të tjerët është maksimumi i të gjithë procesit.