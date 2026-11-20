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
Shënim: rezultat i hershëm, me max_err=0.05. Rezultatet përfundimtare janë te seksioni "Eksperimentet (Java 7)". 

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

Reads të simuluara nga phiX174 real me 16 mutacione të njohura (10 SNP, 3 INS, 3 DEL), përveç `scale` që përdor gjenome të rastësishëm. Parametrat: k=15, min_shared=3, band=5, min_overlap=40, align_band=8, tip_reads=4, bubble_len=15, max_err = max(0.08, 4 × gabimi i simuluar).

### Diagnostikimi i fragmentimit

Në rezultatet e para, contig-u ndahej te coverage 50x dhe te gjenomet e mëdha (genome fraction 0.95 te 50x dhe 0.31 te 20k reads). Hetimi me `experiments/diagnose.py`:
1. Mbulimi i vërtetë në vendet e ndarjes ishte 22-48x, pra nuk ishte problem i të dhënave.
2. Pas thjeshtimit mbeteshin 4-6 nyje me degëzim, që e priznin contig-un.
3. Skajet që mungonin nuk ishin hequr nga transitive reduction, por nuk ekzistonin që në fazën e overlap-it.
4. Çiftet ishin kandidatë të gjetur nga indeksi (49 dhe 29 vota), por alinjimi jepte error rate 0.053 dhe 0.060, mbi pragun `max_err = 0.05`. Brezi i alinjimit (8 apo 16) nuk ndryshonte rezultatin.
5. Ngritja e `max_err` në 0.08 zgjidhi degëzimet (4 → 0), komponentet (20 → 1) dhe prodhoi një contig të vetëm 5350 baza te phix50.
6. `bubble_len` 6 → 15 ndihmoi pjesërisht (contig-u më i gjatë te scale5000: 11859 → 17531) para se të korrigjohej pragu. `tip_reads` 4 → 10 nuk ndryshoi asgjë.

Përfundim: me 1% gabim për bazë, overlap-et e vërteta kanë luhatje të rastit dhe disa kalojnë 5%. Pragu duhet vendosur me diferencë mbi gabimin e pritur.

### Ablation (30x, 1% gabim)

| Metoda | Contig më i gjatë | Genome fraction | Reads të përfshira | TP / FP / FN | Koha |
|---|---|---|---|---|---|
| OLC i plotë | 5,357 | 0.9994 | 100% | 16 / 0 / 0 | 7.7 s |
| Pa transitive reduction | 151 | 0.0278 | 3.79% | 0 / 0 / 16 | 4.1 s |
| Pa quality scores | 5,357 | 0.9994 | 100% | 16 / 0 / 0 | 7.8 s |
| Greedy (baseline) | 5,357 | 0.9994 | 100% | 16 / 0 / 0 | 7.8 s |

Interpretim: transitive reduction është thelbësore (pa të vetëm 151 baza). Quality scores nuk ndryshojnë rezultatin në 30x. Greedy del njësoj si OLC në këtë gjenom. Me pragun e vjetër 0.05 greedy jepte vetëm 1436 baza dhe 6/16 mutacione, por kjo ishte pasojë e skajeve që mungonin, jo e algoritmit. Hipoteza e mundshme: phiX174 nuk ka repeats, ndaj grafi është pothuajse linear dhe greedy mjafton.Testuar me vone te seksioni Repeats:me reads te gjata greedy dhe olc_best japin te njejtin rezultat,me reads te shkurtra olc eshte me i sakte

### Zgjedhja e k (30x, 1% gabim)

| k | Kandidatë | Overlap-e | Genome fraction | Reads të përfshira | TP / FP / FN | Koha |
|---|---|---|---|---|---|---|
| 11 | 26,993 | 26,465 | 0.9994 | 100% | 16 / 0 / 0 | 7.9 s |
| 15 | 26,555 | 26,428 | 0.9994 | 100% | 16 / 0 / 0 | 7.7 s |
| 21 | 26,008 | 25,949 | 0.9994 | 100% | 16 / 0 / 0 | 7.5 s |
| 31 | 23,242 | 23,230 | 0.5879 | 58.9% | 7 / 0 / 9 | 5.3 s |

Interpretim: k = 11-21 japin rezultat të njëjtë. K=31 dëmtohet sepse një k-mer i gjatë ka më shumë gjasa të përmbajë gabim, ndaj çifte të vërteta humbasin dhe grafi ndahet. Zgjedhja k=15 është e justifikuar.

### Gabimi i sekuencimit (30x)

| Gabimi | Contig më i gjatë | Genome fraction | Reads të përfshira | TP / FP / FN |
|---|---|---|---|---|
| 0.1% | 5,378 | 0.9989 | 100% | 16 / 0 / 0 |
| 1% | 5,357 | 0.9994 | 100% | 16 / 0 / 0 |
| 2% | 5,335 | 0.9974 | 100% | 16 / 0 / 0 |
| 5% | 3,679 | 0.6933 | 69.2% | 11 / 0 / 5 |

Interpretim: pragu i overlap-it duhet të ndjekë gabimin (4 × gabimi; me pragun 2.5× rezultati te 2% ishte contig 470 dhe TP 1). Te 5% rindërtimi mbulon 69% të gjenomit, ky është kufiri aktual. Nuk e kemi hetuar shkakun te 5%.

### Coverage (1% gabim)

| Coverage | Reads | Contig më i gjatë | Genome fraction | TP / FP / FN | Koha |
|---|---|---|---|---|---|
| 10x | 430 | 5,360 | 0.9972 | 16 / 0 / 0 | 2.1 s |
| 20x | 861 | 5,361 | 0.9994 | 16 / 0 / 0 | 4.3 s |
| 30x | 1,292 | 5,357 | 0.9994 | 16 / 0 / 0 | 7.6 s |
| 50x | 2,154 | 5,350 | 0.9996 | 16 / 0 / 0 | 17.0 s |

Interpretim: rindërtimi është i plotë që nga 10x. Rezultatet e vjetra jo-monotone (50x më keq se 30x) vinin nga pragu 0.05.

### Skalueshmëria (gjenom i rastësishëm, 30x, 1% gabim)

| Reads | Gjenomi | Brute force | Kandidatë | Zvogëlimi | Koha | Memoria | Genome fraction | TP / FP / FN |
|---|---|---|---|---|---|---|---|---|
| 999 | 4.2 kb | 498,501 | 20,722 | 24x | 6.0 s | 115 MB | 0.9990 | 16 / 0 / 0 |
| 2,000 | 8.3 kb | 1,999,000 | 41,066 | 49x | 11.9 s | 164 MB | 0.9986 | 26 / 0 / 0 |
| 4,999 | 20.8 kb | 12,492,501 | 102,703 | 122x | 30.7 s | 311 MB | 0.9996 | 67 / 0 / 0 |
| 10,001 | 41.7 kb | 50,005,000 | 204,915 | 244x | 63.2 s | 553 MB | 0.9996 | 133 / 8 / 4 |
| 20,002 | 83.3 kb | 200,030,001 | 409,820 | 488x | 132.5 s | 1,034 MB | 0.9998 | 271 / 8 / 5 |

Interpretim: kandidatët rriten afërsisht lineare (rreth 20 për read), brute force rritet me N². Koha dhe memoria rriten afërsisht lineare. Rindërtimi mbetet i plotë (genome fraction mbi 0.998) deri te 20k reads.

### Kufizime të njohura

- **8 FP** te 10k dhe 20k reads (rreth 3% e varianteve) dhe 4-5 FN. Shkaku nuk është verifikuar.
- **K=31** dëmton rindërtimin (genome fraction 0.59).
- **Gabimi 5%** jep genome fraction 0.69.
- **Repeats:** testuar te seksionet e fundit (gjenom sintetik, repeat me 3 kopje identike). Rezultatet varen nga gjatësia e reads kundrejt repeat-it.
- Rezultatet e seksionit Java 7 janë nga një seed i vetëm (42); eksperimentet me repeats përdorin 3 deri 10 seed-a.

## Eksperimenti me repeats (gjenom i rastësishëm 10 kb, 30x, 1% gabim, 3 seed-a: 42, 43, 44)

Gjenomi: 10,000 baza të rastësishme; segmenti 2000..2000+R kopjohet te pozicionet 5000 dhe 8000 (3 kopje gjithsej). Reads 100-150 bp. Metodat: `olc` (transitive reduction + simplify), `olc_best` (olc + best-overlap `resolve_branches`), `greedy`. Skripti: `experiments/repeats.py`; të dhënat: `results/benchmark_repeats.csv`.

| R | Metoda | Contig-e | Contig më i gjatë | Genome fraction | Edit distance | Koha |
|---|---|---|---|---|---|---|
| 0 | olc | 1.3 | 9943 | 0.9986 | 0.0 | 14.8 s |
| 0 | olc_best | 1.3 | 9616 | 0.9658 | 0.0 | 14.3 s |
| 0 | greedy | 4.0 | 9616 | 0.9658 | 0.0 | 14.4 s |
| 80 | olc | 11.3 | 3068 | 0.3081 | 509.0 | 10.8 s |
| 80 | olc_best | 1.3 | 9616 | 0.9658 | 0.0 | 15.3 s |
| 80 | greedy | 4.0 | 9616 | 0.9658 | 0.0 | 15.4 s |
| 300 | olc | 6.7 | 2913 | 0.2926 | 1.0 | 12.3 s |
| 300 | olc_best | 3.3 | 4451 | 0.4471 | 1170.7 | 13.1 s |
| 300 | greedy | 5.7 | 4457 | 0.4476 | 705.7 | 13.1 s |
| 1000 | olc | 6.7 | 2229 | 0.2238 | 3.0 | 14.5 s |
| 1000 | olc_best | 3.3 | 4500 | 0.4522 | 701.0 | 23.6 s |
| 1000 | greedy | 6.7 | 3778 | 0.3797 | 84.7 | 16.9 s |

Mesataret janë mbi 3 seed-a. Gjenom fraction për seed-in individual te R=0 / R=80 për `olc_best`: 0.9987 / 0.9006 / 0.998 (seed 42 / 43 / 44).

### Çfarë është e provuar

1. **R=80:** `olc` fragmentohet (0.31, 10-13 contig-e), ndërsa `olc_best` e rikthen gjenomin njësoj si te R=0. Problemi i `olc` ishte degëzimi në graf te repeat-i dhe `SIMPLIFY_ROUNDS=0` nuk e ndryshoi rezultatin.
2. **R=300 dhe R=1000:** asnjë metodë nuk e rindërton gjenomin. Genome fraction 0.22-0.45. Edit distance i lartë te `olc_best` (1170 te R=300) dhe greedy (705) tregon se contig-u më i gjatë përmban gabime strukturore, ndaj fraction nuk lexohet si përmirësim kundrejt `olc` (0.29).
3. **`olc_best` krahasuar me greedy:** rezultate të ngjashme te R=300 (0.4471 kundrejt 0.4476) dhe R=1000 (0.4522 kundrejt 0.3797, me edit 701 kundrejt 85).
4. **Gabim i vërtetë i gjetur dhe i rregulluar:** versioni i parë i `resolve_branches` llogariste best-in mbi grafin para se të hiqeshin skajet e best-out, dhe mund të linte një nyje pa asnjë hyrje. Testi `tests/test_resolve.py` e kap gabimin: versioni i vjetër heq skajin e vlefshëm (u->v), versioni me dy kalime e ruan. Ky rregullim nuk e ndryshoi seed 43.

### Seed 43, R=0: pse `olc_best` jep 0.9006 (diagnostikuar me `diagnose_best.py` dhe `diagnose_pair.py`)

- `resolve_branches` hoqi vetëm 1 skaj nga 1608 (me binjakun 2 në graf; grafi mbeti simetrik).
- `olc_best` jep 2 contig-e që mbivendosen ~83 bp (pozicione përafërsisht 3-8971 dhe 8888-9990): **gjenomi mbulohet plotësisht**. Rënia 0.999 -> 0.9006 vjen nga metrika që merr vetëm contig-un më të gjatë.
- Skaji i hequr (`read_13 -> read_1993`, overlap 112 bp) ishte **i vërtetë**. Humbi kundër `read_735 -> read_1993` (score 120 kundrejt 109).
- Skaji `read_13 <-> read_735` mungon sepse alinjimi jep edit distance 12 në ~121 bp (raporti 0.099 > MAX_ERR 0.08). `read_735` ka 10 gabime në ~121 bp (~8%) kundrejt 1% të simuluar, pra read jashtëzakonisht i zhurmshëm.
- Përfundim: nuk është bug, por kufizim i best-overlap: kur zgjedh vetëm një skaj për nyje, humbet rrugën nëse fituesi ka hyrje të dobët.

### Hipoteza të pa provuara

- Kufiri te R≥300 është gjatësia e read-eve (100-150 bp, më e shkurtër se repeat-i plus flanqet). Test i propozuar: R=300 me reads 500-700 bp.
- `olc` nuk thyhet te seed 43 sepse `simplify` e trajton `read_735` si tip dhe e hedh. Nuk është verifikuar.

### Kufizime të këtij eksperimenti

- Vetëm 3 seed-a dhe gjenom sintetik i rastësishëm (pa regjione me kompleksitet të ulët).
- Metrika: contig më i gjatë / gjenomi; `pipeline` poliron vetëm `contigs[0]`. Mbulimi i bashkuar i të gjithë contig-eve ose N50 do të ishte më i drejtë.
- `edit_distance` matet mbi contig-un e parë, ndaj nuk është i krahasueshëm mes rreshtave të fragmentuar.
- `olc_best` i shkurtër te seed 43 (R=0 dhe R=80) është i njëjtë me greedy; pra te ky seed best-overlap dhe greedy kanë të njëjtin kufizim.

## Eksperimentet me repeats (gjenom i rastësishëm 10 kb, 30x, 1% gabim)

Gjenomi: 10,000 baza të rastësishme; segmenti 2000..2000+R kopjohet te pozicionet 5000 dhe 8000 (3 kopje identike gjithsej). Metodat:
- `olc`: transitive reduction + simplify (tips dhe bubbles);
- `olc_best`: `olc` + `resolve_branches` (best-overlap në dy kalime);
- `greedy`: baseline.

Skriptet: `experiments/repeats.py`, `experiments/read_length.py`, `experiments/repeats_seeds.py`. Të dhënat: `results/benchmark_repeats.csv`, `results/benchmark_readlen.csv`, `results/benchmark_repeats_seeds.csv`.

### Metrikat

- **covered_fraction:** pjesa e gjenomit e mbuluar nga bashkimi i të gjitha contig-eve që vendosen në referencë me identitet >= 0.95 (`assembly_metrics` te `src/evaluate.py`). Contig-et me identitet nën 0.95 (misassembly, repeat i kolapsuar, ose i pavendosur) nuk numërohen. Pragu 0.95 është zgjedhje e jona (contig-et e papoliruara kanë ~2% gabim), jo vlerë e matur. Një contig pjesërisht i saktë por me një misassembly hidhet i tëri, ndaj metrika mund të nënvlerësojë pak.
- **contigs_bad:** numri i contig-eve që nuk plotësojnë kushtin e mësipërm.
- **N50 / NA50:** N50 mbi të gjithë contig-et; NA50 vetëm mbi contig-et e sakta.
- **genome_fraction** (contig më i gjatë / gjenomi) dhe **edit_distance** matin vetëm `contigs[0]`. Genome fraction shpërblen contig-e të gjata por të gabuara dhe ndëshkon rindërtime të ndara në contig-e që mbivendosen. Prandaj tabelat më poshtë nuk i përdorin dhe interpretimi bazohet te metrikat e mësipërme.
- `time_seconds` ndryshon nga një ekzekutim te tjetri (p.sh. `olc` R=0: 24.3 s në një ekzekutim, 14.9 s te `olc_best`, ndërsa algoritmi është thelbësisht i njëjtë), ndaj nuk përdoret për të krahasuar metodat.

### 1. Repeats, reads 100-150 bp, 3 seed-a (42, 43, 44)

| R | Metoda | Contig-e | covered_fraction | N50 | NA50 | contigs_bad |
|---|---|---|---|---|---|---|
| 0 | olc | 1.3 | 0.9943 | 9943 | 9943 | 0.3 |
| 0 | olc_best | 1.3 | 0.9951 | 9616 | 9616 | 0.0 |
| 0 | greedy | 4.0 | 0.9950 | 9616 | 9616 | 0.0 |
| 80 | olc | 11.3 | 0.9977 | 3053 | 3053 | 0.3 |
| 80 | olc_best | 1.3 | 0.9951 | 9616 | 9616 | 0.0 |
| 80 | greedy | 4.0 | 0.9950 | 9616 | 9616 | 0.0 |
| 300 | olc | 6.7 | 0.9815 | 2898 | 2898 | 0.3 |
| 300 | olc_best | 3.3 | 0.6961 | 2910 | 2644 | 0.7 |
| 300 | greedy | 5.7 | 0.5460 | 2909 | 2643 | 1.0 |
| 1000 | olc | 6.7 | 0.8419 | 2135 | 2170 | 0.3 |
| 1000 | olc_best | 3.3 | 0.6699 | 4204 | 3163 | 0.3 |
| 1000 | greedy | 6.7 | 0.7561 | 3240 | 3195 | 0.3 |

Për R=300 me 10 seed-a, mesatarja e greedy del 0.66 (jo 0.55); me 3 seed-a mesataret janë të paqëndrueshme (shih tabelën 3).

### 2. Gjatësia e reads, 3 seed-a (42, 43, 44)

| R | Reads | Metoda | Contig-e | covered_fraction | N50 | NA50 | contigs_bad |
|---|---|---|---|---|---|---|---|
| 300 | 100-150 | olc | 6.7 | 0.9815 | 2898 | 2898 | 0.3 |
| 300 | 100-150 | olc_best | 3.3 | 0.6961 | 2910 | 2644 | 0.7 |
| 300 | 100-150 | greedy | 5.7 | 0.5460 | 2909 | 2643 | 1.0 |
| 300 | 500-700 | olc | 7.3 | 0.9962 | 2248 | 2248 | 0.0 |
| 300 | 500-700 | olc_best | 1.0 | 0.9936 | 9936 | 9936 | 0.0 |
| 300 | 500-700 | greedy | 1.0 | 0.9936 | 9936 | 9936 | 0.0 |
| 1000 | 100-150 | olc | 6.7 | 0.8419 | 2135 | 2170 | 0.3 |
| 1000 | 100-150 | olc_best | 3.3 | 0.6699 | 4204 | 3163 | 0.3 |
| 1000 | 100-150 | greedy | 6.7 | 0.7561 | 3240 | 3195 | 0.3 |
| 1000 | 1200-1500 | olc | 10.0 | 0.9939 | 2267 | 2267 | 0.0 |
| 1000 | 1200-1500 | olc_best | 1.0 | 0.9918 | 9918 | 9918 | 0.0 |
| 1000 | 1200-1500 | greedy | 1.3 | 0.9919 | 9919 | 9919 | 0.0 |

Koha e `olc_best` me reads 1200-1500 (183 s mesatarisht, kundrejt 50 s të greedy) është e pashpjeguar; nuk është diagnostikuar.

### 3. Repeats, reads 100-150 bp, 10 seed-a (42 deri 51), mesatare ± devijim standard

| R | Metoda | covered_fraction | contigs_bad | N50 | NA50 | Seed-a me covered >= 0.95 |
|---|---|---|---|---|---|---|
| 300 | olc | 0.9831 ± 0.0033 | 0.3 ± 0.48 | 2905 | 2905 | 10 / 10 |
| 300 | olc_best | 0.6603 ± 0.1822 | 1.0 ± 0.67 | 2912 | 2563 | 2 / 10 |
| 300 | greedy | 0.6611 ± 0.1836 | 1.0 ± 0.47 | 2914 | 2670 | 2 / 10 |
| 1000 | olc | 0.8443 ± 0.0031 | 0.3 ± 0.48 | 2148 | 2159 | 0 / 10 |
| 1000 | olc_best | 0.6703 ± 0.2280 | 0.6 ± 0.84 | 3706 | 3355 | 0 / 10 |
| 1000 | greedy | 0.7404 ± 0.1743 | 0.5 ± 0.71 | 3418 | 3391 | 0 / 10 |

### Çfarë është e provuar (sipas të dhënave më sipër)

1. **Gjatësia e reads vendos nëse repeat-i zgjidhet.** Kur reads janë më të gjata se repeat-i plus flanqet (500-700 bp për R=300; 1200-1500 bp për R=1000), `olc_best` dhe greedy rindërtojnë gjenomin me 1 contig, mbulim 0.99+ dhe 0 contig të gabuar në të 3 seed-at. Me reads 100-150 bp asnjë metodë nuk e rindërton gjenomin kur R >= 300 (0/10 seed-a mbi 0.95 te R=1000).
2. **Me reads më të shkurtra se repeat-i, `olc` është më i besueshëm se `olc_best` dhe greedy.** Te R=300 (10 seed-a): `olc` mbulim i saktë 0.983 ± 0.003 dhe 10/10 seed-a mbi 0.95; `olc_best` dhe greedy 0.66 ± 0.18 dhe vetëm 2/10. Dallimi është shumë më i madh se luhatja mes seed-ave. Rezultati i `olc_best` është dypjesësh (rreth 0.98 ose 0.50-0.68), ndaj mesatarja nuk e përfaqëson asnjë rast tipik.
3. **Me reads të gjata, greedy dhe `olc_best` japin rezultate të njëjta.** Kjo ka kuptim: të dyja zgjedhin skajin me score më të lartë për çdo degëzim. `olc` është konservator: nuk e lidh repeat-in dhe mbetet i saktë (0 contig të gabuar, mbulim 0.99+), por i thyer (N50 rreth 2250).
4. **Te R=1000 me reads të shkurtra nuk ka renditje të besueshme mes `olc_best` (0.67) dhe greedy (0.74):** devijimi standard (0.17-0.23) është më i madh se diferenca. Vetëm `olc` është i qëndrueshëm (0.844 ± 0.003).
5. **R=80 (repeat më i shkurtër se reads):** `olc` jep 11 contig-e të sakta (mbulim 0.9977, N50 3053); `olc_best` dhe greedy japin 1 deri 4 contig-e të gjata. Reads e kalojnë repeat-in, ndaj best-overlap e lidh pa gabime (`contigs_bad` 0).
6. **Gabim i vërtetë i gjetur dhe i rregulluar:** versioni i parë i `resolve_branches` llogariste best-in mbi grafin para se të hiqeshin skajet e best-out, dhe mund të linte një nyje pa hyrje. Testi `tests/test_resolve.py` e kap gabimin: versioni i vjetër heq skajin e vlefshëm (u->v), versioni me dy kalime e ruan (verifikuar duke ekzekutuar të dy versionet).

### Seed 43, R=0: pse `genome_fraction` i `olc_best` ishte 0.9006 (diagnostikuar me `diagnose_best.py` dhe `diagnose_pair.py`)

- `resolve_branches` hoqi vetëm 1 skaj nga 1608 (me binjakun, 2 në graf; grafi mbeti simetrik).
- `olc_best` jep 2 contig-e që mbivendosen ~83 bp (afërsisht 3-8971 dhe 8888-9990): gjenomi mbulohet plotësisht (covered_fraction 0.9951, `contigs_bad` 0). Rënia 0.999 -> 0.9006 ishte artefakt i metrikës "contig më i gjatë".
- Skaji i hequr (`read_13 -> read_1993`, overlap 112 bp) ishte i vërtetë. Humbi kundër `read_735 -> read_1993` (score 120 kundrejt 109).
- Skaji `read_13 <-> read_735` mungon sepse alinjimi jep edit distance 12 në ~121 bp (raporti 0.099 > `MAX_ERR` 0.08). `read_735` ka 10 gabime në ~121 bp (~8%) kundrejt 1% të simuluar.
- Përfundim: kufizim i best-overlap me një read shumë të zhurmshëm, jo bug.

### Hipoteza të pa provuara

- `olc_best` dhe greedy dështojnë te reads të shkurtra sepse zgjedhin një degë te repeat-i dhe e kolapsojnë; kjo përputhet me `contigs_bad` dhe covered_fraction, por nuk është verifikuar duke shqyrtuar vetë contig-et.
- `olc` te R=0 ka `contigs_bad` 0.3 (kontig 137 bp i seed 43, që mendojmë se është `read_735`). Nuk është verifikuar.
- Pse `simplify` e lë repeat-in të thyer ndërsa `resolve_branches` jo: nuk është diagnostikuar.
- Pse `olc` te R=1000 arrin vetëm 0.84 mbulim (e kemi parë, por shkaku nuk është hetuar).
- `MAX_OCC=100` te indeksi nuk është testuar drejtpërdrejt.

### Kufizime të këtyre eksperimenteve

- Gjenom sintetik i rastësishëm, 3 kopje identike të repeat-it (pa divergjencë, siç ka shpesh te gjenomet reale), pa regjione me kompleksitet të ulët.
- Model i thjeshtë gabimi (1%, 90% zëvendësime); reads të simuluara, jo reale.
- Tabelat 1 dhe 2: 3 seed-a (mesataret janë të zhurmshme, sidomos te R=1000); tabela 3: 10 seed-a vetëm për reads të shkurtra.
- Pragu i identitetit 0.95 në `assembly_metrics` është zgjedhje e jona.
- `edit_distance` dhe `genome_fraction` nga `benchmark_*.csv` mbeten të përfshira te skedarët CSV, por matin vetëm contig-un e parë dhe nuk duhen përdorur për krahasim.