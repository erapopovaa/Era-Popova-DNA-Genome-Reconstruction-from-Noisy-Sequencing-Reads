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