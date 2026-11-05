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