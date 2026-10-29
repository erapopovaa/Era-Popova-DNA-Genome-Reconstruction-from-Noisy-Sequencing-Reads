"""Parametrat bazë të pipeline-it (vlerat e verifikuara në eksperimentin e parë me phiX174)."""

# Gjenerimi i kandidatëve
K = 15             # madhësia e k-mer
MAX_OCC = 100      # k-mers me më shumë shfaqje hidhen (repeats)
MIN_SHARED = 3     # k-mers minimalë të përbashkët në të njëjtën diagonale
BAND = 5           # gjerësia e brezit të diagonales (toleranca për indels)
MIN_OVERLAP = 40   # gjatësia minimale e overlap-it

# Simulatori
COVERAGE = 30
ERROR_RATE = 0.01
SEED = 42