import sys
import os
from itertools import product

if len(sys.argv) == 2:
    if os.path.isfile(sys.argv[1]):
        filename = sys.argv[1]
    else:
        sys.exit("Exiting!!! The filename does not exist. Please enter a valid file location.")
else:
    sys.exit("Exiting!!! Please enter a filename as input...")
    
with open(filename, "r") as f:
    value = f.readline().strip().split()
    k, m, n = list(map(float, value))
    
""" k : homozygous dominant (AA)
    m : heterozygous (Aa)
    n : homozygous recessive (aa)

So, we have these 9 possible crosses here,
('AA', 'AA')
('AA', 'Aa')
('AA', 'aa')
('Aa', 'AA')
('Aa', 'Aa')
('Aa', 'aa')
('aa', 'AA')
('aa', 'Aa')
('aa', 'aa')"""

probability = 0
pp = k+m+n
crosses = list(product(['AA', 'Aa', 'aa'], ['AA', 'Aa', 'aa']))

for cross in crosses:
    if cross in [('AA', 'AA')]:
        probability += 1.0 * (k*(k-1))
    elif cross in [('AA', 'Aa'), ('Aa', 'AA')]:
        probability += 1.0 * (k*m)
    elif cross in [('AA', 'aa'), ('aa', 'AA')]:
        probability += 1.0 * (k*n)
    elif cross in [('Aa', 'Aa')]:
        probability += 0.75 * (m*(m-1))
    elif cross in [('Aa', 'aa'), ('aa', 'Aa')]:
        probability += 0.5 * (m*n)
    else:
        probability += 0

print(round(probability/(pp*(pp-1)), 5))