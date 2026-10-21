"""Write NZ_CP029256 (from the run's genomic.gbff) and T002_i2 as FASTA references."""
import os, re, shutil, sys
S = '/ebio/abt3_scratch/jmarsh/software/science/snv_build'
GB = '/ebio/abt3_scratch/jmarsh/nanopore/20261012_cm_conjugation/no_sample_id/20261008_1522_P2S-02398-A_PBM28996_eebc9374/genomic.gbff'
os.makedirs(f'{S}/refs', exist_ok=True)
recs, acc, seq, inseq = [], None, [], False
for line in open(GB):
    if line.startswith('VERSION'):
        acc = line.split()[1].split('.')[0]
    elif line.startswith('ORIGIN'):
        inseq, seq = True, []
    elif line.startswith('//'):
        if acc and seq:
            recs.append((acc, ''.join(seq).upper()))
        inseq, acc = False, None
    elif inseq:
        seq.append(re.sub(r'[^A-Za-z]', '', line))
print('gbff records:', [(a, len(s)) for a, s in recs])
main = [r for r in recs if r[0] == 'NZ_CP029256']
assert main, 'NZ_CP029256 not in gbff'
with open(f'{S}/refs/NZ_CP029256.fasta', 'w') as fh:
    fh.write(f'>NZ_CP029256 circular=true\n')
    s = main[0][1]
    for i in range(0, len(s), 80):
        fh.write(s[i:i + 80] + '\n')
shutil.copyfile('/ebio/abt3_scratch/jmarsh/software/science/sv_build/truth/T002_i2.fasta', f'{S}/refs/T002_i2.fasta')
print('refs written')
