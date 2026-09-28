import itertools
from fractions import Fraction
from pathlib import Path
exec(open(Path(__file__).with_name("e5_ideal.py"), encoding="utf-8").read().split("random.seed(5)")[0])

third=Fraction(1,3)
stats=dict(total=0, target_fail=0, zpair_fail=0, both_fail=0, rescued=0)
zpair_fail_cases=[]
min_delta=Fraction(2); min_case=None
min_bestz=Fraction(2); minz_case=None

for n in (5,6):
    npairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1<<len(npairs)):
        rel=[npairs[k] for k in range(len(npairs)) if (mask>>k)&1]
        L=closure(n,rel)
        if not is_poset(L,n): continue
        if antichain_width(L,n)<2: continue  # need non-chain P
        for z in range(n):
            Qz=[v for v in range(n) if v!=z]
            Ql,_=induced(L,n,Qz)
            if antichain_width(Ql,n-1)>2: continue
            r=analyze(L,n,z)
            stats['total']+=1
            if not r['target']: stats['target_fail']+=1; print("TARGET FAIL!", n, rel, z)
            if r['delta']<min_delta: min_delta=r['delta']; min_case=(n,rel,z)
            if r['best_z']<min_bestz: min_bestz=r['best_z']; minz_case=(n,rel,z,r['zmax'])
            if not r['zpair']:
                stats['zpair_fail']+=1
                if not r['qpair']: stats['both_fail']+=1
                else: stats['rescued']+=1
                zpair_fail_cases.append((n,rel,z,r))

print("stats:",stats)
print("min delta over class:",min_delta,"case:",min_case)
print("min best-z-pair value:",min_bestz,"case:",minz_case)
print()
for (n,rel,z,r) in zpair_fail_cases[:12]:
    print("zpair-fail: n=%d z=%d zmax=%s nD=%d nU=%d delta=%s best_z=%s rel=%s"%(
        n,z,r['zmax'],r['nD'],r['nU'],r['delta'],r['best_z'],rel))
