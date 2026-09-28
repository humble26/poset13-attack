import itertools
from fractions import Fraction
from pathlib import Path
exec(open(Path(__file__).with_name("e5_ideal.py"), encoding="utf-8").read().split("random.seed(5)")[0])
third_=Fraction(1,3)

def full_report(less,n,z):
    Q=[v for v in range(n) if v!=z]
    Ql,idx=induced(less,n,Q); m=n-1
    Dmask=sum(1<<idx[v] for v in Q if less[v][z])
    Umask=sum(1<<idx[v] for v in Q if less[z][v])
    # chain partition (same greedy as e7)
    inc=[[not(Ql[i][j] or Ql[j][i]) for j in range(m)] for i in range(m)]
    order=sorted(range(m), key=lambda v: sum(1 for u in range(m) if Ql[u][v]))
    color=[0]*m
    for v in order:
        preds=[u for u in range(m) if Ql[u][v]]
        color[v]= 1-color[preds[0]] if preds else 0
    A=[v for v in range(m) if color[v]==0]; Bq=[v for v in range(m) if color[v]==1]
    A.sort(key=lambda v: sum(1 for u in range(m) if Ql[u][v]))
    Bq.sort(key=lambda v: sum(1 for u in range(m) if Ql[u][v]))
    aidx={v:k for k,v in enumerate(A)}; bidx={v:k for k,v in enumerate(Bq)}
    pts=[]
    for S in ideals(Ql,m):
        if (S & Dmask)==Dmask and (S & Umask)==0:
            i=sum(1 for k in range(len(A)) if (S>>aidx[A[k]])&1)
            j=sum(1 for k in range(len(Bq)) if (S>>bidx[Bq[k]])&1)
            Ilist=[Q[x] for x in range(m) if (S>>x)&1]
            Clist=[Q[x] for x in range(m) if not (S>>x)&1]
            Il,_=induced(less,n,Ilist); Cl,_=induced(less,n,Clist)
            pts.append((i,j,count_ext(Il,len(Ilist)),count_ext(Cl,len(Clist))))
    eP=sum(E*F for (_,_,E,F) in pts)
    print("  box ideals (i,j,E,F,EF,eF/eP):")
    for (i,j,E,F) in sorted(pts): print("   (%d,%d) E=%d F=%d EF=%d %s"%(i,j,E,F,E*F,Fraction(E*F,eP)))
    # H sweep
    rowmass={}
    for (i,j,E,F) in pts: rowmass[j]=rowmass.get(j,0)+E*F
    H=Fraction(0)
    for j in sorted(rowmass):
        H+=Fraction(rowmass[j],eP)
        print("   H(%d)=%s  rowmass/eP=%s"%(j,H,Fraction(rowmass[j],eP)))
    # z-pair p-values directly
    print("  z-pairs:", end=" ")
    for v in Q:
        if not(less[v][z] or less[z][v]):
            print("(%d,%d):p=%s"%(z,v,pval(less,n,z,v)), end=" ")
    print()

# find an only_col case
found=0
for n in (5,6):
    npairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1<<len(npairs)):
        rel=[npairs[k] for k in range(len(npairs)) if (mask>>k)&1]
        L=closure(n,rel)
        if not is_poset(L,n): continue
        if antichain_width(L,n)<2: continue
        for z in range(n):
            Qz=[v for v in range(n) if v!=z]
            Ql,_=induced(L,n,Qz)
            if antichain_width(Ql,n-1)>2: continue
            r=analyze(L,n,z)
            if r['zpair']: continue
            # recompute heavy lines
            Ql2,idx=induced(L,n,Qz); mm=n-1
            Dmask=sum(1<<idx[v] for v in Qz if L[v][z])
            Umask=sum(1<<idx[v] for v in Qz if L[z][v])
            order=sorted(range(mm), key=lambda v: sum(1 for u in range(mm) if Ql2[u][v]))
            color=[0]*mm
            for v in order:
                pr=[u for u in range(mm) if Ql2[u][v]]
                color[v]= 1-color[pr[0]] if pr else 0
            A=[v for v in range(mm) if color[v]==0]; Bq=[v for v in range(mm) if color[v]==1]
            A.sort(key=lambda v: sum(1 for u in range(mm) if Ql2[u][v]))
            Bq.sort(key=lambda v: sum(1 for u in range(mm) if Ql2[u][v]))
            aidx={v:k for k,v in enumerate(A)}; bidx={v:k for k,v in enumerate(Bq)}
            pts=[]
            for S in ideals(Ql2,mm):
                if (S & Dmask)==Dmask and (S & Umask)==0:
                    i=sum(1 for k in range(len(A)) if (S>>aidx[A[k]])&1)
                    j=sum(1 for k in range(len(Bq)) if (S>>bidx[Bq[k]])&1)
                    Il=[Qz[x] for x in range(mm) if (S>>x)&1]
                    Cl=[Qz[x] for x in range(mm) if not (S>>x)&1]
                    Il2,_=induced(L,n,Il); Cl2,_=induced(L,n,Cl)
                    pts.append((i,j,count_ext(Il2,len(Il)),count_ext(Cl2,len(Cl))))
            eP=sum(E*F for (_,_,E,F) in pts)
            colmass={}; rowmass={}
            for (i,j,E,F) in pts:
                colmass[i]=colmass.get(i,0)+E*F; rowmass[j]=rowmass.get(j,0)+E*F
            hc=any(Fraction(c,eP)>third_ for c in colmass.values())
            hr=any(Fraction(r,eP)>third_ for r in rowmass.values())
            if hc and not hr:
                found+=1
                if found<=2:
                    print("ONLY-COL CASE n=%d z=%d rel=%s"%(n,z,rel))
                    full_report(L,n,z)
                    print()
            if found>=2: break
        if found>=2: break
    if found>=2: break
