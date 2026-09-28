"""Reproducible finite part of the norm-five classification (no dependencies)."""
from itertools import product, permutations, combinations
from reduction_model import det_integer, minor_certificate

A49 = ((2,1,1,1),(-1,2,-1,1),(-1,1,2,-1),(-1,-1,1,2))
A7 = ((2,1,1,1),(-1,2,-1,1),(-1,1,2,-1),(1,1,-1,2))
A17 = ((2,1,1,1),(-1,2,-1,1),(1,-1,2,1),(1,1,-1,2))
REPS = {49:A49, 7:A7, 17:A17}
B = tuple(A49[i] for i in (1,0,3,2))

def observability_matrix(M,j):
    """Rows e_j^T M^k, k=0,...,n-1 (single nonlinear coordinate j)."""
    n=len(M); O=[tuple(int(i==j) for i in range(n))]
    for k in range(n-1):
        O.append(tuple(sum(O[-1][i]*M[i][s] for i in range(n)) for s in range(n)))
    return O

def normal_form(a,b,e,h):
    return ((2,1,1,1),(a,2,b,-b),(-a*b*e,e,2,-e),(a*b*h,h,-h,2))

def normalized_orbit(M):
    # After aligning unique +/-2 entries with the diagonal and making it +2,
    # all residual signed row/column permutations are simultaneous permutations
    # and sign conjugations. First-row normalization removes the latter.
    out=set()
    for pi in permutations(range(4)):
        B=tuple(tuple(M[pi[i]][pi[j]] for j in range(4)) for i in range(4))
        d=(1,)+B[0][1:]
        out.add(tuple(tuple(d[i]*B[i][j]*d[j] for j in range(4)) for i in range(4)))
    return out

def characteristic_polynomial(M):
    # Coefficients in increasing degree, from the Leibniz formula.
    n=len(M); out=[0]*(n+1)
    for pi in permutations(range(n)):
        terms=[(-1)**sum(pi[i]>pi[j] for i in range(n) for j in range(i+1,n))]
        for i,j in enumerate(pi):
            ans=[0]*(len(terms)+1)
            for k,x in enumerate(terms):
                ans[k]-=x*M[i][j]
                if i==j:ans[k+1]+=x
            terms=ans
        out=[a+b for a,b in zip(out,terms)]
    return out

def prime_divisors(n):
    n=abs(n); out=[]; d=2
    while d*d<=n:
        if n%d==0:
            out.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:out.append(n)
    return out

def classify():
    all_normal=set()
    positions=[(i,j) for i in range(1,4) for j in range(4) if i!=j]
    for signs in product((-1,1),repeat=9):
        A=[[2 if i==j else 1 for j in range(4)] for i in range(4)]
        for (i,j),s in zip(positions,signs):A[i][j]=s
        A=tuple(map(tuple,A))
        cert=minor_certificate(A)
        if all(m['det']!=0 for m in cert if len(m['rows'])==2):
            all_normal.add(A)
    parametric={normal_form(*s) for s in product((-1,1),repeat=4)}
    assert all_normal==parametric and len(all_normal)==16
    assert all(all(m['det'] for m in minor_certificate(M)) for M in all_normal)
    classes=[]; seen=set()
    for det,M in REPS.items():
        orbit=normalized_orbit(M)
        assert not (seen&orbit)
        assert all(abs(det_integer(B))==det for B in orbit)
        seen|=orbit
        cert=minor_certificate(M)
        classes.append(dict(label=f'A{det}',matrix=M,det=det_integer(M),
                            normalized_members=len(orbit),
                            minor_abs={k:sorted({abs(m['det']) for m in cert if len(m['rows'])==k}) for k in range(1,5)},
                            bad_primes=sorted({p for m in cert for p in prime_divisors(m['det'])}),
                            characteristic_polynomial=characteristic_polynomial(M),
                            certificate=cert))
    assert seen==all_normal
    return dict(normalized_sign_patterns_examined=512,normal_forms=16,classes=classes,
                parameter_table=[dict(signs=s,det=det_integer(normal_form(*s)),
                                      representative=f'A{abs(det_integer(normal_form(*s)))}')
                                 for s in product((-1,1),repeat=4)])

def classify_positive():
    """Complete finite enumeration for positive integer row sums at most seven.

    Distinct rows are required. Rows are sorted to quotient row permutations;
    a second canonicalization quotients column permutations.
    """
    rows=[r for r in product(range(1,5),repeat=4) if sum(r)<=7]
    compat={(i,j):all(rows[i][a]*rows[j][b]!=rows[i][b]*rows[j][a]
                     for a,b in combinations(range(4),2))
            for i,j in combinations(range(len(rows)),2)}
    classes={}; survivors=[]
    for ids in combinations(range(len(rows)),4):
        if not all(compat[i,j] for i,j in combinations(ids,2)):continue
        M=tuple(rows[i] for i in ids)
        cert=minor_certificate(M)
        if not all(m['det'] for m in cert):continue
        survivors.append(M)
        key=min(tuple(sorted(tuple(row[j] for j in pi) for row in M))
                for pi in permutations(range(4)))
        classes.setdefault(key,[]).append(M)
    assert len(rows)==35 and len(survivors)==30 and len(classes)==4
    records=[]
    for M,members in classes.items():
        cert=minor_certificate(M)
        records.append(dict(matrix=M,det=det_integer(M),
                            row_unordered_members=len(members),
                            labeled_members=24*len(members),
                            row_sums=sorted(map(sum,M)),
                            minor_abs={k:sorted({abs(m['det']) for m in cert if len(m['rows'])==k}) for k in range(1,5)},
                            bad_primes=sorted({p for m in cert for p in prime_divisors(m['det'])}),
                            certificate=cert))
    return dict(row_candidates=35,row_unordered_matrices_examined=52360,
                rational_MDS_survivors=30,classes=records)

if __name__=='__main__':
    import json
    print(json.dumps(classify(),indent=2))
