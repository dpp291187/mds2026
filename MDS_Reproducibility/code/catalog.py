"""Complete reduced cancellation-free catalog and a bounded-load synthesis certificate.

The catalog is exhaustive ONLY in the grammar documented in the article:
unsigned gate identities are coefficient vectors modulo global sign; each vector
has one producer; no dead nodes or coefficient cancellation; binary +/- gates
and doublings; at most 13 arithmetic gates. It is NOT all integer SLPs.
The separate five-normalization certificate permits cancellation and output reuse.
"""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
import json, hashlib
from reduction_model import signed_m4, circuit, optimize, brute, assess
ROOT=Path(__file__).resolve().parents[1]
A=signed_m4().matrix()
PROBE=(10000,10001,10010,10100)


def l1(v): return sum(map(abs,v))

def canonical(v):
    s=1 if next(x for x in v if x)>0 else -1
    return tuple(s*x for x in v),s


def subforms(v):
    for aa in product(*(range(abs(x)+1) for x in v)):
        a=tuple((1 if x>=0 else -1)*z for x,z in zip(v,aa))
        if any(a): yield a


@lru_cache(None)
def trees(v):
    if l1(v)==1: return (frozenset(),)
    out=set()
    for a in subforms(v):
        b=tuple(x-y for x,y in zip(v,a))
        if not any(b) or a>b: continue
        ac,sa=canonical(a); bc,sb=canonical(b)
        gate=tuple(sorted(((ac,sa),(bc,sb))))
        for ta in trees(ac):
            for tb in trees(bc):
                d=dict(ta)
                if any(k in d and d[k]!=g for k,g in tb): continue
                d.update(tb);d[v]=gate
                out.add(frozenset(d.items()))
    return tuple(sorted(out,key=lambda t:sorted(t)))


def enumerate_catalog():
    roots=[canonical(r)[0] for r in A]
    families=[trees(r) for r in roots]
    assert list(map(len,families))==[60]*4
    combos={frozenset()};stage_counts=[]
    for i,ls in enumerate(families):
        nxt=set()
        for t in combos:
            d=dict(t)
            for tt in ls:
                if any(k in d and d[k]!=g for k,g in tt): continue
                z=t|tt
                if len(z)+3-i<=13:nxt.add(z)
        combos=nxt;stage_counts.append(len(combos))
    return sorted(combos,key=lambda t:sorted(t)),stage_counts


def materialize(gates,index):
    # Positive evaluations orient every wire. Because every A-row evaluates
    # positively at PROBE, each output has exactly the required sign, and
    # no gate has two negative operands. Thus all gates are actual +/- or 2u.
    ori=lambda v: 1 if sum(a*b for a,b in zip(v,PROBE))>0 else -1
    names={tuple(int(i==j) for i in range(4)):f'x{j}' for j in range(4)}
    ops=[]
    for v,children in sorted(gates,key=lambda g:(l1(g[0]),g[0])):
        assert sum(a*b for a,b in zip(v,PROBE))!=0
        name=f'g{len(ops)}';terms={}
        for u,s in children:
            key=names[u];coef=ori(v)*s*ori(u)
            terms[key]=terms.get(key,0)+coef
        assert any(a>0 for a in terms.values())
        assert all(a in (-1,1,2) for a in terms.values())
        ops.append((name,list(terms.items())));names[v]=name
    for r in A:
        v,s=canonical(r); assert ori(v)==s
    c=circuit(f'A49_CF_{index:05d}',4,ops,[names[canonical(r)[0]] for r in A])
    assert c.matrix()==A
    return c


def canceling_fifteen():
    return circuit('A49_canceling_15',4,[
        ('s12',[('x1',1),('x2',1)]),
        ('s123',[('s12',1),('x3',1)]),
        ('w',[('x0',1),('s123',-1)]),
        ('d0',[('x0',2)]),('t0',[('d0',1),('x0',1)]),
        ('y0',[('t0',1),('w',-1)]),
        ('d2',[('x2',2)]),('t1',[('x1',1),('d2',-1)]),
        ('y1',[('t1',1),('w',-1)]),
        ('d3',[('x3',2)]),('t2',[('x2',1),('d3',-1)]),
        ('y2',[('t2',1),('w',-1)]),
        ('d1',[('x1',2)]),('t3',[('x3',1),('d1',-1)]),
        ('y3',[('t3',1),('w',-1)]),
    ],['y0','y1','y2','y3'])


@lru_cache(None)
def coefficients(k):
    return tuple(a for a in product(range(-4,5),repeat=k) if l1(a)<=4)


def representations(target,bases):
    out=[]
    for a in coefficients(len(bases)):
        u=tuple(target[j]-sum(t*b[j] for t,b in zip(a,bases)) for j in range(4))
        if l1(a)+l1(u)<=4:out.append((a,u))
    return out


def five_reduction_certificate():
    # First useful non-output normalization must have an integer coefficient
    # vector of l1 <=4. Vectors 0 and +/- e_i cannot enable the first output.
    universe=[v for v in product(range(-4,5),repeat=4)
              if 2<=l1(v)<=4 and next(x for x in v if x)>0]
    assert len(universe)==156
    records=[];successful=[]
    for w in universe:
        reached={0}
        for mask in range(16):
            if mask not in reached:continue
            bases=[w]+[A[i] for i in range(4) if mask>>i&1]
            for i in range(4):
                if mask>>i&1:continue
                if representations(A[i],bases): reached.add(mask|1<<i)
        records.append(dict(w=w,reachable_masks=sorted(reached)))
        if 15 in reached:successful.append(w)
    assert successful==[(1,-1,-1,-1),(1,-1,1,1),(1,1,-1,1),(1,1,1,-1)]
    unique=[]
    for w in successful:
        vectors=[];rowdata=[]
        for i,r in enumerate(A):
            bases=[w]+[A[j] for j in range(4) if j!=i]
            rep=representations(r,bases)
            assert len(rep)==1
            a,u=rep[0]
            assert a[1:]==(0,0,0) and abs(a[0])==1 and l1(u)==3
            v=u+(a[0],); vectors.append(v)
            rowdata.append(dict(output=i,w_coefficient=a[0],input_coefficients=u))
        # Gate sharing must correspond to a common signed subform with load>=2.
        # Pre-normalization w gates contain only x-atoms, not the reset w atom.
        sets=[{canonical(v)[0] for v in subforms(t) if l1(v)>=2} for t in vectors]
        pre={canonical(v+(0,))[0] for v in subforms(w) if 2<=l1(v)<4}
        assert all(not (pre&s) for s in sets)
        assert all(not (sets[i]&sets[j]) for i in range(4) for j in range(i))
        unique.append(dict(w=w,outputs=rowdata,proper_shared_forms=0))
    c=canceling_fifteen();assert c.matrix()==A
    z=optimize(c,4); b=brute(c,4)
    assert z['reductions']==b['reductions']==5
    assert sum(not v.output for v in c.nodes)==15
    # Independent exact integer polynomial checks, no symbolic dependency.
    B=[list(A[i]) for i in [1,0,3,2]]
    mm=lambda X,Y:[[sum(a*b for a,b in zip(r,col)) for col in zip(*Y)] for r in X]
    B2=mm(B,B);B4=mm(B2,B2)
    assert all(B4[i][j]-10*B2[i][j]+49*(i==j)==0 for i in range(4) for j in range(4))
    assert B2[0][1]==2
    return dict(universe_size=len(universe),all_first_forms=records,successful_forms=unique,
                max_reachable_histogram=dict(Counter(max(m.bit_count() for m in r['reachable_masks']) for r in records)),
                minimum_arithmetic_operations_for_five_normalizations_at_H4=15,
                witness_15=dict(nodes=[dict(name=v.name,terms=v.terms,output=v.output) for v in c.nodes],outputs=c.outputs,schedule=z),
                B=B,B_squared=B2,polynomial_B=[1,0,-10,0,49],polynomial_B_squared=[1,-10,49])


def main():
    cert=five_reduction_certificate()
    (ROOT/'data/five_reduction_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    print('Universal H4 five-normalization certificate: 156 first forms, 4 successful, minimum 15 operations.',flush=True)
    allg,stages=enumerate_catalog();assert len(allg)==13851
    hist=Counter(map(len,allg));assert hist=={12:243,13:13608}
    results=[];counts=Counter();best={};crosschecks=0
    path=ROOT/'data/cancellation_free_catalog.jsonl'
    with path.open('w') as f:
        for i,g in enumerate(allg):
            c=materialize(g,i);n=len(g);row=dict(id=i,operations=n,costs={})
            witnesses={}
            for H in [3,4,5]:
                r=optimize(c,H);assert r is not None
                row['costs'][str(H)]=r['reductions']
                counts[(n,H,r['reductions'])]+=1
                key=(n,H)
                if key not in best or r['reductions']<best[key]['cost']:
                    best[key]=dict(id=i,cost=r['reductions'],schedule=r,
                                   nodes=[dict(name=v.name,terms=v.terms,output=v.output) for v in c.nodes],outputs=c.outputs)
                if i%1000==0 and H==4:
                    assert brute(c,H)['reductions']==r['reductions'];crosschecks+=1
                witnesses[str(H)]=dict(cost=r['reductions'],selected=r['selected'])
            f.write(json.dumps(dict(id=i,gates=sorted(g),results=witnesses),separators=(',',':'))+'\n')
            results.append(row)
            if (i+1)%2000==0:print('Scheduled',i+1,'catalog DAGs.',flush=True)
    out=dict(grammar='unique-vector cancellation-free binary signed addition DAGs modulo wire sign; <=13 gates',
             per_output_trees=[60]*4,stage_counts=stages,catalog_size=len(allg),operation_histogram=dict(hist),
             cost_histogram=[dict(operations=n,H=H,cost=cost,count=c) for (n,H,cost),c in sorted(counts.items())],
             best=[dict(operations=n,H=H,**v) for (n,H),v in sorted(best.items())],
             dp_instances=3*len(allg),brute_crosschecks=crosschecks,
             catalog_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),results=results)
    (ROOT/'data/catalog_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Catalog complete:',hist,'DP cases',out['dp_instances'],'brute',crosschecks,flush=True)
    print('Minima',[(n,H,v['cost']) for (n,H),v in sorted(best.items())],flush=True)

if __name__=='__main__':main()
