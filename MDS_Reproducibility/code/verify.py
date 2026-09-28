"""Run exact certificate checks and write the reproducible results.json.

Run: python verify.py
Only Python's standard library is used. Runtime is recorded, not benchmarked.
"""
import json, random, time
from pathlib import Path
from itertools import combinations, permutations, product
from reduction_model import *
from classify import (classify, classify_positive, REPS, prime_divisors, B,
                      observability_matrix, characteristic_polynomial)

def det_leibniz(M):
    n=len(M); total=0
    for pi in permutations(range(n)):
        term=(-1)**sum(pi[i]>pi[j] for i in range(n) for j in range(i+1,n))
        for i,j in enumerate(pi):term*=M[i][j]
        total+=term
    return total

def main():
    start=time.perf_counter(); rng=random.Random(26092026)
    signed=classify(); positive=classify_positive()
    matrices=[r['matrix'] for r in signed['classes']+positive['classes']]
    checks=0
    for M in matrices:
        for m in minor_certificate(M):
            submatrix=[[M[i][j] for j in m['cols']] for i in m['rows']]
            assert m['det']==det_leibniz(submatrix);checks+=1
    for _ in range(400):
        n=rng.randrange(1,5)
        M=[[rng.randrange(-5,6) for j in range(n)] for i in range(n)]
        assert det_integer(M)==det_leibniz(M);checks+=1

    # Exact mathematical identities used in the exposition.
    A=REPS[49]
    assert all(sum(A[k][i]*A[k][j] for k in range(4))==7*int(i==j)
               for i in range(4) for j in range(4))
    assert all(sum(A[i][k]*A[k][j] for k in range(4))-4*A[i][j]+7*int(i==j)==0
               for i in range(4) for j in range(4))
    assert signed_m4().matrix()==A
    obs_dets=[det_integer(observability_matrix(B,j)) for j in range(4)]
    assert obs_dets==[140]*4
    assert all(det_leibniz(observability_matrix(B,j))==140 for j in range(4))
    assert all(det_integer(observability_matrix(A,j))==0 for j in range(4))

    dags=[hl_m4(),p3_m4(),signed_m4()]
    tables=[]; dp_brute=0; dominance_checks=0
    for c in dags:
        for H in [2,3,4,5,6,7,8,12,16]:
            r=optimize(c,H);b=brute(c,H)
            assert (None if r is None else r['reductions'])==(None if b is None else b['reductions'])
            dp_brute+=1
            assert r['reductions']==optimize(c,H,dominance=True)['reductions'];dominance_checks+=1
            tables.append(dict(circuit=c.name,H=H,**r))
    for sample in range(80):
        n=3; nodes=[]
        for i in range(7):
            ids=rng.sample(range(n+i),rng.choice((1,2)))
            terms=tuple((j,rng.choice((-2,-1,1,2))) for j in ids)
            nodes.append(Node(f'r{i}',terms))
        nodes.extend(Node(f'out{k}',((n+6-k,1),),True) for k in range(2))
        c=Circuit(f'random{sample}',n,tuple(nodes),(n+7,n+8));c.validate()
        for H in [2,3,4,5,6]:
            r=optimize(c,H);b=brute(c,H);d=optimize(c,H,dominance=True)
            cost=lambda z:None if z is None else z['reductions']
            assert cost(r)==cost(b)==cost(d)
            dp_brute+=1;dominance_checks+=1

    graphs=0
    for n in range(2,5):
        es=list(combinations(range(n),2))
        for mask in range(1<<len(es)):
            edges=[e for k,e in enumerate(es) if mask>>k&1]
            c=vertex_cover_circuit(n,edges)
            r=optimize(c,3)
            assert r['reductions']==len(edges)+min_vertex_cover(n,edges)
            graphs+=1

    # Exhaustive semantic and per-node raw-range checking, not sampling.
    p=13;q=p//2;semantic_vectors=0
    for c in dags:
        M=c.matrix()
        for H in [3,5]:
            r=optimize(c,H);pick=set(r['selected'])
            for x in product(range(-q,q+1),repeat=4):
                out,vals=evaluate(c,p,pick,x,centered=True)
                expected=[sum(a*b for a,b in zip(row,x))%p for row in M]
                assert [v%p for v in out]==expected
                for i,v in enumerate(c.nodes,c.inputs):
                    raw=sum(a*vals[j] for j,a in v.terms)
                    cert=r['details'][i-c.inputs]
                    assert abs(raw)<=cert['raw_load']*q<=H*q
                    assert abs(vals[i])<=cert['post_load']*q
                semantic_vectors+=1

    result=dict(model='unit-cost compositional absolute-load certificate; fixed DAG',
                signed_classification=signed,positive_classification=positive,
                exact_reduction_tables=tables,
                circuit_matrices={c.name:c.matrix() for c in dags},
                observable_norm_five=dict(matrix=B,observability_determinants=obs_dets,
                                          characteristic_polynomial=characteristic_polynomial(B),
                                          certificate=minor_certificate(B)),
                verification=dict(independent_determinants=checks,
                                  dp_vs_exhaustive_subsets=dp_brute,
                                  dominance_checks=dominance_checks,
                                  vertex_cover_graphs=graphs,
                                  exhaustive_centered_vectors=semantic_vectors,
                                  semantic_prime=p),
                elapsed_seconds=time.perf_counter()-start)
    target=Path(__file__).resolve().parents[1] / 'data' / 'results.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['verification'],indent=2))
    print('seconds',round(result['elapsed_seconds'],3))
    print('written',target)

if __name__=='__main__':main()
