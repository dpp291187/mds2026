"""Independent set-sum check of the 156-form normalization certificate.

Unlike catalog.py, this checker never enumerates residual coefficients or solves
for u. It builds the radius-four sumset of all signed canonical source vectors.
"""
from pathlib import Path
from itertools import product
import json
from reduction_model import signed_m4, evaluate, assess
from catalog import canceling_fifteen
ROOT=Path(__file__).resolve().parents[1]
A=signed_m4().matrix();ZERO=(0,0,0,0)
E=[tuple(int(i==j) for i in range(4)) for j in range(4)]

def radius_four(bases):
    steps=set(bases)|{tuple(-a for a in b) for b in bases}
    reached={ZERO};front={ZERO}
    for _ in range(4):
        nxt={tuple(a+b for a,b in zip(x,y)) for x in front for y in steps}-reached
        reached.update(nxt);front=nxt
    return reached

def main():
    d=json.loads((ROOT/'data/five_reduction_certificate.json').read_text())
    count=0
    for record in d['all_first_forms']:
        w=tuple(record['w']); masks={0}
        for mask in range(16):
            if mask not in masks:continue
            points=radius_four(E+[w]+[A[i] for i in range(4) if mask>>i&1])
            count+=1
            for i in range(4):
                if not(mask>>i&1) and A[i] in points:masks.add(mask|1<<i)
        assert sorted(masks)==record['reachable_masks']
    c=canceling_fifteen();selected=d['witness_15']['schedule']['selected']
    schedule=assess(c,4,selected);assert schedule['reductions']==5
    q=6;semantic=0
    for x in product(range(-q,q+1),repeat=4):
        y,values=evaluate(c,13,selected,x,True)
        expected=[sum(a*b for a,b in zip(r,x))%13 for r in A]
        assert [a%13 for a in y]==expected
        for k,v in enumerate(c.nodes,c.inputs):
            raw=sum(a*values[j] for j,a in v.terms)
            assert abs(raw)<=q*schedule['details'][k-c.inputs]['raw_load']
            assert abs(values[k])<=q*schedule['loads'][k]
        semantic+=1
    out=dict(independent_sumset_instances=count,first_forms=156,
             matching_reachability=True,centered_semantic_inputs=semantic,
             witness_operations=15,witness_reductions=5,capacity=4)
    (ROOT/'data/independent_synthesis_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))
if __name__=='__main__':main()
