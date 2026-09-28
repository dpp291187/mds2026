"""All-class witnesses, exact finite ranges, field tests and block composition.

Python standard library only. No heuristic circuit-search result is used as an
arithmetic lower bound. All semantic ranges below are exact for the specified
prime, circuit and selected schedule, not for all possible schedules.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import gcd, isqrt
from pathlib import Path
import json
from classify import A49,A7,A17,B,REPS
from class_extension import a7_eleven,a17_eleven
from catalog import canceling_fifteen
from reduction_model import hl_m4,p3_m4,signed_m4,optimize,brute,assess,evaluate
ROOT=Path(__file__).resolve().parents[1]
ZERO=(0,0,0,0)
E=[tuple(int(i==j) for i in range(4)) for j in range(4)]


def radius4(bases):
    # An independent multiset sumset implementation, with no residual solving.
    steps=sorted(set(bases)|{tuple(-a for a in b) for b in bases})
    out={ZERO}
    for count in range(1,5):
        for terms in combinations_with_replacement(steps,count):
            out.add(tuple(map(sum,zip(*terms))))
    return out


def check_reachability():
    data=json.loads((ROOT/'data/all_class_reachability.json').read_text())
    out={}
    for label,A in REPS.items():
        count=0;hist=Counter()
        for row in data[str(label)]['all_first_forms']:
            reached={0};w=tuple(row['w'])
            for mask in range(16):
                if mask not in reached:continue
                points=radius4(E+[w]+[A[i] for i in range(4) if mask>>i&1]);count+=1
                for i in range(4):
                    if not(mask>>i&1) and A[i] in points:reached.add(mask|1<<i)
            assert sorted(reached)==row['reachable_masks']
            hist[max(m.bit_count() for m in reached)]+=1
        out[str(label)]=dict(independent_sumsets=count,first_forms=156,histogram=dict(hist))
        print('Independent reachability',label,count,dict(hist),flush=True)
    return out


def exact_ranges(c,H,p=13):
    best=optimize(c,H);assert best
    q=(p-1)//2;selected=best['selected'];A=c.matrix()
    raw_peak=[0]*len(c.nodes);post_peak=[0]*len(c.nodes)
    attaining_raw=[None]*len(c.nodes);attaining_post=[None]*len(c.nodes)
    inputs=0
    for x in product(range(-q,q+1),repeat=c.inputs):
        y,values=evaluate(c,p,selected,x,True)
        assert [a%p for a in y]==[sum(a*b for a,b in zip(row,x))%p for row in A]
        for k,v in enumerate(c.nodes):
            raw=sum(a*values[j] for j,a in v.terms);post=values[c.inputs+k]
            assert abs(raw)<=q*best['details'][k]['raw_load']
            assert abs(post)<=q*best['details'][k]['post_load']
            if abs(raw)>raw_peak[k]:raw_peak[k]=abs(raw);attaining_raw[k]=x
            if abs(post)>post_peak[k]:post_peak[k]=abs(post);attaining_post[k]=x
        inputs+=1
    rows=[]
    for k,v in enumerate(c.nodes):
        rb=q*best['details'][k]['raw_load'];pb=q*best['details'][k]['post_load']
        rows.append(dict(node=v.name,output_sink=v.output,raw_certificate=rb,exact_raw_peak=raw_peak[k],raw_attaining_input=attaining_raw[k],post_certificate=pb,exact_post_peak=post_peak[k],post_attaining_input=attaining_post[k]))
    ar=[r for r in rows if not r['output_sink']]
    tight=sum(r['raw_certificate']==r['exact_raw_peak'] for r in ar)
    ratio=min(Fraction(r['exact_raw_peak'],r['raw_certificate']) for r in ar)
    report=dict(circuit=c.name,H=H,p=p,inputs=inputs,operations=len(ar),normalizations=best['reductions'],selected=selected,arithmetic_raw_tight=tight,minimum_raw_utilization=str(ratio),maximum_raw_certificate=max(r['raw_certificate'] for r in ar),maximum_exact_raw_peak=max(r['exact_raw_peak'] for r in ar),nodes=rows)
    print('Ranges',c.name,'H',H,'cost',best['reductions'],'tight',f'{tight}/{len(ar)}','min',str(ratio),flush=True)
    return report


def factor(n):
    out={};d=2
    while d*d<=n:
        while n%d==0:out[d]=out.get(d,0)+1;n//=d
        d=3 if d==2 else d+2
    if n>1:out[n]=out.get(n,0)+1
    return out


def primality_certificate(p):
    fac=factor(p-1)
    # Full factorization + an element of order p-1 gives Lucas's prime test.
    for q in fac:assert factor(q)=={q:1}
    for a in range(2,10000):
        if pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fac):
            return dict(factorization=fac,primitive_root=a)
    raise AssertionError('No primality certificate found')


def square_root(n,p):
    n%=p
    if pow(n,(p-1)//2,p)!=1:return None
    if p%4==3:return pow(n,(p+1)//4,p)
    q=p-1;s=0
    while q%2==0:q//=2;s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1:z+=1
    m=s;c=pow(z,q,p);t=pow(n,q,p);r=pow(n,(q+1)//2,p)
    while t!=1:
        i=1;v=t*t%p
        while v!=1:v=v*v%p;i+=1
        b=pow(c,1<<(m-i-1),p)
        r=r*b%p;t=t*b*b%p;c=b*b%p;m=i
    assert r*r%p==n
    return min(r,p-r)


def rank_mod(M,p):
    a=[[v%p for v in row] for row in M];r=0
    for j in range(len(a[0])):
        i=next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None:continue
        a[r],a[i]=a[i],a[r];inv=pow(a[r][j],-1,p)
        a[r]=[v*inv%p for v in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][j];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
        r+=1
    return r


def field_checks():
    primes=[('BabyBear',15*2**27+1),('Mersenne31',2**31-1),('Goldilocks',2**64-2**32+1),('KoalaBear',2**31-2**24+1)]
    B2=[[sum(B[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    out=[]
    for name,p in primes:
        cert=primality_certificate(p);leg=pow(p-6,(p-1)//2,p)
        chi=1 if leg==1 else -1;assert leg in (1,p-1)
        assert (chi==1)==(p%24 in (1,5,7,11))
        root=square_root(-6,p);eig=[]
        if root is not None:
            assert (root*root+6)%p==0
            for value in [(5+2*root)%p,(5-2*root)%p]:
                rank=rank_mod([[B2[i][j]-value*(i==j) for j in range(4)] for i in range(4)],p)
                assert rank==2;eig.append(dict(value=value,eigenspace_dimension=4-rank))
        out.append(dict(name=name,p=p,residue_mod24=p%24,minus_six_character=chi,square_root=root,eigenvalues=eig,primality_certificate=cert))
    return out


def block_checks(fields):
    out=[]
    matrices=dict(REPS);matrices['P3']=p3_m4().matrix()
    for k in range(2,7):
        for name,A in matrices.items():
            C=[[A[i%4][j%4]*(1+(i//4==j//4)) for j in range(4*k)] for i in range(4*k)]
            L=max(sum(map(abs,row)) for row in A)
            assert max(sum(map(abs,row)) for row in C)==(k+1)*L
            assert C[0][0]*C[1][4]-C[0][4]*C[1][0]==0
        widths=[]
        for field in fields:
            p=field['p'];q=(p-1)//2
            w5=1+(q*(k+1)*5).bit_length();w7=1+(q*(k+1)*7).bit_length()
            widths.append(dict(field=field['name'],p=p,signed_block_width=w5,positive_block_width=w7,saving=w7-w5))
        out.append(dict(blocks=k,dimension=4*k,signed_norm=(k+1)*5,positive_norm=(k+1)*7,widths=widths))
    return out


def main():
    reach=check_reachability()
    witnesses=[];checks=0
    for c,A in [(a7_eleven(),A7),(a17_eleven(),A17)]:
        assert c.matrix()==A and sum(not v.output for v in c.nodes)==11
        schedules={}
        for H,want in [(3,None),(4,5),(5,4)]:
            r=optimize(c,H);b=brute(c,H);assert r['reductions']==b['reductions'];checks+=1
            if want is not None:assert r['reductions']==want
            schedules[str(H)]=r
        witnesses.append(dict(circuit=c.name,matrix=A,operations=11,nodes=[dict(name=v.name,terms=v.terms,output=v.output) for v in c.nodes],outputs=c.outputs,schedules=schedules))
    # Keep original curves; add the explicit new witnesses to the capacity sweep.
    grid=[]
    for c in [a7_eleven(),a17_eleven()]:
        for H in range(2,17):
            r=optimize(c,H);d=optimize(c,H,True);assert r['reductions']==d['reductions']
            grid.append(dict(circuit=c.name,H=H,cost=r['reductions'],selected=r['selected']))
    (ROOT/'data/class_witnesses.json').write_text(json.dumps(dict(witnesses=witnesses,subset_comparisons=checks,grid=grid),indent=2)+'\n')
    ranges=[]
    for c in [hl_m4(),p3_m4(),signed_m4(),a7_eleven(),a17_eleven()]:
        for H in [3,4,5]:ranges.append(exact_ranges(c,H))
    ranges.append(exact_ranges(canceling_fifteen(),4))
    (ROOT/'data/exact_ranges.json').write_text(json.dumps(dict(semantics='centered F13 inputs, selected certificate-optimal schedules, per-node absolute maxima',total_inputs=sum(r['inputs'] for r in ranges),cases=ranges),indent=2)+'\n')
    fields=field_checks();blocks=block_checks(fields)
    bound=[dict(n=n,lower_bound=2*n-(1+isqrt(8*n-7))//2) for n in range(2,9)]
    assert [r['lower_bound'] for r in bound]==[2,4,5,7,9,10,12]
    out=dict(reachability=reach,subset_comparisons=checks,range_cases=len(ranges),range_inputs=sum(r['inputs'] for r in ranges),fields=fields,block_compositions=blocks,general_bounds=bound)
    (ROOT/'data/revision4_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Revision 4 checks complete:',out['range_cases'],out['range_inputs'],flush=True)
    print('Field characters:',[(r['name'],r['minus_six_character']) for r in fields],flush=True)

if __name__=='__main__':main()
