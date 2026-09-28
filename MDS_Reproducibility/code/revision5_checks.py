"""Exact revision-5 observability and output-wiring certificates.

Python standard library only. The 24 row permutations of each supplied matrix
are enumerated exhaustively. The chosen permutation is one explicit design
choice, not a claim of exhaustive circuit synthesis or global arithmetic
optimality. Range attainment is checked only for the selected schedules over
centered F13 inputs. Arithmetic depth excludes normalization and routing delay.
"""
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt
from pathlib import Path
import json

from classify import A49, A7, A17, B, observability_matrix, characteristic_polynomial
from class_extension import a7_eleven
from reduction_model import det_integer, minor_certificate, optimize, evaluate

ROOT = Path(__file__).resolve().parents[1]
ROW_ORDER = (0, 3, 2, 1)
OBSERVABLE_A7 = tuple(A7[i] for i in ROW_ORDER)


def observable_a7_eleven():
    """Return the existing A7 DAG with output tuple (y0,y3,y2,y1).

    Dataclass replacement changes only the name and order of output IDs. The
    original arithmetic nodes, output sinks and all dependencies remain intact.
    """
    original = a7_eleven()
    result = replace(original, name='Observable_A7_eleven',
                     outputs=tuple(original.outputs[i] for i in ROW_ORDER))
    result.validate()
    assert result.nodes is original.nodes
    assert result.matrix() == OBSERVABLE_A7
    return result


def multiply(A, B):
    return tuple(tuple(sum(a*b for a, b in zip(row, col))
                       for col in zip(*B)) for row in A)


def independent_observability(A, j):
    """Construct matrix powers first, independently of the Krylov-row helper."""
    n = len(A)
    power = tuple(tuple(int(i == k) for k in range(n)) for i in range(n))
    rows = []
    for _ in range(n):
        rows.append(power[j])
        power = multiply(power, A)
    return tuple(rows)


def leibniz_determinant(A):
    n = len(A)
    total = 0
    for pi in permutations(range(n)):
        term = (-1)**sum(pi[i] > pi[j] for i in range(n)
                         for j in range(i+1, n))
        for i, j in enumerate(pi):
            term *= A[i][j]
        total += term
    return total


def prime_divisors(n):
    n = abs(n)
    factors = []
    divisor = 2
    while divisor*divisor <= n:
        if n % divisor == 0:
            factors.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        factors.append(n)
    return factors


def rational_rank(A):
    a = [[Fraction(v) for v in row] for row in A]
    rank = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        leading = a[rank][j]
        a[rank] = [v/leading for v in a[rank]]
        for i in range(len(a)):
            if i != rank:
                leading = a[i][j]
                a[i] = [v-leading*w for v, w in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def observation_record(A):
    determinants = []
    for j in range(4):
        independent = independent_observability(A, j)
        supplied = tuple(observability_matrix(A, j))
        assert independent == supplied
        d = leibniz_determinant(independent)
        assert d == det_integer(supplied)
        determinants.append(d)
    # Newton identities check the existing polynomial computation independently.
    powers = [tuple(tuple(int(i == j) for j in range(4)) for i in range(4))]
    for _ in range(4):
        powers.append(multiply(powers[-1], A))
    traces = [sum(P[i][i] for i in range(4)) for P in powers]
    coefficients = [Fraction(1)]
    for k in range(1, 5):
        coefficients.append(-sum(coefficients[k-i]*traces[i]
                                 for i in range(1, k+1))/k)
    assert all(c.denominator == 1 for c in coefficients)
    coefficients = list(map(int, coefficients))
    assert coefficients == list(reversed(characteristic_polynomial(A)))
    return dict(matrix=A, determinant=leibniz_determinant(A),
                characteristic_coefficients_descending=coefficients,
                observation_determinants=determinants,
                all_coordinate_observable_over_Q=all(determinants),
                singular_coordinates=[j for j, d in enumerate(determinants) if not d],
                prime_divisors_of_nonzero_observation_determinants=sorted(
                    {p for d in determinants if d for p in prime_divisors(d)}))


def check_observability():
    bases = {name: observation_record(A) for name, A in
             [('A49', A49), ('A7', A7), ('A17', A17), ('B', B),
              ('observable_A7', OBSERVABLE_A7)]}
    expected = dict(A49=[0, 0, 0, 0], A7=[-24, -24, -24, 0],
                    A17=[24, 24, 24, 24], B=[140, 140, 140, 140],
                    observable_A7=[-36, -4, -4, 28])
    for name, ds in expected.items():
        assert bases[name]['observation_determinants'] == ds
    census = {}
    for name, A, expected_count, expected_uniform in [
            ('A7', A7, 22, 4), ('A17', A17, 24, 7)]:
        rows = []
        for pi in permutations(range(4)):
            row = observation_record(tuple(A[i] for i in pi))
            row['row_order'] = pi
            row['no_observation_prime_above_seven'] = (
                row['all_coordinate_observable_over_Q'] and
                all(p <= 7 for p in row[
                    'prime_divisors_of_nonzero_observation_determinants']))
            rows.append(row)
        count = sum(r['all_coordinate_observable_over_Q'] for r in rows)
        uniform = sum(r['no_observation_prime_above_seven'] for r in rows)
        assert (count, uniform) == (expected_count, expected_uniform)
        census[name] = dict(examined=len(rows), observable_over_Q=count,
                            observable_without_new_prime_above_seven=uniform,
                            observable_coordinate_histogram=dict(Counter(
                                4-len(r['singular_coordinates']) for r in rows)),
                            records=rows)
    return dict(base_matrices=bases, permutation_census=census,
                observability_matrix_method_comparisons=212,
                observation_determinant_method_comparisons=212,
                characteristic_polynomial_method_comparisons=53)


def check_invariant_plane():
    U = ((1, 0), (0, 1), (1, 1), (0, 0))
    C = ((3, 2), (-2, 1))
    O = independent_observability(A7, 3)
    assert multiply(A7, U) == multiply(U, C)
    assert multiply(O, U) == ((0, 0),)*4
    assert rational_rank(U) == 2 and rational_rank(O) == 2
    assert O == ((0, 0, 0, 1), (1, 1, -1, 2),
                 (4, 4, -4, 7), (15, 15, -15, 26))
    assert O[2] == tuple(4*x-y for x, y in zip(O[1], O[0]))
    assert O[3] == tuple(15*x-4*y for x, y in zip(O[1], O[0]))
    # The 2x2 minor of rows 0,1 and columns 3,0 equals 1. The integer
    # row relations show rank at most 2 in every characteristic, hence exactly 2.
    assert leibniz_determinant([[O[i][j] for j in (3, 0)]
                               for i in (0, 1)]) == 1
    assert leibniz_determinant(C) == 7
    return dict(basis_columns=U, restricted_matrix=C, restricted_determinant=7,
                observation_matrix=O, rank=2, rank_valid_over_every_field=True,
                kernel_parameterization='(a,b,a+b,0)',
                action_parameterization='(3a+2b,-2a+b,a+3b,0)',
                nonzero_persistence_condition='characteristic not equal to 7',
                admissible_MDS_prime_condition='p > 7')


def check_minors():
    certs = {}
    comparisons = 0
    for name, A in [('A7', A7), ('observable_A7', OBSERVABLE_A7)]:
        cert = minor_certificate(A)
        for m in cert:
            submatrix = [[A[i][j] for j in m['cols']] for i in m['rows']]
            assert leibniz_determinant(submatrix) == m['det']
            comparisons += 1
        assert len(cert) == 69 and all(m['det'] for m in cert)
        certs[name] = cert
    spectra = {name: {str(k): sorted(abs(m['det']) for m in cert
                                   if len(m['rows']) == k)
                       for k in range(1, 5)} for name, cert in certs.items()}
    assert spectra['A7'] == spectra['observable_A7']
    bad = sorted({p for m in certs['observable_A7']
                  for p in prime_divisors(m['det'])})
    row_norms = [sum(map(abs, row)) for row in OBSERVABLE_A7]
    assert row_norms == [5]*4 and bad == [2, 3, 5, 7]
    return dict(row_norms=row_norms, norm=5, bad_primes=bad,
                minor_method_comparisons=comparisons,
                equal_absolute_minor_multisets=True,
                absolute_minor_multisets=spectra['observable_A7'],
                selected_matrix_minor_certificate=certs['observable_A7'])


def independent_brute(c, H):
    """Enumerate every optional subset without the DP or assess routine."""
    internal = [i for i, node in enumerate(c.nodes, c.inputs) if not node.output]
    best = None
    feasible = 0
    for bits in range(1 << len(internal)):
        chosen = {i for k, i in enumerate(internal) if bits >> k & 1}
        load = [1]*c.inputs
        count = 0
        valid = True
        for i, node in enumerate(c.nodes, c.inputs):
            raw = sum(abs(a)*load[j] for j, a in node.terms)
            if raw > H:
                valid = False
                break
            normalize = raw > 1 if node.output else i in chosen
            load.append(1 if normalize else raw)
            count += int(normalize)
        if valid:
            feasible += 1
            if best is None or count < best:
                best = count
    return dict(H=H, enumerated_subsets=1 << len(internal),
                feasible_subsets=feasible, optimum=best)


def check_circuit():
    original = a7_eleven()
    selected = observable_a7_eleven()
    assert selected.nodes == original.nodes
    assert selected.outputs == tuple(original.outputs[j] for j in ROW_ORDER)
    depths = [0]*selected.inputs
    for node in selected.nodes:
        depths.append(max(depths[j] for j, _ in node.terms) + int(not node.output))
    operations = sum(not node.output for node in selected.nodes)
    depth = max(depths[j] for j in selected.outputs)
    assert (operations, depth) == (11, 4)
    sweep = []
    for H in range(2, 17):
        a = optimize(original, H)
        ad = optimize(original, H, dominance=True)
        b = optimize(selected, H)
        bd = optimize(selected, H, dominance=True)
        assert a == b and ad == bd
        assert b['reductions'] == bd['reductions']
        sweep.append(dict(H=H, normalizations=b['reductions'],
                          plain_schedule=b, dominance_schedule=bd,
                          all_original_results_identical=True))
    assert [r['normalizations'] for r in sweep] == [11, 7, 5] + [4]*12
    brute_results = []
    for H, expected in [(3, 7), (4, 5), (5, 4)]:
        b = independent_brute(selected, H)
        assert b['optimum'] == optimize(selected, H)['reductions'] == expected
        brute_results.append(b)
    return dict(name=selected.name, matrix=selected.matrix(),
                operations=operations, arithmetic_depth=depth,
                arithmetic_depth_excludes_normalization=True,
                row_order=ROW_ORDER, returned_output_names=['y0','y3','y2','y1'],
                original_output_ids=original.outputs, selected_output_ids=selected.outputs,
                arithmetic_nodes_and_sinks_unchanged=True,
                nodes=[dict(name=n.name, terms=n.terms, output_sink=n.output)
                       for n in selected.nodes],
                capacity_sweep=sweep, capacity_count=15, DP_runs=60,
                original_vs_permuted_full_result_equalities=30,
                plain_vs_dominance_cost_comparisons=15,
                independent_subset_checks=brute_results,
                independent_subsets_examined=sum(r['enumerated_subsets'] for r in brute_results),
                global_arithmetic_optimality_claim=False)


def check_centered_ranges():
    c = observable_a7_eleven()
    p = 13
    q = 6
    A = c.matrix()
    cases = []
    for H in [4, 5]:
        r = optimize(c, H)
        raw_peak = [0]*len(c.nodes)
        retained_peak = [0]*len(c.nodes)
        raw_witness = [None]*len(c.nodes)
        retained_witness = [None]*len(c.nodes)
        inputs = 0
        for x in product(range(-q, q+1), repeat=4):
            outputs, values = evaluate(c, p, r['selected'], x, centered=True)
            expected = [sum(a*b for a, b in zip(row, x)) % p for row in A]
            expected = [v-p if v > q else v for v in expected]
            assert outputs == expected
            for k, node in enumerate(c.nodes):
                raw = sum(a*values[j] for j, a in node.terms)
                retained = values[c.inputs+k]
                assert abs(raw) <= q*r['details'][k]['raw_load']
                assert abs(retained) <= q*r['details'][k]['post_load']
                if abs(raw) > raw_peak[k]:
                    raw_peak[k] = abs(raw)
                    raw_witness[k] = x
                if abs(retained) > retained_peak[k]:
                    retained_peak[k] = abs(retained)
                    retained_witness[k] = x
            inputs += 1
        assert inputs == 28561
        rows = []
        for k, node in enumerate(c.nodes):
            raw_bound = q*r['details'][k]['raw_load']
            retained_bound = q*r['details'][k]['post_load']
            assert raw_peak[k] == raw_bound and retained_peak[k] == retained_bound
            rows.append(dict(name=node.name, output_sink=node.output,
                             raw_certificate=raw_bound, exact_raw_peak=raw_peak[k],
                             raw_attaining_input=raw_witness[k],
                             retained_certificate=retained_bound,
                             exact_retained_peak=retained_peak[k],
                             retained_attaining_input=retained_witness[k]))
        cases.append(dict(H=H, prime=p, inputs=inputs,
                          normalizations=r['reductions'], selected_internal_nodes=r['selected'],
                          arithmetic_raw_bounds_attained=11,
                          all_raw_bounds_attained=len(c.nodes),
                          all_retained_bounds_attained=len(c.nodes),
                          node_bounds=rows))
        print('Centered F13 output/range checks:', H, inputs, flush=True)
    return dict(scope='Selected certificate-optimal schedules at H=4 and H=5 only',
                cases=cases, input_schedule_evaluations=57122,
                output_coordinate_equalities=228488,
                raw_bound_inequality_checks=856830,
                retained_bound_inequality_checks=856830,
                arithmetic_raw_bounds_attained=22, all_raw_bounds_attained=30,
                all_retained_bounds_attained=30)


def check_named_fields(minors, observation):
    saved = json.loads((ROOT/'data/revision4_checks.json').read_text())['fields']
    expected = {'BabyBear':2013265921, 'Mersenne31':2147483647,
                'Goldilocks':18446744069414584321, 'KoalaBear':2130706433}
    rows = []
    for field in saved:
        name, p = field['name'], field['p']
        assert expected[name] == p
        cert = field['primality_certificate']
        factors = {int(q): e for q, e in cert['factorization'].items()}
        product_value = 1
        for q, e in factors.items():
            assert q >= 2 and all(q % d for d in range(2, isqrt(q)+1))
            product_value *= q**e
        assert product_value == p-1
        a = cert['primitive_root']
        assert pow(a, p-1, p) == 1
        assert all(gcd(pow(a, (p-1)//q, p)-1, p) == 1 for q in factors)
        residues = [m['det'] % p for m in minors['selected_matrix_minor_certificate']]
        obs_residues = [d % p for d in observation[
            'base_matrices']['observable_A7']['observation_determinants']]
        assert all(residues) and all(obs_residues)
        a17_obs = [d % p for d in observation['base_matrices']['A17']['observation_determinants']]
        assert all(a17_obs)
        rows.append(dict(name=name, prime=p, reused_primality_certificate_reverified=True,
                         selected_matrix_MDS=True, nonzero_minors=len(residues),
                         selected_observation_determinants_mod_p=obs_residues,
                         selected_all_coordinate_observable=True,
                         A17_all_coordinate_observable=True,
                         A17_MDS=p not in {2,3,5,7,17}))
    assert set(expected) == {r['name'] for r in rows}
    return dict(fields=rows, reverified_primality_certificates=4,
                selected_nonzero_minor_tests=276,
                selected_nonzero_observation_tests=16,
                A17_nonzero_observation_tests=16)


def main():
    observations = check_observability()
    plane = check_invariant_plane()
    minors = check_minors()
    witness = check_circuit()
    ranges = check_centered_ranges()
    fields = check_named_fields(minors, observations)
    output = dict(
        methodology=dict(
            row_permutation_convention='A_pi[i][j] = A[pi[i]][j], zero-based indexing',
            finite_census='All 24 row permutations of each A7 and A17; exhaustive',
            candidate_choice='Explicit swap of output coordinates 1 and 3; no global circuit search claim',
            output_wiring_cost='Free output relabeling in the stated DAG model',
            security_scope='Only invariant subspaces contained in a single coordinate hyperplane',
            recommendation_scope='H in {4,5}; compare attained operation counts and arithmetic depths',
            global_arithmetic_minimum_proved=False),
        observability=observations, A7_invariant_plane=plane,
        selected_matrix=minors, circuit=witness, centered_F13=ranges, named_fields=fields)
    target = ROOT/'data/revision5_checks.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    print('Revision 5 checks complete:', ranges['input_schedule_evaluations'],
          'input/schedule evaluations;', observations[
              'observation_determinant_method_comparisons'],
          'independent observation determinants;', witness['DP_runs'], 'DP runs.', flush=True)


if __name__ == '__main__':
    main()
