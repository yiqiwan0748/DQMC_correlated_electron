"""Fixed QuSpin dimer exercise; numerical calls occur ONLY in a root-admitted PC unit.

Importing this module performs no physical computation and imports no scientific
package. ``run`` receives actual source-qualified exports from the new worker.
The analytic reference is evaluated inside that future unit, never at import.
"""
import math
from fractions import Fraction


MODEL = {'sites': 2, 't': 1, 'U': 4, 'mu': 0, 'Nup': 1, 'Ndown': 1,
         'dimension': 4, 'boundary': 'open', 'interaction': 'uncentered_nup_ndown',
         'energy_offset': 0}
ENERGY_LIMIT = Fraction(1, 10**11)
HERMITIAN_LIMIT = Fraction(1, 10**12)


def need(value, message):
    if not value:
        raise ValueError(message)


def fraction_record(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator)}


def rational(value):
    need(type(value) is float and math.isfinite(value), 'Finite stored binary64 scalar required')
    numerator, denominator = value.as_integer_ratio()
    return Fraction(numerator, denominator)


def interval_distance(value, interval):
    """Conservative distance to every value of a closed rational interval."""
    low, high = interval
    need(type(value) is Fraction and type(low) is type(high) is Fraction and low <= high,
         'Ordered exact rational reference interval required')
    return max(abs(value-low), abs(value-high))


def analytic_intervals():
    # No float sqrt, midpoint fit, observed spectrum, or target matrix is input.
    q = 10**80
    k = math.isqrt(32*q*q)
    need(k*k <= 32*q*q < (k+1)*(k+1), 'Integer square-root enclosure failed')
    intervals = [(Fraction(4*q-k-1, 2*q), Fraction(4*q-k, 2*q)),
                 (Fraction(0), Fraction(0)), (Fraction(4), Fraction(4)),
                 (Fraction(4*q+k, 2*q), Fraction(4*q+k+1, 2*q))]
    return intervals, {'q': str(q), 'floor_sqrt_32_q_squared': str(k),
                       'integer_enclosure_checked': True,
                       'intervals': [[fraction_record(x) for x in pair] for pair in intervals],
                       'scope': 'Analytic spectrum of the fixed uncentered two-site Hubbard sector only.'}


def stored_array(array, np, shape):
    need(type(array) is np.ndarray and tuple(array.shape) == tuple(shape),
         'Actual NumPy array has unexpected exact type or shape')
    need(array.dtype.kind == 'f' and array.dtype.itemsize == 8 and str(array.dtype) == 'float64',
         'Finite real native binary64 arrays required; no dtype cast or complex discard')
    if len(shape) == 1:
        data = [float(array[i]) for i in range(shape[0])]
    else:
        need(len(shape) == 2, 'Only registered vector and matrix outputs supported')
        data = [[float(array[i, j]) for j in range(shape[1])] for i in range(shape[0])]
    flat = data if len(shape) == 1 else [x for row in data for x in row]
    need(all(math.isfinite(x) for x in flat), 'Nonfinite numerical output')
    return data


def physical_masks(label):
    need(type(label) is int and 0 <= label < 16, 'Registered four-orbital integer label required')
    reverse_two = lambda value: ((value & 1) << 1) | ((value >> 1) & 1)
    return reverse_two(label >> 2), reverse_two(label & 3)


def point_diagnostics(H, E, V):
    # Exact arithmetic of the STORED binary64 values. This is not an enclosure
    # for a true normalized eigenvector or a ground-state/gap certificate.
    h = [[rational(x) for x in row] for row in H]
    e = [rational(x) for x in E]
    v = [[rational(x) for x in row] for row in V]
    residual = [[sum(h[i][k]*v[k][j] for k in range(4))-v[i][j]*e[j]
                 for j in range(4)] for i in range(4)]
    gram = [[sum(v[k][i]*v[k][j] for k in range(4))-Fraction(i == j)
             for j in range(4)] for i in range(4)]
    return {'stored_dyadic_residual_matrix': [[fraction_record(x) for x in row] for row in residual],
            'stored_dyadic_Gram_defect_matrix': [[fraction_record(x) for x in row] for row in gram],
            'residual_max_abs': fraction_record(max(abs(x) for row in residual for x in row)),
            'Gram_max_abs': fraction_record(max(abs(x) for row in gram for x in row)),
            'qualification': 'POINT_DIAGNOSTICS_NOT_EIGENVECTOR_OR_GLOBAL_GAP_CERTIFICATES'}


def run(np, spinful_class, hamiltonian_class, invoke, save, check_time):
    """No package is selected here; the admitted worker supplies pinned exports.

    ``invoke`` records and checks actual exported class or bound-method entry
    identities before each registered call. The caller checks fresh complete
    maps/modules/providers/API witnesses before AND after this function.
    """
    check_time()
    basis = invoke('basis_constructor', spinful_class, (2,),
                   {'Nf': (1, 1), 'double_occupancy': True}, spinful_class)
    need(type(basis) is spinful_class and int(basis.Ns) == 4, 'Exact fixed-spin dimension4 basis required')
    states = basis.states
    need(type(states) is np.ndarray and tuple(states.shape) == (4,) and
         states.dtype.kind in ('u', 'i') and 1 <= states.dtype.itemsize <= 8,
         'Actual stored labels must be a complete finite native integer array')
    labels = [int(value) for value in states]
    masks = [physical_masks(value) for value in labels]
    need(len(labels) == len(set(labels)) == 4 and set(masks) == {(1,1), (1,2), (2,1), (2,2)},
         'Complete four distinct Nup1/Ndown1 physical mask pairs required')
    save('basis_labels', {'labels': labels, 'physical_masks': [list(x) for x in masks],
                          'basis_dimension': int(basis.Ns), 'stored_label_dtype': str(states.dtype), 'model': MODEL})
    flags = {'basis': basis, 'dtype': np.float64, 'check_herm': True,
             'check_symm': False, 'check_pcon': False}
    numbers = []
    for species, symbol, mask_index in [('up', 'n|', 0), ('down', '|n', 1)]:
        for site in range(2):
            label = 'number_'+species+'_'+str(site)
            operator = invoke(label+'_constructor', hamiltonian_class,
                              ([[symbol, [[1.0, site]]]], []), flags, hamiltonian_class)
            need(type(operator) is hamiltonian_class, 'Exact actual number-operator class required')
            raw = invoke(label+'_toarray', operator.toarray, (), {}, hamiltonian_class, operator)
            matrix = stored_array(raw, np, (4,4))
            save(label, {'matrix': matrix, 'species': species, 'site': site})
            for i in range(4):
                for j in range(4):
                    expected = float((masks[i][mask_index] >> site) & 1) if i == j else 0.0
                    need(rational(matrix[i][j]) == rational(expected),
                         'Actual number operator disagrees with high-bit to little-site mask decoding')
            numbers.append({'species': species, 'site': site, 'matrix': matrix})
            check_time()
    # Coefficients follow the independently constructed uncentered-chain oracle:
    # direct Kij=-t in BOTH spin strings, and U*n_up*n_down (never z|z).
    static = [['+-|', [[-1.0, 0, 1], [-1.0, 1, 0]]],
              ['|+-', [[-1.0, 0, 1], [-1.0, 1, 0]]],
              ['n|n', [[4.0, 0, 0], [4.0, 1, 1]]]]
    Hobject = invoke('Hamiltonian_constructor', hamiltonian_class, (static, []), flags, hamiltonian_class)
    need(type(Hobject) is hamiltonian_class, 'Exact constructed Hamiltonian type required')
    Hraw = invoke('Hamiltonian_toarray', Hobject.toarray, (), {}, hamiltonian_class, Hobject)
    H = stored_array(Hraw, np, (4,4))
    save('Hamiltonian', {'matrix': H, 'static_terms': static, 'model': MODEL})
    hermitian = max(abs(rational(H[i][j])-rational(H[j][i])) for i in range(4) for j in range(4))
    need(hermitian <= HERMITIAN_LIMIT, 'Registered Hermitian point-matrix gate failed')
    output = invoke('Hamiltonian_eigh', Hobject.eigh, (), {}, hamiltonian_class, Hobject)
    need(type(output) in (tuple, list) and len(output) == 2, 'Complete actual eigh outputs required')
    E, V = stored_array(output[0], np, (4,)), stored_array(output[1], np, (4,4))
    save('complete_eigensystem', {'energies': E, 'vectors_by_column': V})
    need(all(rational(E[i]) <= rational(E[i+1]) for i in range(3)), 'Actual complete levels are not already sorted')
    intervals, analytic = analytic_intervals()
    distances = [interval_distance(rational(value), pair) for value, pair in zip(E, intervals)]
    analytic_gates = [{'index': i, 'exact_distance_upper': fraction_record(value),
                       'limit': fraction_record(ENERGY_LIMIT), 'pass': value <= ENERGY_LIMIT}
                      for i, value in enumerate(distances)]
    diagnostics = point_diagnostics(H, E, V)
    payload = {'schema': 'PC_FIXED_DIMER_NUMERICAL_OUTPUT_V1', 'model': MODEL,
               'basis_labels': labels, 'physical_masks': [list(x) for x in masks],
               'actual_number_operators': numbers, 'H': H, 'E': E, 'V_by_column': V,
               'analytic_reference': analytic, 'analytic_energy_gates': analytic_gates,
               'Hermitian_max_abs': fraction_record(hermitian),
               'Hermitian_limit': fraction_record(HERMITIAN_LIMIT), 'diagnostics': diagnostics,
               'scientific_credit': 0, 'X08_status': 'NOT_RUN',
               'eigenvector_global_gap_ground_state_qualification': 'NOT_CERTIFIED',
               'scope': 'Fixed QuSpin dimer functional pilot only. No native solver or X08 convergence claim.'}
    save('science_result', payload)
    need(all(item['pass'] for item in analytic_gates), 'One or more frozen analytic spectrum gates failed')
    check_time()
    return payload
