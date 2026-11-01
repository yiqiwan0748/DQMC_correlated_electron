"""Exact full-only ordered P13 gate metadata schedule, derived from frozen call sites.
No model imports, files, scalar certificates, admission or scientific execution.
Gate names intentionally repeat for each of the4096 actual pipe entries per U.
"""
SOURCE_BINDINGS = {'comparator.py': {'bytes': 49045, 'sha256': 'ed2f0172df1fd361fde3cc6951feacc68118308c2f5c1f500c9919f5453a0ff0'}, 'protocol_helpers.py': {'bytes': 1523, 'sha256': '950525ba84b706d9e55824cebfb344484295a44f904ef3f652d364efd98b34a3'}, 'native_producer_core.jl': {'bytes': 9710, 'sha256': 'fda7ec19080ef823f6000693fa78c7ca0a8da24ee0f204340fe85ad658515ff9'}, 'tenpy_producer_core.py': {'bytes': 14971, 'sha256': '65835eb811f9f836cbccfc17a0e94d9de83056aadf48aff3baa8aad0f96ba99a'}}
SOURCE_INVENTORY_SHA256 = "4d6b29741c71f65bbe5fae63354545e785961723cef889a6717300882cda369a"
GATE_CALL_SITES = [217, 228, 282, 287, 292, 295, 303, 304, 306, 311, 320, 322, 330, 334, 337, 345, 355, 358, 360, 361, 362, 366, 373, 374, 398, 401, 405, 428, 432, 434, 458, 462, 466, 469, 471, 472, 473, 475, 477, 481, 485, 488, 489, 494, 503, 508, 509, 512, 521, 523, 528, 558, 562, 563, 566, 571, 581, 583, 586, 598, 604, 620]

def _site_217(key,name,site):
    return {'name': f'{name}:site{site}:local state {key}', 'kind': 'physical'}

def _site_228(name,site):
    return {'name': f'{name}:site{site}:complete exact local CAR/number/spin', 'kind': 'physical'}

def _site_282(label):
    return {'name': label + ':fixed full six-leg physical pipe', 'kind': 'eligibility'}

def _site_287(label):
    return {'name': label + ':pipe slices/qflat dimensions', 'kind': 'eligibility'}

def _site_292(label,site):
    return {'name': label + f':incoming leg{site} slices', 'kind': 'eligibility'}

def _site_295(label,site):
    return {'name': label + f':incoming leg{site} physical charges', 'kind': 'eligibility'}

def _site_303(label):
    return {'name': label + ':whole pipe charge fusion matches actual qflat', 'kind': 'eligibility'}

def _site_304(label):
    return {'name': label + ':all4096 physical digit tuples form a full pipe bijection', 'kind': 'eligibility'}

def _site_306(label):
    return {'name': label + ':actual balanced charge mask exactly400', 'kind': 'eligibility'}

def _site_311(label):
    return {'name': label + ':actual canonical projection reproduced from raw pipe', 'kind': 'eligibility'}

def _site_320(backend):
    return {'name': backend + ':complete Ui8/Uf4 pair', 'kind': 'eligibility'}

def _site_322(backend):
    return {'name': backend + ':complete actual operator table pairs', 'kind': 'eligibility'}

def _site_330(label):
    return {'name': label + ':complete real sorted400 eigensystem', 'kind': 'eligibility'}

def _site_334(label):
    return {'name': label + ':complete physical Hamiltonian entries', 'kind': 'physical', 'limit': 1e-11}

def _site_337(label):
    return {'name': label + ':actual full MPO all400x400 entries', 'kind': 'physical', 'limit': 1e-11}

def _site_345():
    return {'name': 'tenpy:initial/final actual projection and local phase compatibility', 'kind': 'eligibility'}

def _site_355(backend):
    return {'name': backend + ':complete original initial400vector', 'kind': 'eligibility'}

def _site_358(backend):
    return {'name': backend + ':finite positive reported original eigenvector norm', 'kind': 'eligibility'}

def _site_360(backend):
    return {'name': backend + ':original eigenvector norm association', 'kind': 'precision', 'limit': 1e-12}

def _site_361(backend):
    return {'name': backend + ':actual used initial norm', 'kind': 'precision', 'limit': 1e-12}

def _site_362(backend):
    return {'name': backend + ':actual initial normalization formula', 'kind': 'precision', 'limit': 1e-10}

def _site_366(backend):
    return {'name': backend + ':qualified initial same-sector gap', 'kind': 'precision', 'required_minimum': 0.001}

def _site_373(backend):
    return {'name': backend + ':actual initial H residual', 'kind': 'precision', 'limit': 1e-10}

def _site_374(backend):
    return {'name': backend + ':Hermitian Rayleigh imaginary enclosure contains zero', 'kind': 'eligibility'}

def _site_398(label):
    return {'name': label + ':complete post-entropy MPS export stage', 'kind': 'eligibility'}

def _site_401(label):
    return {'name': label + ':post-MPS finite positive scalar and physical bond budget', 'kind': 'eligibility'}

def _site_405(label):
    return {'name': label + ':actual allB/unit TeNPy endpoint association', 'kind': 'eligibility'}

def _site_428(label):
    return {'name': label + ':all4096 original complete amplitudes retained', 'kind': 'eligibility'}

def _site_432(label):
    return {'name': label + ':independent all4096 tensor/export identity', 'kind': 'precision', 'limit': 1e-10}

def _site_434(label):
    return {'name': label + ':actual post-MPS raw norm', 'kind': 'precision', 'limit': 1e-12}

def _site_458(label):
    return {'name': label + ':exact original time and arm', 'kind': 'eligibility'}

def _site_462(label):
    return {'name': label + ':registered native cutoff/bond/normalization call', 'kind': 'eligibility'}

def _site_466(label):
    return {'name': label + ':registered actual from_full cutoff/external budget', 'kind': 'eligibility'}

def _site_469(label):
    return {'name': label + ':all400 actual raw/normalized input states', 'kind': 'eligibility'}

def _site_471(label):
    return {'name': label + ':finite positive original raw evolution norm', 'kind': 'eligibility'}

def _site_472(label):
    return {'name': label + ':original raw norm association', 'kind': 'precision', 'limit': 1e-12}

def _site_473(label):
    return {'name': label + ':actual raw evolved norm', 'kind': 'precision', 'limit': 1e-12}

def _site_475(label):
    return {'name': label + ':actual normalized MPS input norm', 'kind': 'precision', 'limit': 1e-12}

def _site_477(label):
    return {'name': label + ':actual raw-to-input normalization formula', 'kind': 'precision', 'limit': 1e-10}

def _site_481(label):
    return {'name': label + ':actual direct complete spectral action identity', 'kind': 'precision', 'limit': 1e-10}

def _site_485(label):
    return {'name': label + ':actual all-bond rank observed within registered budgets', 'kind': 'eligibility'}

def _site_488(label):
    return {'name': label + ':post-normalize norm API association', 'kind': 'precision', 'limit': 1e-12}

def _site_489(label):
    return {'name': label + ':actual user Schmidt normcheck association', 'kind': 'precision', 'limit': 1e-12}

def _site_494(label):
    return {'name': label + ':complete normalized pre/post MPS state distance', 'kind': 'precision', 'limit': 1e-10}

def _site_503(label):
    return {'name': label + ':finite observed Schmidt values and physical rank', 'kind': 'eligibility'}

def _site_508(label):
    return {'name': label + ':observed Schmidt sum normalization', 'kind': 'precision', 'limit': 1e-12}

def _site_509(label):
    return {'name': label + ':new observed Schmidt purity vs full tensor state', 'kind': 'precision', 'limit': 1e-10}

def _site_512(label):
    return {'name': label + ':actual reported Schmidt probability sum', 'kind': 'precision', 'limit': 1e-12}

def _site_521(label):
    return {'name': label + ':finite original actual library S2', 'kind': 'eligibility'}

def _site_523(label):
    return {'name': label + ':actual library entropy vs own normalized full state plus budgets', 'kind': 'physical', 'limit': 1e-10}

def _site_528(label):
    return {'name': label + ':primary total entropy radius', 'kind': 'precision', 'limit': 4e-08}

def _site_558(backend):
    return {'name': backend + ':source-defined complete producer scope', 'kind': 'eligibility'}

def _site_562(backend,stage):
    return {'name': backend + ':saved actual runtime ' + stage, 'kind': 'eligibility'}

def _site_563():
    return {'name': 'tenpy:source-recorded initial/final physical mapping equality', 'kind': 'eligibility'}

def _site_566(backend):
    return {'name': backend + ':complete unique registered time-arm coverage', 'kind': 'eligibility'}

def _site_571(backend,k):
    return {'name': backend + f':k{k}:same original pre-conversion evolution in both arms', 'kind': 'eligibility'}

def _site_581(k):
    return {'name': f'k{k}:complete primary S2 cross-backend', 'kind': 'physical', 'limit': 1e-07}

def _site_583(k):
    return {'name': f'k{k}:complete primary purity cross-backend', 'kind': 'physical', 'limit': 1e-08}

def _site_586(backend,k):
    return {'name': f'{backend}:k{k}:primary/coarse S2 stability', 'kind': 'physical', 'limit': 1e-07}

def _site_598():
    return {'name': 'fixed nine-point primary growth slope cross-backend', 'kind': 'physical', 'limit': 1e-06}

def _site_604(backend):
    return {'name': backend + ':nontrivial registered full curve', 'kind': 'exercise', 'required_lower': 0.0001}

def _site_620():
    return {'name': 'full exact registered comparison coverage', 'kind': 'eligibility'}

def expected_full_gate_schedule():
    """Return ordered {name,kind,limit?/required_minimum?/required_lower?} metadata only."""
    out = []
    for backend in ('native','tenpy'):
        out.append(_site_558(backend))
        for stage in ('before','after'):
            out.append(_site_562(backend,stage))
        if backend == 'tenpy':
            out.append(_site_563())
        out.append(_site_320(backend))
        out.append(_site_322(backend))
        for U in (8.0,4.0):
            label = f'{backend}:U{U}'
            name = label
            for site in range(6):
                for key in ((0,0),(1,0),(0,1),(1,1)):
                    out.append(_site_217(key,name,site))
                out.append(_site_228(name,site))
            if backend == 'tenpy':
                out.append(_site_282(label))
                out.append(_site_287(label))
                for site in range(6):
                    out.append(_site_292(label,site))
                    out.append(_site_295(label,site))
                for _pipe_tuple in range(4096):
                    out.append(_site_303(label))
                out.append(_site_304(label))
                out.append(_site_306(label))
                out.append(_site_311(label))
            out.append(_site_330(label))
            out.append(_site_334(label))
            out.append(_site_337(label))
        if backend == 'tenpy':
            out.append(_site_345())
        out.append(_site_355(backend))
        out.append(_site_358(backend))
        out.append(_site_360(backend))
        out.append(_site_361(backend))
        out.append(_site_362(backend))
        out.append(_site_366(backend))
        out.append(_site_373(backend))
        out.append(_site_374(backend))
        out.append(_site_566(backend))
        for k in range(33):
            out.append(_site_571(backend,k))
        for arm in ('coarse','primary'):
            for k in range(33):
                label = f'{backend}:k{k}:{arm}'
                out.append(_site_458(label))
                if backend == 'native':
                    out.append(_site_462(label))
                else:
                    out.append(_site_466(label))
                out.append(_site_469(label))
                out.append(_site_471(label))
                out.append(_site_472(label))
                out.append(_site_473(label))
                out.append(_site_475(label))
                out.append(_site_477(label))
                out.append(_site_481(label))
                out.append(_site_398(label))
                out.append(_site_401(label))
                if backend == 'tenpy':
                    out.append(_site_405(label))
                out.append(_site_428(label))
                out.append(_site_432(label))
                out.append(_site_434(label))
                out.append(_site_485(label))
                if backend == 'native':
                    out.append(_site_488(label))
                    out.append(_site_489(label))
                out.append(_site_494(label))
                out.append(_site_503(label))
                out.append(_site_508(label))
                out.append(_site_509(label))
                if backend == 'tenpy':
                    out.append(_site_512(label))
                out.append(_site_521(label))
                out.append(_site_523(label))
                if arm == 'primary':
                    out.append(_site_528(label))
    for k in range(33):
        out.append(_site_581(k))
        out.append(_site_583(k))
        for backend in ('native','tenpy'):
            out.append(_site_586(backend,k))
    out.append(_site_598())
    for backend in ('native','tenpy'):
        out.append(_site_604(backend))
    out.append(_site_620())
    return out
