"""Prospective P13 full-state/entropy comparator; no scientific imports/run.

A reviewed allocation wrapper must inject qualified NumPy and mandatory runtime
and checkpoint callbacks. All results remain uncredited until independent audit.
Uses conditional P17 balls, not directed-rounding interval arithmetic.
"""
import copy
import math
import hashlib
import inspect
from pathlib import Path
from itertools import product
from protocol_helpers import L,CUT,ARMS,coordinates,regroup_sign,canonical_site_linear_index
from inherited_p17_certificate import Arithmetic,RealBall,ComplexBall,eigensystem_certificate
from evolution_bounds import (spectral_action_ball,eigensystem_operator_error,initial_ground_state_distance,
    log_positive,entropy_from_purity_with_state_budget,growth_rate_ball)

LIMITS={'matrix':1e-11,'gap':1e-3,'norm':1e-12,'initial_residual':1e-10,'spectral_identity':1e-10,
    'conversion':1e-10,'library_entropy':1e-10,'primary_radius':4e-8,'cross_entropy':1e-7,
    'cross_purity':1e-8,'slope':1e-6,'arm_entropy':1e-7,'curve_range':1e-4}


def callable_witness(name,function):
 actual=getattr(function,'__func__',function)
 module=inspect.getmodule(actual)
 try:source=inspect.getsourcefile(actual)
 except (TypeError,OSError):source=None
 try:line=inspect.getsourcelines(actual)[1] if source else None
 except (TypeError,OSError):line=None
 module_path=getattr(module,'__file__',None)if module is not None else None
 path=source or module_path
 result={'name':name,'callable':callable(actual),'callable_type':type(actual).__module__+'.'+type(actual).__qualname__,'module':getattr(actual,'__module__',None),'qualname':getattr(actual,'__qualname__',getattr(actual,'__name__',None)),'python_source_available':source is not None,'source_line':line}
 if path is not None:
  real=Path(path).resolve(strict=True);data=real.read_bytes()
  result.update(file=str(real),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
 else:result['file_resolution_scope']='No standalone Python source path; qualified array-extension/provider records remain mandatory.'
 return result

def api_witnesses(np):
    entries=[('Arithmetic.__init__',Arithmetic.__init__),('eigensystem_certificate',eigensystem_certificate),
        ('spectral_action_ball',spectral_action_ball),('log_positive',log_positive),('growth_rate_ball',growth_rate_ball),
        ('numpy.asarray',np.asarray),('numpy.dot',np.dot),('numpy.matmul',np.matmul)]
    result=[callable_witness(name,function)for name,function in entries]
    extension=np._core._multiarray_umath
    path=Path(extension.__file__).resolve(strict=True);data=path.read_bytes()
    result.append({'name':'actual_NumPy_array_extension','module':extension.__name__,'file':str(path),
        'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),
        'scope':'Actual compiled array provider identity; no invented standalone Python source for ufunc.'})
    return result


def unpack(value,np):
    if not isinstance(value,dict) or set(value)!={'shape','order','real','imag'}:
        raise RuntimeError('Exact finite complex pack schema required')
    shape=value['shape']
    if (not isinstance(shape,list) or any(type(x)is not int or x<0 or x>4096 for x in shape) or
        value['order']not in('F','C') or math.prod(shape)>160000 or
        len(value['real'])!=math.prod(shape) or len(value['imag'])!=math.prod(shape)):
        raise RuntimeError('Unqualified complex pack dimensions/order')
    real=np.asarray(value['real'],dtype=np.float64).reshape(shape,order=value['order'])
    imag=np.asarray(value['imag'],dtype=np.float64).reshape(shape,order=value['order'])
    result=np.empty(shape,dtype=np.complex128);result.real=real;result.imag=imag
    if not np.all(np.isfinite(result)):raise RuntimeError('Nonfinite original physical export')
    return result


def physical_entries():
    coords=coordinates();index={pair:i for i,pair in enumerate(coords)};entries=[]
    for column,(up,down)in enumerate(coords):
        entries.append((column,column,(up&down).bit_count()))
        occupied=up|(down<<L)
        for offset in(0,L):
            for j in range(L-1):
                for destination,source in((offset+j,offset+j+1),(offset+j+1,offset+j)):
                    if (occupied>>source)&1 and not((occupied>>destination)&1):
                        after=occupied^(1<<source)
                        sign=(-1 if(occupied&((1<<source)-1)).bit_count()%2 else 1)
                        sign*=(-1 if(after&((1<<destination)-1)).bit_count()%2 else 1)
                        target=after^(1<<destination)
                        row=index[target&63,target>>L]
                        entries.append((row,column,-sign))
    return entries


def validate_registration(raw,backend,pilot):
    """Pure metadata coverage/controls, independent of any numeric package."""
    schemas={'native':'P13_NATIVE_MPO_ENTROPY_DENSE_ADAPTER_CORE_V2','tenpy':'P13_TENPY_EXACTDIAG_ENTROPY_PRODUCER_CORE_V2'}
    if backend not in schemas or type(pilot)is not bool:raise RuntimeError('Unknown backend/mode')
    if (raw.get('schema')!=schemas[backend]or raw.get('status')!='PRODUCER_COMPLETE_NOT_A_SCIENTIFIC_PASS'or
        type(raw.get('scientific_credit'))is not int or raw['scientific_credit']!=0 or
        type(raw.get('pilot_only'))is not bool or raw['pilot_only']!=pilot or
        raw.get('coordinates')!=[list(x)for x in coordinates()]):raise RuntimeError('Original producer registration differs')
    for stage in('before','after'):
        if not isinstance(raw.get('runtime_'+stage),dict)or raw['runtime_'+stage].get('status')!='PASS':
            raise RuntimeError('Complete actual producer runtime observation missing')
    if len(raw.get('models',[]))!=2 or [x.get('U')for x in raw['models']]!=[8.0,4.0]:
        raise RuntimeError('Both registered original Hamiltonians required')
    times=range(5 if pilot else 33);rows=raw.get('samples',[]);keys=[]
    for row in rows:
        k=row['time_numerator'];arm=row['arm']
        if type(k)is not int or k not in times or arm not in ARMS or type(row['time_denominator'])is not int or row['time_denominator']!=8 or row['time']!=k/8:
            raise RuntimeError('Sample has unregistered time/arm')
        if backend=='native':
            if row['maxdim']!=ARMS[arm]['native_maxdim']or row['cutoff']!=ARMS[arm]['native_cutoff']:
                raise RuntimeError('Native original controls differ')
        elif row['external_bond_budget']!=ARMS[arm]['tenpy_budget']or row['from_full_cutoff']!=ARMS[arm]['tenpy_cutoff']:
            raise RuntimeError('TeNPy original controls differ')
        keys.append((k,arm))
    if len(keys)!=2*len(times)or len(set(keys))!=len(keys)or set(keys)!={(k,arm)for k in times for arm in ARMS}:
        raise RuntimeError('Incomplete/repeated registered time-arm samples')
    return dict(zip(keys,rows))


def pipe_forward_factory(legs,slices,qmap):
    """Exact raw-table replay of LegPipe.map_incoming_flat (no library call)."""
    if (not isinstance(slices,list)or not slices or slices[0]!=0 or any(type(x)is not int for x in slices)or
        any(a>=b for a,b in zip(slices,slices[1:]))):raise RuntimeError('Invalid outgoing pipe slices')
    for leg in legs:
        ends=leg['slices']
        if (not isinstance(ends,list)or len(ends)<2 or ends[0]!=0 or any(type(x)is not int for x in ends)or
            any(a>=b for a,b in zip(ends,ends[1:]))):raise RuntimeError('Invalid incoming pipe slices')
    by_q={}
    for row in qmap:
        if not isinstance(row,list)or len(row)!=3+len(legs)or any(type(x)is not int for x in row):raise RuntimeError('Invalid raw q_map row')
        key=tuple(row[3:])
        if key in by_q or any(not 0<=q<len(legs[j]['slices'])-1 for j,q in enumerate(key)):
            raise RuntimeError('Repeated/out-of-range incoming q_map tuple')
        by_q[key]=row
    if len(by_q)!=math.prod(len(leg['slices'])-1 for leg in legs):raise RuntimeError('Incomplete raw q_map block grid')
    def forward(digits):
        if len(digits)!=len(legs):raise RuntimeError('Incoming digit tuple length differs')
        qindices=[];within=0;stride=1
        for site in range(len(legs)-1,-1,-1):
            ends=legs[site]['slices'];digit=digits[site]
            if type(digit)is not int or not 0<=digit<ends[-1]:raise RuntimeError('Incoming physical digit outside leg')
            qindex=next(i for i,(lo,hi)in enumerate(zip(ends,ends[1:]))if lo<=digit<hi)
            qindices.insert(0,qindex);within+=stride*(digit-ends[qindex]);stride*=ends[qindex+1]-ends[qindex]
        row=by_q[tuple(qindices)];lo,hi,out=row[:3]
        if not(0<=out<len(slices)-1 and hi-lo==stride and 0<=lo<hi<=slices[out+1]-slices[out]):
            raise RuntimeError('Raw pipe block displacement/size is invalid')
        return slices[out]+lo+within
    return forward


class Comparison:
    def __init__(self,np,a,report,on_record):
        self.np=np;self.a=a;self.report=report;self.on_record=on_record
    def checkpoint(self,phase):
        self.report['phase']=phase;self.on_record('checkpoint',copy.deepcopy(self.report))
    def gate(self,name,ok,kind='physical',**details):
        row={'name':name,'status':'PASS'if bool(ok)else'FAIL','kind':kind,**details}
        self.report['gates'].append(row);return bool(ok)
    def need(self,name,ok,**details):
        if not self.gate(name,ok,'eligibility',**details):raise RuntimeError(name)
    def bound(self,x):
        return self.a.scalar_upper(self.a.up(x))
    def upper_abs(self,x):
        return self.a.inflate(self.np.abs(x.c)+x.r,128)
    def complex_abs(self,x):
        return self.a.inflate(self.upper_abs(x.real)+self.upper_abs(x.imag),128)
    def maximum(self,x):
        return self.bound(self.np.max(x))
    def cb(self,x):
        return ComplexBall.exact(self.a,x)
    def reshape(self,x,shape):
        return x.reshape(shape)
    def error(self,name,value,limit,kind='precision',**details):
        scalar=self.a.check(value)
        if scalar.shape!=()or scalar<0:raise RuntimeError('Nonnegative scalar error required')
        self.gate(name,scalar<=limit,kind,value_upper=self.bound(scalar),limit=limit,
            qualified_longdouble_value=str(scalar),**details)
        return scalar<=limit
    def serialize_ball(self,value):
        center=float(value.c)
        if not math.isfinite(center):raise RuntimeError('Nonfinite serialized center')
        radius=self.upper_abs(self.a.ball(center)-value)
        return {'center':center,'radius':self.bound(radius)}
    def physical_model(self,U):
        np=self.np;matrix=np.zeros((400,400),dtype=np.float64)
        for row,column,value in physical_entries():
            matrix[row,column]+=U*value if row==column else value
        # U8/U4 and all hopping/occupation integers are exactly represented;
        # each finite matrix entry is an exact small integer before ball work.
        return self.a.ball(matrix)
    def local_layout(self, local, name):
        """Exact 0,+/-1 local CAR is required; no fitted labels or phase rounding."""
        np = self.np
        if len(local) != L:
            raise RuntimeError('All six actual local operator tables are mandatory')
        layouts = []
        for site, table in enumerate(local):
            op = {key: unpack(value, np) for key, value in table.items()}
            required = {'Id','Nup','Ndn','Ntot','Sz','Cdagup','Cdagdn','Cup','Cdn','F'}
            if set(op) != required or any(x.shape != (4,4) for x in op.values()):
                raise RuntimeError('Incomplete physical local operator export')
            if any(not np.all(np.isfinite(x)) for x in op.values()):
                raise RuntimeError('Nonfinite local operator')
            occupations = {}
            for j in range(4):
                u, d = op['Nup'][j,j], op['Ndn'][j,j]
                if u not in (0,1) or d not in (0,1):
                    raise RuntimeError('Actual occupation labels are not exact zero/one')
                key = (int(u.real), int(d.real))
                if key in occupations:
                    raise RuntimeError('Repeated actual local occupation label')
                occupations[key] = j
            if set(occupations) != {(0,0),(1,0),(0,1),(1,1)}:
                raise RuntimeError('Incomplete occupation bijection')
            vac = np.zeros(4); vac[occupations[0,0]] = 1
            vectors = {(0,0):vac, (1,0):op['Cdagup']@vac,
                       (0,1):op['Cdagdn']@vac, (1,1):op['Cdagup']@op['Cdagdn']@vac}
            phases = {key:value[occupations[key]] for key,value in vectors.items()}
            if any(x not in (-1,1) for x in phases.values()):
                raise RuntimeError('Unqualified actual local fermionic phase')
            for key, vector in vectors.items():
                ideal = np.zeros(4,dtype=np.complex128); ideal[occupations[key]]=phases[key]
                self.gate(f'{name}:site{site}:local state {key}', np.array_equal(vector,ideal))
            ident=np.eye(4);zero=np.zeros((4,4));cu=op['Cup'];cd=op['Cdn']
            conditions = [op['Id']==ident, op['Cup']==op['Cdagup'].conj().T,
                op['Cdn']==op['Cdagdn'].conj().T,
                cu@op['Cdagup']+op['Cdagup']@cu==ident,
                cd@op['Cdagdn']+op['Cdagdn']@cd==ident,
                cu@cd+cd@cu==zero,cu@op['Cdagdn']+op['Cdagdn']@cu==zero,
                op['Nup']==op['Cdagup']@cu,op['Ndn']==op['Cdagdn']@cd,
                op['Ntot']==op['Nup']+op['Ndn'],
                op['Sz']==(op['Nup']-op['Ndn'])/2,
                op['F']==ident-2*op['Ntot']+4*(op['Nup']@op['Ndn'])]
            self.gate(f'{name}:site{site}:complete exact local CAR/number/spin',
                      all(np.all(x) for x in conditions))
            # Requiring exact sparse integer/half-integer local identities makes
            # these small binary64 controls exact; no numerical phase fit occurs.
            if not all(np.all(x) for x in conditions):
                raise RuntimeError('Actual local CAR/operator identity failed')
            layouts.append((occupations,phases))
        return layouts

    def mapping(self, coords, layouts):
        np=self.np
        digits=[];phases=[]
        for u,d in coords:
            occ=[((u>>j)&1,(d>>j)&1) for j in range(L)]
            digits.append([layouts[j][0][occ[j]] for j in range(L)])
            p=regroup_sign(u,d)
            for j in range(L): p*=layouts[j][1][occ[j]]
            phases.append(p)
        return np.asarray(digits,dtype=np.int64),np.asarray(phases,dtype=np.complex128)

    def contract_mpo(self, exported, coords, layouts):
        """All entries, including zeros, with batched split-real error balls."""
        np=self.np;a=self.a
        if exported.get('axes') != ['left','right','bra','ket'] or len(exported['sites'])!=L:
            raise RuntimeError('Unqualified actual MPO axes')
        tensors=[unpack(x,np) for x in exported['sites']]
        for j,t in enumerate(tensors):
            if t.ndim!=4 or t.shape[2:]!=(4,4) or (j and tensors[j-1].shape[1]!=t.shape[0]):
                raise RuntimeError('Inconsistent actual MPO tensor shape')
        left=exported.get('left_boundary_index',0)
        right=exported.get('right_boundary_index',0)
        if not (isinstance(left,int) and 0<=left<tensors[0].shape[0] and
                isinstance(right,int) and 0<=right<tensors[-1].shape[1]):
            raise RuntimeError('Invalid actual MPO boundary identity index')
        if 'left_boundary_index' not in exported and tensors[0].shape[0]!=1:
            raise RuntimeError('Missing native singleton left boundary')
        if 'right_boundary_index' not in exported and tensors[-1].shape[1]!=1:
            raise RuntimeError('Missing native singleton right boundary')
        digits,phases=self.mapping(coords,layouts);n=len(coords)
        cr=np.zeros((n,n),dtype=a.ld);ci=cr.copy();rr=cr.copy();ri=cr.copy()
        for row in range(n):
            initial=np.zeros((n,tensors[0].shape[0]));initial[:,left]=1
            state=self.cb(initial)
            for j,tensor in enumerate(tensors):
                selected=np.transpose(tensor[:,:,digits[row,j],digits[:,j]],(2,0,1))
                product=self.reshape(state,(n,state.real.c.shape[1],1))*self.cb(selected)
                state=product.sum(axis=1)
            value=state[:,right]*self.cb(phases[row].conjugate()*phases)
            cr[row]=value.real.c;ci[row]=value.imag.c;rr[row]=value.real.r;ri[row]=value.imag.r
        return ComplexBall(RealBall(a,cr,rr),RealBall(a,ci,ri))

    def qualify_pipe(self,model,layouts,label):
        """Replay raw LegPipe tables with exact integers; no private sort guess."""
        pipe=model['actual_pipe_mapping'];legs=pipe['incoming_legs']
        self.need(label+':fixed full six-leg physical pipe',pipe['nlegs']==6 and pipe['ind_len']==4096 and
            pipe['subshape']==[4]*6 and len(legs)==6 and pipe['qconj']in(-1,1) and
            pipe['chinfo_mod']==[1,1] and pipe['chinfo_names']==['N','2*Sz'])
        def integers(values):return isinstance(values,list)and all(type(x)is int for x in values)
        slices=pipe['slices'];qflat=pipe['qflat'];qmap=pipe['q_map']
        self.need(label+':pipe slices/qflat dimensions',integers(slices)and slices[0]==0 and slices[-1]==4096 and
            all(a<b for a,b in zip(slices,slices[1:]))and len(qflat)==4096 and
            all(integers(row)and len(row)==2 for row in qflat))
        for site,leg in enumerate(legs):
            ends=leg['slices'];flat=leg['qflat'];occ=layouts[site][0]
            self.need(label+f':incoming leg{site} slices',leg['ind_len']==4 and leg['qconj']in(-1,1)and
                integers(ends)and ends[0]==0 and ends[-1]==4 and all(a<b for a,b in zip(ends,ends[1:]))and
                len(flat)==4 and all(integers(row)and len(row)==2 for row in flat))
            self.need(label+f':incoming leg{site} physical charges',
                all(flat[index]==[up+down,up-down]for(up,down),index in occ.items())and
                all(all(flat[j]==flat[lo]for j in range(lo,hi))for lo,hi in zip(ends,ends[1:])))
        forward=pipe_forward_factory(legs,slices,qmap)
        outgoing=[]
        for digits in product(range(4),repeat=6):
            flat=forward(digits);outgoing.append(flat)
            charge=[pipe['qconj']*sum(legs[j]['qconj']*legs[j]['qflat'][digits[j]][q]for j in range(6))for q in range(2)]
            self.need(label+':whole pipe charge fusion matches actual qflat',0<=flat<4096 and qflat[flat]==charge)
        self.need(label+':all4096 physical digit tuples form a full pipe bijection',sorted(outgoing)==list(range(4096)))
        selected=[i for i,q in enumerate(qflat)if q==[6,0]]
        self.need(label+':actual balanced charge mask exactly400',len(selected)==400 and
            model['actual_charge_sector']==[6,0]and model['actual_charge_names']==['N','2*Sz']and
            model['actual_charge_mask']==selected and model['actual_pipe_charge_table']==qflat)
        digits,phases=self.mapping(coordinates(),layouts)
        derived=[selected.index(forward(tuple(int(x)for x in row)))for row in digits]
        self.need(label+':actual canonical projection reproduced from raw pipe',
            sorted(derived)==list(range(400))and model['actual_canonical_to_selected_permutation']==derived and
            model['actual_canonical_local_digits']==digits.tolist()and
            self.np.array_equal(unpack(model['actual_canonical_phases'],self.np),phases))
        return {'full_pipe_bijection_size':4096,'balanced_projection_size':400,
            'projection':derived,'mask':selected,'digits':digits.tolist(),'phases':model['actual_canonical_phases']}

    def qualify_models(self,raw,backend):
        np=self.np;a=self.a;models=raw['models']
        self.need(backend+':complete Ui8/Uf4 pair',len(models)==2 and [x['U']for x in models]==[8.0,4.0])
        local=(raw['local_operators_initial_and_final']if backend=='tenpy'else[raw['local_operators']]*2)
        self.need(backend+':complete actual operator table pairs',len(local)==2)
        result=[];projection=[]
        for index,(U,model,tables)in enumerate(zip((8.0,4.0),models,local)):
            label=f'{backend}:U{U}';self.checkpoint(label+':before_full_model_certificate')
            layouts=self.local_layout(tables,label)
            if backend=='tenpy':projection.append(self.qualify_pipe(model,layouts,label))
            H=self.physical_model(U);actual=unpack(model['Hamiltonian'],np)
            V=unpack(model['eigenvectors'],np);values=np.asarray(model['eigenvalues'],dtype=np.float64)
            self.need(label+':complete real sorted400 eigensystem',actual.shape==V.shape==(400,400)and
                values.shape==(400,)and np.all(np.isfinite(values))and np.all(V.imag==0)and np.all(actual.imag==0)and
                np.all(values[1:]>=values[:-1]))
            matrix_error=self.np.max(self.complex_abs(self.cb(actual)-ComplexBall(H,a.ball(np.zeros((400,400))))))
            self.error(label+':complete physical Hamiltonian entries',matrix_error,LIMITS['matrix'],'physical')
            represented=self.contract_mpo(model['represented_mpo'],coordinates(),layouts)
            mpo_error=self.np.max(self.complex_abs(represented-ComplexBall(H,a.ball(np.zeros((400,400))))))
            self.error(label+':actual full MPO all400x400 entries',mpo_error,LIMITS['matrix'],'physical')
            certificate=eigensystem_certificate(a,H,V.real,values)
            record={'backend':backend,'U':U,'dimension':400,'certificate':{k:int(v)if k=='dimension'else self.bound(v)for k,v in certificate.items()},
                'matrix_error_upper':self.bound(matrix_error),'actual_MPO_error_upper':self.bound(mpo_error)}
            self.report['model_certificates'].append(record)
            result.append({'H':H,'V':V.real,'values':values,'certificate':certificate,'layouts':layouts})
            self.checkpoint(label+':after_full_model_certificate')
        if backend=='tenpy':
            self.need('tenpy:initial/final actual projection and local phase compatibility',projection[0]==projection[1])
        return result

    def normalize(self,vector):
        norm2=vector.norm_squared()
        if norm2.lower()<=0:raise RuntimeError('Actual vector norm interval not strictly positive')
        return vector.divide_real(norm2.sqrt()),norm2

    def initial_certificate(self,raw,model,backend):
        np=self.np;a=self.a;initial=unpack(raw['actual_initial_state'],np)
        self.need(backend+':complete original initial400vector',initial.shape==(400,))
        vector=self.cb(initial);column=self.cb(model['V'][:,0]);norm2=vector.norm_squared()
        original_norm=raw['original_initial_eigenvector_norm']
        self.need(backend+':finite positive reported original eigenvector norm',type(original_norm)in(int,float)and math.isfinite(original_norm)and original_norm>0)
        column_norm=column.norm_squared().sqrt()
        self.error(backend+':original eigenvector norm association',self.upper_abs(column_norm-original_norm),LIMITS['norm'])
        self.error(backend+':actual used initial norm',self.upper_abs(norm2-1),LIMITS['norm'])
        self.error(backend+':actual initial normalization formula',
            (vector-column.divide_real(a.ball(original_norm))).norm_upper(),LIMITS['conversion'])
        normalization=(vector-column).norm_upper()
        ground=initial_ground_state_distance(a,model['certificate'],model['values'])
        self.gate(backend+':qualified initial same-sector gap',ground['gap_lower']>=LIMITS['gap'],'precision',
            gap_lower_decimal=str(ground['gap_lower']),required_minimum=LIMITS['gap'])
        H=ComplexBall(model['H'],a.ball(np.zeros((400,400))))
        normalized_initial,_=self.normalize(vector)
        action=H.matmul(normalized_initial.reshape((400,1)))
        energy=normalized_initial.conjugate().reshape((1,400)).matmul(action)[0,0]
        residual=(action-normalized_initial.reshape((400,1))*energy.real).norm_upper()
        self.error(backend+':actual initial H residual',residual,LIMITS['initial_residual'])
        distance=(a.ball(ground['recorded_ground_distance_upper'])+a.ball(normalization)).upper()
        record={'backend':backend,'actual_initial_norm_squared':self.serialize_ball(norm2),
            'actual_initial_norm_upper':self.bound(norm2.sqrt().upper()),'original_eigenvector_norm_reported':original_norm,
            'recorded_column_to_used_initial_error_upper':self.bound(normalization),'ground_identification_upper':self.bound(ground['recorded_ground_distance_upper']),
            'initial_physical_state_distance_upper':self.bound(distance),'initial_residual_upper':self.bound(residual),
            'initial_residual_scope':'Exact-normalized actual initial ket; original used initial remains unchanged in spectral action and norm/error budgets.',
            'phase_scope':'Ground-state ray identified; no raw cross-backend phase fitting or final eigenvector alignment.'}
        self.report['initial_state_certificates'].append(record)
        return {'vector':initial,'distance':distance,'norm_upper':norm2.sqrt().upper()}

    def reconstruct_all_sites(self,sample,layouts,label):
        np=self.np;a=self.a;export=sample['state']
        self.need(label+':complete post-entropy MPS export stage',export.get('axes')==['left','right','physical']and
            export.get('export_stage')in('AFTER_ACTUAL_USER_ENTROPY_CALL','AFTER_ACTUAL_TENPY_ENTROPY_CALL')and len(export['sites'])==6)
        tensors=[unpack(row,np)for row in export['sites']];scalar=export['external_norm_scalar']
        self.need(label+':post-MPS finite positive scalar and physical bond budget',type(scalar)in(int,float)and math.isfinite(scalar)and scalar>0 and
            all(t.ndim==3 and t.shape[2]==4 and 0<t.shape[0]<=64 and 0<t.shape[1]<=64 for t in tensors)and
            tensors[0].shape[0]==tensors[-1].shape[1]==1 and all(tensors[j].shape[1]==tensors[j+1].shape[0]for j in range(5)))
        if 'actual_canonical_forms'in export:
            self.need(label+':actual allB/unit TeNPy endpoint association',export['actual_canonical_forms']==[[0.0,1.0]]*6 and
                export['actual_left_schmidt']==export['actual_right_schmidt']==[1.0]and scalar==1.0 and
                export['scalar_placement']=='UNSCALED actual B tensors; actual psi.norm multiplied at END of complete contraction')
        all_coords=tuple((u,d)for u in range(64)for d in range(64))
        digits,group_phases=self.mapping(all_coords,layouts)
        site_indices=np.asarray([canonical_site_linear_index(u,d)for u,d in all_coords],dtype=np.int64)
        actual_indices=np.asarray([sum(int(row[j])*4**j for j in range(6))for row in digits],dtype=np.int64)
        # Product of local canonical creation phases maps actual site basis to
        # canonical site basis. Grouped-spin reordering is used only for the
        # 400-state spectral input embedding, never to drop outside-sector data.
        local_phases=group_phases*np.asarray([regroup_sign(u,d)for u,d in all_coords])
        cr=np.zeros(4096,dtype=a.ld);ci=cr.copy();rr=cr.copy();ri=cr.copy()
        for start in range(0,4096,128):
            stop=min(start+128,4096);batch=digits[start:stop];rows=stop-start
            x=self.cb(np.ones((rows,1)))
            for j,tensor in enumerate(tensors):
                selected=np.transpose(tensor[:,:,batch[:,j]],(2,0,1))
                x=(x.reshape((rows,x.real.c.shape[1],1))*self.cb(selected)).sum(axis=1)
            value=x[:,0]*a.ball(scalar)*self.cb(local_phases[start:stop].conjugate())
            inds=site_indices[start:stop]
            cr[inds]=value.real.c;ci[inds]=value.imag.c;rr[inds]=value.real.r;ri[inds]=value.imag.r
        full=ComplexBall(RealBall(a,cr,rr),RealBall(a,ci,ri))
        complete=unpack(export['complete_site_amplitudes'],np)
        self.need(label+':all4096 original complete amplitudes retained',complete.shape==(4096,))
        reordered=np.zeros(4096,dtype=np.complex128)
        reordered[site_indices]=complete[actual_indices]*local_phases.conjugate()
        identity=(full-self.cb(reordered)).norm_upper()
        self.error(label+':independent all4096 tensor/export identity',identity,LIMITS['conversion'])
        normalized,norm2=self.normalize(full)
        self.error(label+':actual post-MPS raw norm',self.upper_abs(norm2-1),LIMITS['norm'])
        return normalized,norm2,[t.shape[1]for t in tensors[:-1]],identity,full

    def embed_sector(self,vector):
        np=self.np;real=np.zeros(4096);imag=real.copy()
        for index,(u,d)in enumerate(coordinates()):
            location=canonical_site_linear_index(u,d);value=regroup_sign(u,d)*vector[index]
            real[location]=value.real;imag[location]=value.imag
        return self.cb(real+1j*imag)

    def purity(self,normalized_full):
        a=self.a
        matrix=ComplexBall(RealBall(a,normalized_full.real.c.reshape((64,64),order='F'),normalized_full.real.r.reshape((64,64),order='F')),
            RealBall(a,normalized_full.imag.c.reshape((64,64),order='F'),normalized_full.imag.r.reshape((64,64),order='F')))
        rho=matrix.matmul(ComplexBall(matrix.real.T,-matrix.imag.T))
        # Exact normalized tensor outer product is Hermitian/PSD; Frobenius^2
        # equals Tr(rho^2). No symmetrization/eigenvalue clipping occurs.
        result=rho.norm_squared()
        if result.lower()<=0:raise RuntimeError('Fresh purity lower endpoint is nonpositive')
        return result

    def sample(self,row,backend,initial,final_model,k,arm):
        np=self.np;a=self.a;label=f'{backend}:k{k}:{arm}'
        self.checkpoint(label+':before_state_certificate')
        self.need(label+':exact original time and arm',type(row['time_numerator'])is int and row['time_numerator']==k and
            type(row['time_denominator'])is int and row['time_denominator']==8 and row['time']==k/8 and row['arm']==arm)
        budget=ARMS[arm]['native_maxdim'if backend=='native'else'tenpy_budget']
        if backend=='native':
            self.need(label+':registered native cutoff/bond/normalization call',row['maxdim']==budget and row['cutoff']==ARMS[arm]['native_cutoff']and
                row['actual_normalize_call']=='LinearAlgebra.normalize!(psi) before actual user entropy'and
                math.isfinite(row['pre_entropy_MPS_norm_before_normalize'])and row['pre_entropy_MPS_norm_before_normalize']>0)
        else:
            self.need(label+':registered actual from_full cutoff/external budget',row['external_bond_budget']==budget and
                row['from_full_cutoff']==ARMS[arm]['tenpy_cutoff'])
        raw=unpack(row['actual_raw_evolved_state'],np);used=unpack(row['actual_normalized_MPS_input_sector_state'],np)
        self.need(label+':all400 actual raw/normalized input states',raw.shape==used.shape==(400,))
        raw_ball=self.cb(raw);input_ball=self.cb(used);raw_norm2=raw_ball.norm_squared();norm=row['raw_norm']
        self.need(label+':finite positive original raw evolution norm',type(norm)in(int,float)and math.isfinite(norm)and norm>0)
        self.error(label+':original raw norm association',self.upper_abs(raw_norm2.sqrt()-norm),LIMITS['norm'])
        self.error(label+':actual raw evolved norm',self.upper_abs(raw_norm2-1),LIMITS['norm'])
        input_norm2=input_ball.norm_squared()
        self.error(label+':actual normalized MPS input norm',self.upper_abs(input_norm2-1),LIMITS['norm'])
        normalization_identity=(input_ball-raw_ball.divide_real(a.ball(norm))).norm_upper()
        self.error(label+':actual raw-to-input normalization formula',normalization_identity,LIMITS['conversion'])
        raw_to_used=(input_ball-raw_ball).norm_upper()
        spectral=spectral_action_ball(a,final_model['V'],final_model['values'],initial['vector'],k/8).reshape((400,))
        spectral_error=(raw_ball-spectral).norm_upper()
        self.error(label+':actual direct complete spectral action identity',spectral_error,LIMITS['spectral_identity'])
        op_error=eigensystem_operator_error(a,final_model['certificate'],k/8)
        op_term=(a.ball(op_error)*a.ball(initial['norm_upper'])).upper()
        post_normalized,norm2,bonds,export_error,post_raw=self.reconstruct_all_sites(row,final_model['layouts'],label)
        self.need(label+':actual all-bond rank observed within registered budgets',
            row['observed_max_bond']==max(bonds)and max(bonds)<=64 and max(bonds)<budget)
        if backend=='native':
            self.error(label+':post-normalize norm API association',self.upper_abs(norm2.sqrt()-row['pre_entropy_MPS_norm_after_normalize']),LIMITS['norm'])
            self.error(label+':actual user Schmidt normcheck association',self.upper_abs(norm2-row['library_normcheck']),LIMITS['norm'])
        embedded=self.embed_sector(used)
        conversion=(post_normalized-embedded).norm_upper()
        raw_conversion=(post_raw-embedded).norm_upper()
        post_normalization_change=(post_normalized-post_raw).norm_upper()
        self.error(label+':complete normalized pre/post MPS state distance',conversion,LIMITS['conversion'])
        outside=np.asarray([canonical_site_linear_index(u,d)for u in range(64)for d in range(64)
            if u.bit_count()!=3 or d.bit_count()!=3],dtype=np.int64)
        leakage=post_normalized[outside].norm_squared().upper()
        # The full4096-vector conversion bound includes this outside-sector
        # weight. It is recorded, never projected away or forced to exact zero.
        distance=(a.ball(initial['distance'])+a.ball(op_term)+a.ball(spectral_error)+a.ball(raw_to_used)+a.ball(conversion)).upper()
        purity=self.purity(post_normalized);own_entropy=-log_positive(a,purity)
        schmidt=np.asarray(row['post_call_schmidt_values_new_observer'],dtype=np.float64)
        self.need(label+':finite observed Schmidt values and physical rank',schmidt.ndim==1 and 1<=len(schmidt)<=64 and
            np.all(np.isfinite(schmidt))and np.all(schmidt>=0))
        p=a.ball(schmidt)*a.ball(schmidt);probability_sum=p.sum()
        if probability_sum.lower()<=0:raise RuntimeError('Observed Schmidt sum is nonpositive')
        schmidt_purity=(p*p).sum()/(probability_sum*probability_sum)
        self.error(label+':observed Schmidt sum normalization',self.upper_abs(probability_sum-1),LIMITS['norm'])
        self.error(label+':new observed Schmidt purity vs full tensor state',self.upper_abs(schmidt_purity-purity),LIMITS['library_entropy'])
        omission=a.ld(0);norm_allowance=a.ld(0);drop_purity=a.ld(0)
        if backend=='tenpy':
            self.error(label+':actual reported Schmidt probability sum',self.upper_abs(probability_sum-row['library_schmidt_probability_sum']),LIMITS['norm'])
            # Actual TeNPy threshold is its represented binary64 literal1e-30.
            # A normalized pure-state purity floor is not used to clip data.
            drop_purity=(a.ball(64)*a.ball(1e-30)*a.ball(1e-30)).upper()
            lower=purity.lower()
            if lower<=drop_purity:raise RuntimeError('Omitted probability purity budget has no positive denominator')
            omission=(a.ball(drop_purity)/(a.ball(lower)-a.ball(drop_purity))).upper()
            norm_allowance=(self.upper_abs(log_positive(a,probability_sum))*2)
        library=row['library_S2']
        self.need(label+':finite original actual library S2',type(library)in(int,float)and math.isfinite(library))
        library_error=(a.ball(self.upper_abs(a.ball(library)-own_entropy))+a.ball(omission)+a.ball(norm_allowance)).upper()
        self.error(label+':actual library entropy vs own normalized full state plus budgets',library_error,LIMITS['library_entropy'],'physical')
        physical_purity=RealBall(a,purity.c,a.inflate(purity.r+(a.ball(distance)*4).upper()))
        physical_entropy=entropy_from_purity_with_state_budget(a,purity,distance)
        library_radius=self.upper_abs(physical_entropy-a.ball(library))
        library_physical=RealBall(a,a.ld(library),a.inflate(library_radius))
        if arm=='primary':self.error(label+':primary total entropy radius',library_physical.r,LIMITS['primary_radius'])
        record={'backend':backend,'time_numerator':k,'time_denominator':8,'arm':arm,'library_S2_original':library,
            'raw_evolved_norm_squared':self.serialize_ball(raw_norm2),'raw_norm_reported_original':norm,
            'normalized_input_norm_squared':self.serialize_ball(input_norm2),'post_MPS_raw_norm_squared':self.serialize_ball(norm2),
            'actual_initial_norm_upper_multiplied_into_operator_error':self.bound(initial['norm_upper']),
            'initial_ground_and_column_normalization_distance_upper':self.bound(initial['distance']),
            'final_spectral_operator_error_upper':self.bound(op_error),'operator_error_times_actual_initial_norm_upper':self.bound(op_term),
            'actual_spectral_state_identity_error_upper':self.bound(spectral_error),'raw_to_used_normalization_change_upper':self.bound(raw_to_used),
            'normalization_formula_identity_error_upper':self.bound(normalization_identity),'complete_MPS_conversion_error_upper':self.bound(conversion),
            'raw_post_MPS_to_input_distance_diagnostic_upper':self.bound(raw_conversion),
            'raw_post_MPS_to_exact_normalized_state_distance_upper':self.bound(post_normalization_change),
            'normalization_budget_rule':'Physical conversion is measured directly to exact-normalized post MPS; raw normalization change is retained without double counting.',
            'complete_tensor_export_identity_error_upper':self.bound(export_error),'outside_balanced_sector_weight_upper_retained':self.bound(leakage),
            'total_physical_state_distance_upper':self.bound(distance),'bond_ranks_original':bonds,
            'Schmidt_probability_sum':self.serialize_ball(probability_sum),'TeNPy_omitted_purity_allowance':self.bound(drop_purity),
            'TeNPy_omitted_entropy_allowance':self.bound(omission),'TeNPy_Schmidt_normalization_entropy_allowance':self.bound(norm_allowance),
            'library_own_state_entropy_discrepancy_plus_budgets_upper':self.bound(library_error),
            'fresh_normalized_MPS_purity':self.serialize_ball(purity),'physical_purity':self.serialize_ball(physical_purity),
            'fresh_normalized_MPS_entropy':self.serialize_ball(own_entropy),'physical_entropy_centered_on_actual_library_value':self.serialize_ball(library_physical),
            'certificate_scope':'Conditional floating balls for actual full states; raw data unchanged, no directed intervals or eigenvector phase fitting.'}
        if backend=='native':record.update(pre_entropy_MPS_norm_before_normalize_original=row['pre_entropy_MPS_norm_before_normalize'],
            pre_entropy_MPS_norm_after_normalize_original=row['pre_entropy_MPS_norm_after_normalize'],actual_normalize_call=row['actual_normalize_call'])
        self.report['state_certificates'].append(record);self.checkpoint(label+':after_state_certificate')
        return library_physical,physical_purity

    def run(self,native,reference,pilot):
        np=self.np;a=self.a;times=range(5 if pilot else 33);raws={'native':native,'tenpy':reference};models={};initial={};curves={};purities={}
        schemas={'native':'P13_NATIVE_MPO_ENTROPY_DENSE_ADAPTER_CORE_V2','tenpy':'P13_TENPY_EXACTDIAG_ENTROPY_PRODUCER_CORE_V2'}
        for backend,raw in raws.items():
            registered_rows=validate_registration(raw,backend,pilot)
            self.need(backend+':source-defined complete producer scope',raw.get('schema')==schemas[backend]and
                raw.get('status')=='PRODUCER_COMPLETE_NOT_A_SCIENTIFIC_PASS'and raw.get('scientific_credit')==0 and
                type(raw.get('pilot_only'))is bool and raw['pilot_only']==pilot and raw['coordinates']==[list(x)for x in coordinates()])
            for stage in('before','after'):
                self.need(backend+':saved actual runtime '+stage,isinstance(raw.get('runtime_'+stage),dict)and raw['runtime_'+stage].get('status')=='PASS')
            if backend=='tenpy':self.need('tenpy:source-recorded initial/final physical mapping equality',raw.get('initial_final_physical_basis_mapping_equal')is True)
            models[backend]=self.qualify_models(raw,backend);initial[backend]=self.initial_certificate(raw,models[backend][0],backend)
            keys=[(row['time_numerator'],row['arm'])for row in raw['samples']]
            self.need(backend+':complete unique registered time-arm coverage',len(keys)==2*len(times)and len(set(keys))==len(keys)and
                set(keys)=={(k,arm)for k in times for arm in ARMS})
            by=registered_rows
            for k in times:
                left=by[k,'coarse'];right=by[k,'primary']
                self.need(backend+f':k{k}:same original pre-conversion evolution in both arms',
                    left['actual_raw_evolved_state']==right['actual_raw_evolved_state']and left['raw_norm']==right['raw_norm']and
                    left['actual_normalized_MPS_input_sector_state']==right['actual_normalized_MPS_input_sector_state'])
            for arm in ARMS:
                curves[backend,arm]=[];purities[backend,arm]=[]
                for k in times:
                    entropy,purity=self.sample(by[k,arm],backend,initial[backend],models[backend][1],k,arm)
                    curves[backend,arm].append(entropy);purities[backend,arm].append(purity)
        cross_s2=0;cross_purity=0;arm_checks=0;slope_checks=0
        for k in times:
            self.error(f'k{k}:complete primary S2 cross-backend',self.upper_abs(curves['native','primary'][k]-curves['tenpy','primary'][k]),LIMITS['cross_entropy'],'physical')
            cross_s2+=1
            self.error(f'k{k}:complete primary purity cross-backend',self.upper_abs(purities['native','primary'][k]-purities['tenpy','primary'][k]),LIMITS['cross_purity'],'physical')
            cross_purity+=1
            for backend in raws:
                self.error(f'{backend}:k{k}:primary/coarse S2 stability',self.upper_abs(curves[backend,'primary'][k]-curves[backend,'coarse'][k]),LIMITS['arm_entropy'],'physical')
                arm_checks+=1
        if pilot:
            self.report['growth_rate_status']='NOT_RUN_PILOT_FIVE_TIMES_DO_NOT_COVER_FIXED_NINE_POINT_WINDOW'
            self.report['nontrivial_full_curve_status']='NOT_RUN_PILOT_ONLY'
        else:
            slopes={}
            for key,values in curves.items():
                slope=growth_rate_ball(a,values);slopes[key]=slope
                self.report['growth_rates'].append({'backend':key[0],'arm':key[1],**self.serialize_ball(slope),
                    'time_numerators':list(range(4,13)),'weight_numerators':[-8,-6,-4,-2,0,2,4,6,8],
                    'weight_denominator':15,'scope':'Fixed nine-point OLS including intercept; finite-window rate, no universal/asymptotic claim.'})
            self.error('fixed nine-point primary growth slope cross-backend',self.upper_abs(slopes['native','primary']-slopes['tenpy','primary']),LIMITS['slope'],'physical')
            slope_checks+=1
            self.report['growth_rate_status']='COMPLETE_FIXED_NINE_POINT_OLS'
            for backend in raws:
                curve=curves[backend,'primary']
                lower=a.down(max(value.lower()for value in curve)-min(value.upper()for value in curve))
                self.gate(backend+':nontrivial registered full curve',lower>=LIMITS['curve_range'],'exercise',
                    range_lower_decimal=str(lower),required_lower=LIMITS['curve_range'])
            self.report['nontrivial_full_curve_status']='CHECKED_FIXED_FULL_CURVE_NO_REPLACEMENT'
        expected={'model_certificates':4,'initial_state_certificates':2,'state_certificates':4*len(times),
            'times_per_backend_arm':len(times),'primary_cross_backend_S2_checks':len(times),
            'primary_cross_backend_purity_checks':len(times),'primary_coarse_curve_checks':2*len(times),
            'growth_rate_certificates':0 if pilot else 4,'growth_window_points':0 if pilot else 9,
            'primary_growth_slope_checks':0 if pilot else 1}
        actual={**expected,'model_certificates':len(self.report['model_certificates']),
            'initial_state_certificates':len(self.report['initial_state_certificates']),
            'state_certificates':len(self.report['state_certificates']),'growth_rate_certificates':len(self.report['growth_rates']),
            'times_per_backend_arm':min(len(value)for value in curves.values()),
            'primary_cross_backend_S2_checks':cross_s2,'primary_cross_backend_purity_checks':cross_purity,
            'primary_coarse_curve_checks':arm_checks,'primary_growth_slope_checks':slope_checks,
            'growth_window_points':0 if pilot else min(len(row['time_numerators'])for row in self.report['growth_rates'])}
        self.report['coverage']=actual;self.report['expected_coverage']=expected
        self.need('full exact registered comparison coverage',actual==expected)
        self.report['comparison_run_completed']=True


def compare(np,native,reference,*,pilot=False,runtime_guard=None,on_record=None):
    """Called only by a new admitted wrapper, after package imports/qualification."""
    if not callable(runtime_guard)or not callable(on_record):
        raise RuntimeError('Actual allocated source/runtime and checkpoint callbacks are mandatory')
    if type(pilot)is not bool:raise TypeError('Explicit Boolean pilot/full mode required')
    report={'schema':'P13_FULL_STATE_ENTROPY_COMPARISON_CORE_V1','status':'UNRESOLVED','catalogue_id':'P13',
        'catalogue_status':'NOT_RUN','scientific_credit':0,'pilot_only':pilot,'comparison_run_completed':False,
        'conditional_floating_certificate_not_interval_arithmetic':True,'no_automatic_catalogue_credit':True,
        'gates':[],'model_certificates':[],'initial_state_certificates':[],'state_certificates':[],'growth_rates':[],
        'growth_rate_status':'NOT_RUN','nontrivial_full_curve_status':'NOT_RUN'}
    on_record('checkpoint',copy.deepcopy(report));a=None;old_policy=None
    try:
        report['runtime_before']=runtime_guard('before',api_witnesses(np))
        if not isinstance(report['runtime_before'],dict)or report['runtime_before'].get('status')!='PASS':raise RuntimeError('Actual comparator before runtime unresolved')
        old_policy=np.seterr(all='raise');a=Arithmetic(np)
        tool=Comparison(np,a,report,on_record);tool.run(native,reference,pilot)
    except Exception as error:
        report['error']={'type':type(error).__name__,'message':str(error)}
    finally:
        try:
            report['runtime_after']=runtime_guard('after',api_witnesses(np))
            if not isinstance(report['runtime_after'],dict)or report['runtime_after'].get('status')!='PASS':raise RuntimeError('Actual comparator after runtime unresolved')
        except Exception as error:report['after_guard_error']={'type':type(error).__name__,'message':str(error)}
        if old_policy is not None:np.seterr(**old_policy)
        if a is not None:
            report['arithmetic']=a.metadata;report['arithmetic_calls']=dict(a.calls)
        failed=[row for row in report['gates']if row['status']=='FAIL'];report['failed_gate_names']=[row['name']for row in failed]
        if report.get('error')or report.get('after_guard_error')or not report['comparison_run_completed']:
            report['status']='UNRESOLVED_PARTIAL_OR_INVALID_EVIDENCE'
        elif any(row['kind']in('precision','eligibility','exercise')for row in failed):
            report['status']='UNRESOLVED_PRECISION_ELIGIBILITY_OR_NONEXERCISING_FIXED_FIXTURE'
        elif failed:report['status']='FAIL_PHYSICAL_GATES_PENDING_INDEPENDENT_AUDIT'
        else:report['status']='PASS_PILOT_ZERO_CREDIT_PENDING_INDEPENDENT_AUDIT'if pilot else'PASS_FULL_PENDING_INDEPENDENT_AUDIT'
        on_record('checkpoint',copy.deepcopy(report))
    return report
