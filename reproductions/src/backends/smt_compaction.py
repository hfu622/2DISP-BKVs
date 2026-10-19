"""Exact integer-contour LRA optimizer for recognized continuous compaction MPS.

Only directly encoded integer supporting planes are accepted; no coefficient
or constant rounding is allowed. Movement windows remain binary64 rationals. No stored layout,
proof, or answer is loaded. Unknown never certifies optimality.
"""
import sys,json,time,os,math,hashlib,gzip
from fractions import Fraction as F
from functools import reduce
from pathlib import Path

class Unsupported(Exception):pass
def require(condition):
    if not condition:raise Unsupported()

def parse(payload):
    rows=[];section=None;integer=False;bounds={};kinds={};objectives={}
    for line in payload.decode('ascii').splitlines():
        p=line.split()
        if p[0] in ['NAME','ROWS','COLUMNS','RHS','BOUNDS','ENDATA']:section=p[0];continue
        if section=='ROWS':
            if p[1]=='OBJ':require(p==['N','OBJ']);continue
            require(int(p[1][1:])==len(rows));rows.append([[],{'L':'<','E':'=','G':'>'}[p[0]],0.])
        elif section=='COLUMNS':
            if p[1]=="'MARKER'":integer=p[2]=="'INTORG'";continue
            i=int(p[0][1:]);a=float(p[2]);kinds[i]=integer
            if p[1]=='OBJ':objectives[i]=a
            else:rows[int(p[1][1:])][0].append((i,a))
        elif section=='RHS':rows[int(p[1][1:])][2]=float(p[2])
        elif section=='BOUNDS':bounds.setdefault(int(p[2][1:]),{})[p[0]]=float(p[3])
        else:raise Unsupported()
    require(set(bounds)==set(range(len(bounds))))
    return rows,bounds,kinds,objectives

def solve(payload,hint,seconds):
    started=time.monotonic();deadline=started+seconds
    rows,bounds,kinds,objectives=parse(payload);n=len(bounds)
    base=next((i for i in range(n) if kinds[i]),n)
    require(base>=3 and base%2==1)
    require(all(not kinds[i] for i in range(base)) and all(kinds[i] for i in range(base,n)))
    require(objectives[0]==1 and all(objectives[i]==0 for i in range(1,n)))
    require(all(bounds[i]=={'LO':0.,'UP':1.} for i in range(base,n)))
    require(all(set(bounds[i])=={'LO','UP'} for i in range(base)))
    sys.path.insert(0,os.environ.get('BKV_Z3_PATH','F:/packing-deps/z3-bkv'))
    import z3
    x=[z3.Real('V'+str(i)) for i in range(base)]
    def exact(value):
        q=F.from_float(value);return z3.Q(q.numerator,q.denominator)
    constraints=[]
    for i in range(base):constraints.extend([x[i]>=exact(bounds[i]['LO']),x[i]<=exact(bounds[i]['UP'])])
    groups=[];facet={};max_coeff=0.;max_constant=0.
    def plane(terms,rhs):
        require(all(math.isfinite(a) and a.is_integer() and abs(a)<2**52 for i,a in terms))
        require(math.isfinite(rhs) and rhs.is_integer() and abs(rhs)<2**52)
        return z3.Sum([int(a)*x[i] for i,a in terms]),int(rhs)
    for terms,sense,rhs in rows:
        binary=[(i,a) for i,a in terms if i>=base];real=[(i,a) for i,a in terms if i<base]
        if not real:
            require(sense=='=' and rhs==1 and all(a==1 for i,a in binary));groups.append([i for i,a in binary]);continue
        if binary:
            require(len(binary)==1 and binary[0][1]>0 and sense=='<')
            z,m=binary[0];require(z not in facet)
            # Subtract exact rationals before converting, so a cancellation
            # cannot hide a noninteger or rounded big-M constant.
            require(m.is_integer() and abs(m)<2**52 and rhs.is_integer() and abs(rhs)<2**52)
            constant=F.from_float(rhs)-F.from_float(m);require(constant.denominator==1 and abs(constant)<2**52)
            expression,gamma=plane(real,float(constant));facet[z]=expression<=gamma
        else:
            expression,gamma=plane(real,rhs)
            constraints.append(expression<=gamma if sense=='<' else expression>=gamma if sense=='>' else expression==gamma)
    require(sorted(i for group in groups for i in group)==list(range(base,n)))
    for group in groups:constraints.append(z3.Or([facet.get(i,z3.BoolVal(True)) for i in group]))
    opt=z3.Optimize();opt.set(timeout=max(1,int((deadline-time.monotonic())*1000)));opt.add(*constraints)
    if hint is not None:
        require(len(hint)==n)
        for i,v in enumerate(x):opt.set_initial_value(v,exact(hint[i]))
    objective=opt.minimize(x[0]);status=opt.check()
    lower=opt.lower(objective);upper=opt.upper(objective)
    numeric=lambda q:z3.is_rational_value(q) or z3.is_int_value(q)
    optimal=status==z3.sat and numeric(lower) and numeric(upper) and z3.is_true(z3.simplify(lower==upper))
    out=dict(applicable=True,status=str(status),optimal=optimal,lower=str(lower),upper=str(upper),runtime=time.monotonic()-started,
        solver_version=z3.get_version_string(),coordinate_columns=base,disjunctions=len(groups),max_coefficient_recovery_error=max_coeff,
        max_constant_recovery_error=max_constant,model_sha256=hashlib.sha256(payload).hexdigest(),bounds='exact binary64 rationals',
        geometry_convention='Direct source-vertex integer supporting planes; no coefficient recovery or rounding; original encoded rows rechecked')
    try:
        model=opt.model();values=[0.]*n
        for i,v in enumerate(x):
            q=model.eval(v,model_completion=True);require(z3.is_rational_value(q));values[i]=q.numerator_as_long()/q.denominator_as_long()
        for group in groups:
            selected=next(i for i in group if z3.is_true(model.eval(facet.get(i,z3.BoolVal(True)),model_completion=True)));values[selected]=1.
        violations=[]
        for terms,sense,rhs in rows:
            lhs=math.fsum(a*values[i] for i,a in terms);violations.append(abs(lhs-rhs) if sense=='=' else max(0,lhs-rhs) if sense=='<' else max(0,rhs-lhs))
        bound_error=max(max(0,bounds[i]['LO']-v,v-bounds[i]['UP']) for i,v in enumerate(values))
        out.update(max_original_row_violation=max(violations),max_original_bound_violation=bound_error)
        if max(violations)<1e-7 and bound_error<1e-7:out.update(values=values,objective=values[0])
        else:out['optimal']=False
    except (z3.Z3Exception,Unsupported,StopIteration):out['optimal']=False
    return out

if __name__=='__main__':
    payload=Path(sys.argv[1]).read_bytes();hint=json.loads(Path(sys.argv[2]).read_text())
    try:out=solve(payload,hint,float(sys.argv[4]))
    except Unsupported:out={'applicable':False,'optimal':False,'status':'unsupported_model'}
    target=Path(sys.argv[3]);target.with_suffix('.tmp').write_text(json.dumps(out,allow_nan=False));os.replace(target.with_suffix('.tmp'),target)
    trace=os.environ.get('BKV_SMT_TRACE_DIR')
    if trace:
        folder=Path(trace);folder.mkdir(parents=True,exist_ok=True);key=hashlib.sha256(payload).hexdigest()
        (folder/(key+'.mps.gz')).write_bytes(gzip.compress(payload,mtime=0))
        (folder/(key+'.json')).write_text(json.dumps(out,indent=2)+'\n')
