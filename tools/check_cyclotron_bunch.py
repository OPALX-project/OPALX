#!/usr/bin/env python3
"""Small accelerating-bunch integration checks; no regression rebaselining.

Use ~/.venv-h6/bin/python. Maintained Cyclotron-RF-Stop supplies maps/launch.
This is a fixed-time code-path check, not a mesh-convergence or extraction study.
"""
from pathlib import Path
import argparse,hashlib,json,os,re,shutil,subprocess
import h5py,numpy as np

def run(binary,fixture,target,solver,charge,ranks):
    """Run the same seeded 3D bunch through the first gap, returning sorted states."""
    shutil.copytree(fixture,target)
    rng=np.random.default_rng(615)
    values=rng.normal(size=(256,6))
    values[:,::2]*=.001;values[:,1::2]*=.0001
    values[:,1::2]+=np.array([-.36920124,0,-.151830777])
    # Exact centering gives each backend/rank count the same reference launch.
    values[:,::2]-=values[:,::2].mean(axis=0)
    values[:,1::2]-=values[:,1::2].mean(axis=0)-np.array([-.36920124,0,-.151830777])
    with (target/'bunch.dat').open('w') as f:
        f.write('256\nx px y py z pz\n');np.savetxt(f,values,fmt='%.17e')
    source=target/'Cyclotron-RF-Stop.in';text=source.read_text()
    text=text.replace('PSDUMPFREQ=0','PSDUMPFREQ=1000000000').replace('FNAME="proton.dat", NPARTDIST=1','FNAME="bunch.dat", NPARTDIST=256')
    text=text.replace('NALLOC=1, BCHARGE=1.602176634e-19',f'NALLOC=256, BCHARGE={charge:.17g}')
    text=re.sub(r'MAXSTEPS=\d+','MAXSTEPS=32',text);text=re.sub(r', EKINSTOP=[^;]+','',text)
    text=text.replace('DT=4.1131951299769659e-11','DT=1.6452780519907864e-10')
    text=text.replace('TYPE=NONE, NX=16, NY=16, NZ=16',f'TYPE={solver}, NX=32, NY=32, NZ=32')
    if solver=='OPEN':
        text=text.replace('FS0: FIELDSOLVER,','SCBIN: BINNING, MAXBINS=1, PARAMETER=GAMMAZ, ADAPTIVEBINNING=FALSE;\nFS0: FIELDSOLVER, BINS=SCBIN,')
        text=text.replace('BCFFTZ=OPEN;','BCFFTZ=OPEN, BBOXINCR=2, GREENSF=STANDARD;')
    source.write_text(text)
    command=['mpirun','-np',str(ranks),str(binary),source.name,'--info','2']
    with (target/'run.log').open('w') as log:
        subprocess.run(command,cwd=target,env=dict(os.environ,OMP_NUM_THREADS='1',OMP_PROC_BIND='false'),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
    with h5py.File(target/'Cyclotron-RF-Stop.h5') as f:
        g=f[max(f,key=lambda k:int(k.split('#')[1]))];order=np.argsort(g['id'][:])
        states={n:g[n][:][order] for n in ['id','x','y','z','px','py','pz']}
        assert len(states['id'])==256
        assert all(np.isfinite(v).all() for v in states.values())
    return states

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--regression-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();binary=args.executable.resolve();args.output.mkdir(parents=True,exist_ok=False)
    fixture=args.regression_root/'RegressionTests/Cyclotron-RF-Stop'
    none=run(binary,fixture,args.output/'none','NONE',0,1)
    # The PIC adapter intentionally rejects exactly zero charge. Use a
    # negligible positive bunch charge; exact zero fields are tested in C++.
    zero=run(binary,fixture,args.output/'near-zero','OPEN',1e-30,1)
    one=run(binary,fixture,args.output/'charged-one','OPEN',3.9486673247778873e-11,1)
    two=run(binary,fixture,args.output/'charged-two','OPEN',3.9486673247778873e-11,2)
    errors={}
    # FROMFILE assigns IDs after rank distribution, so sample-to-ID mapping can
    # differ with rank count. Compare physical state sets, using unique x here.
    for states in (one,two):
        assert len(np.unique(states['id']))==256
        order=np.argsort(states['x'])
        for n in states:states[n]=states[n][order]
    for label,a,b in [('zero_charge',zero,none),('mpi',one,two)]:
        errors[label]={n:float(np.max(np.abs(a[n]-b[n]))) for n in a if n != "id"}
        for n in ['x','y','z','px','py','pz']:
            # Same fixed short physical run: 0.1 nm positions and 1e-10 p/(mc)
            # bound frame roundoff/reductions, not a relaxed physics reference.
            np.testing.assert_allclose(a[n],b[n],rtol=0,atol=1e-10)
    effect=max(float(np.max(np.abs(np.sort(one[n])-np.sort(zero[n])))) for n in ['px','py','pz'])
    assert effect>1e-8,'Charged run did not generate a resolved collective impulse'
    report={'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'maximum_errors':errors,'resolved_self_field_momentum_change':effect}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
