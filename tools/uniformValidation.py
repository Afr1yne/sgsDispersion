#!/usr/bin/env python3
"""Native cloud benchmark with prescribed linear shear and analytic transverse dispersion.
Uses stock kEqn/kinematicParcelFoam: no digital-filter patch or function objects.
"""
import json, math, os, shutil, subprocess, sys
from pathlib import Path
from makeCases import write, listfield, header
from analyse import read_vtp, foam_values, moment

root=Path(os.environ['FOAM_RUN'])/'sgsDispersion/validation/uniform'
G=.1; k=.1; delta=10.; Ck=.094; Ce=.93
nut=Ck*math.sqrt(k)*delta
eps=nut*G*G
sgs=(eps*delta/Ce)**(2/3)
variance=2*sgs/3
tau=2*variance/(4*eps)
count=30000
T=10.; release=2.; end=12.

def field(case,name,dims,internal,patches,kind):
    write(case,'0/'+name,f'dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n'+patches+'\n}\n',kind)

def make(case,h):
    if (case/'system/controlDict').exists(): return
    write(case,'system/blockMeshDict','''scale 1;
vertices ((0 -50 0)(200 -50 0)(200 50 0)(0 50 0)(0 -50 100)(200 -50 100)(200 50 100)(0 50 100));
blocks (hex (0 1 2 3 4 5 6 7) (20 10 10) simpleGrading (1 1 1));
edges ();
boundary
(
 inlet { type patch; faces ((0 4 7 3)); }
 outlet { type patch; faces ((1 2 6 5)); }
 ground { type wall; faces ((0 3 2 1)); }
 top { type patch; faces ((4 5 6 7)); }
 sideA { type cyclic; neighbourPatch sideB; faces ((0 1 5 4)); }
 sideB { type cyclic; neighbourPatch sideA; faces ((3 7 6 2)); }
);
mergePatchPairs ();
''')
    write(case,'system/controlDict',f'''application kinematicParcelFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end};
deltaT {h}; adjustTimeStep false;
writeControl timeStep; writeInterval {round(end/h)};
writeFormat ascii; writePrecision 13; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
libs ("libsgsDispersion.so");
''')
    write(case,'system/fvSchemes','''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss linear; div(phi,k) Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''')
    write(case,'system/fvSolution','''solvers { p { solver GAMG; tolerance 1e-9; relTol 0; smoother GaussSeidel; } U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0; } k { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0; } }
PIMPLE { solvePrimaryRegion false; nOuterCorrectors 1; nCorrectors 1; nNonOrthogonalCorrectors 0; }
''')
    write(case,'constant/dynamicMeshDict','dynamicFvMesh staticFvMesh;')
    write(case,'constant/transportProperties','transportModel Newtonian; nu [0 2 -1 0 0 0 0] 1.5e-5; rhoInf [1 -3 0 0 0 0 0] 1.2;')
    write(case,'constant/turbulenceProperties','simulationType LES; LES { LESModel kEqn; turbulence on; printCoeffs on; delta cubeRootVol; kEqnCoeffs { Ce 0.93; Ck 0.094; } }')
    write(case,'constant/g','dimensions [0 1 -2 0 0 0 0]; value (0 0 0);','uniformDimensionedVectorField')
    u=[f'({G*(5+10*z):.15g} 0 0)' for z in range(10) for _ in range(200)]
    inlet=[f'({G*(5+10*z):.15g} 0 0)' for z in range(10) for _ in range(10)]
    internal='nonuniform List<vector>\n2000\n(\n'+'\n'.join(u)+'\n)'
    ub='inlet { type fixedValue; value nonuniform List<vector>\n100\n(\n'+'\n'.join(inlet)+'\n); }\n'
    ub+='outlet { type zeroGradient; } ground { type fixedValue; value uniform (0 0 0); } top { type fixedValue; value uniform (10 0 0); } sideA { type cyclic; } sideB { type cyclic; }'
    field(case,'U','[0 1 -1 0 0 0 0]',internal,ub,'volVectorField')
    field(case,'p','[0 2 -2 0 0 0 0]','uniform 0','inlet { type zeroGradient; } outlet { type fixedValue; value uniform 0; } ground { type zeroGradient; } top { type zeroGradient; } sideA { type cyclic; } sideB { type cyclic; }','volScalarField')
    for name,value,dims in [('k',k,'[0 2 -2 0 0 0 0]'),('nut',nut,'[0 2 -1 0 0 0 0]')]:
        kind='fixedValue' if name=='k' else 'calculated'
        patches='\n'.join(f'{p} {{ type {kind}; value uniform {value:.15g}; }}' for p in ['inlet','outlet','ground','top'])+'\nsideA { type cyclic; } sideB { type cyclic; }'
        field(case,name,dims,f'uniform {value:.15g}',patches,'volScalarField')
    write(case,'constant/releasePositions',f'{count}\n(\n'+'(40 0 50)\n'*count+')','vectorField')
    write(case,'constant/kinematicCloudProperties',f'''solution
{{ active true; coupled false; transient true; cellValueSourceCorrection false; maxCo 0.3;
 interpolationSchemes {{ rho cell; U cellPoint; muc cell; }} integrationSchemes {{ U analytical; }} }}
constantProperties {{ rho0 1000; }}
subModels
{{ particleForces {{ schillerNaumannDrag; }}
 injectionModels {{ pulse {{ type manualInjection; SOI {release}; massTotal 1e-9; parcelBasisType fixed; nParticle 1; positionsFile releasePositions; U0 (5 0 0);
 sizeDistribution {{ type fixedValue; fixedValueDistribution {{ value 1e-6; }} }} }} }}
 dispersionModel sgsDispersion;
 sgsDispersionCoeffs {{ averagingTime 0.5; C0 4; beta 0; cEpsilon 0.93; stationaryInitial true; }}
 patchInteractionModel localInteraction;
 localInteractionCoeffs {{ patches (ground {{ type stick; }} "(inlet|outlet|top)" {{ type escape; }}); }}
 stochasticCollisionModel none; surfaceFilmModel none;
}}
cloudFunctions {{}}
''')

def run():
    records=[]
    for h in [.2,.1,.05]:
        case=root/('dt_'+str(h).replace('.','p'))
        make(case,h)
        for command in ['blockMesh','checkMesh','kinematicParcelFoam']:
            print(f'{case.name}: {command}',flush=True)
            with (case/('log.'+command)).open('w') as log:
                subprocess.run([command,'-case',str(case)],stdout=log,stderr=subprocess.STDOUT,check=True)
        assert 'Mesh OK' in (case/'log.checkMesh').read_text()
        assert '\nEnd\n' in (case/'log.kinematicParcelFoam').read_text()
        with (case/'log.foamToVTK').open('w') as log:
            subprocess.run(['foamToVTK','-case',str(case),'-ascii','-no-boundary','-no-internal','-time','12'],stdout=log,stderr=subprocess.STDOUT,check=True)
        directory=case/'VTK/lagrangian/kinematicCloud'
        series=json.loads((directory/'kinematicCloud.vtp.series').read_text())
        data=read_vtp(directory/series['files'][-1]['name'])
        assert len(data['Points'])==count and all(v>0 for v in data['active'])
        _,spread=moment([p[1] for p in data['Points']])
        modelVariance=foam_values(case/'12/kinematicCloudSGSVariance')
        modelEpsilon=foam_values(case/'12/kinematicCloudSGSDissipation')
        assert max(abs(v/variance-1) for v in modelVariance)<1e-9
        assert max(abs(v/eps-1) for v in modelEpsilon)<1e-9
        assert max(abs(v-1) for v in foam_values(case/'12/kinematicCloudSGSFraction'))<1e-10
        # Native injection updates yield age=T+h; first positional step has initial
        # zero transverse particle velocity, so stationary transport lasts age-h=T.
        age=sum(data['age'])/count
        duration=age-h
        reference=2*variance*tau*(duration-tau*(1-math.exp(-duration/tau)))
        measured=spread*spread
        error=measured/reference-1
        assert abs(error)<.035, (h,measured,reference,error)
        records.append({'dt_s':h,'particle_count':count,'native_age_s':age,
          'stationary_transport_duration_s':duration,'measured_y_variance_m2':measured,
          'Taylor_reference_m2':reference,'relative_error':error,
          'SGS_variance_m2_s2':variance,'dissipation_m2_s3':eps,'correlation_time_s':tau})
        print(json.dumps(records[-1]),flush=True)
    destination=root.parents[1]/'results/validation'
    destination.mkdir(parents=True,exist_ok=True)
    (destination/'native_uniform.json').write_text(json.dumps({'status':'PASS','scope':'native small-inertia limit, beta=0, prescribed linear shear; transverse homogeneous dispersion','records':records},indent=2))
    print('PASS native production, variance/fraction and transverse Taylor dispersion at three time steps',flush=True)

if __name__=='__main__':run()
