#!/usr/bin/env python3
"""Independent prescribed-profile ABL development cases; not a reproduced experiment."""
import argparse, json, math, random, os
from pathlib import Path

ROOT=Path(os.environ['FOAM_RUN'])/'sgsDispersion'
def header(name,kind='dictionary'):
    return f'FoamFile\n{{ version 2.0; format ascii; class {kind}; object {name}; }}\n'
def write(case,name,text,kind='dictionary'):
    target=case/name; target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(header(target.name,kind)+text,encoding='utf-8')
def listfield(values):
    return f'nonuniform List<scalar>\n{len(values)}\n(\n'+ '\n'.join(f'{x:.12g}' for x in values)+'\n)'
def field(case,name,dims,values,inlet,ground):
    internal=listfield(values) if isinstance(values,list) else f'uniform {values}'
    write(case,'0/'+name,f'dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n'
          f'inlet {{ {inlet} }}\noutlet {{ type zeroGradient; }}\nground {{ {ground} }}\n'
          'top { type zeroGradient; }\nsideA { type cyclic; }\nsideB { type cyclic; }\n}\n',
          'volVectorField' if name=='U' else 'volScalarField')

def make(suite):
    smoke=suite=='smoke'
    nx,ny,nz=(12,6,6) if smoke else (24,12,12)
    lx,ly,lz=240.,120.,120.
    dt=0.25 if smoke else 0.5
    soi,end=(2.,10.) if smoke else (410.,440.)
    count=300 if smoke else 3000
    ustar,z0,kappa,kTotal=0.3,0.1,0.41,1.0
    averaging=lz/ustar
    dx,dy,dz=lx/nx,ly/ny,lz/nz
    delta=(dx*dy*dz)**(1/3)
    heights=[(k+.5)*dz for k in range(nz)]
    vel=lambda z:ustar/kappa*math.log((z+z0)/z0)
    eps=lambda z:ustar**3/(kappa*(z+z0))
    energy=lambda z:(eps(z)*delta/.93)**(2/3)
    rng=random.Random(20261011)
    points=[(40+rng.uniform(-.25,.25),rng.uniform(-.25,.25),50+rng.uniform(-.25,.25)) for _ in range(count)]
    names=('RANS_DRW','LES_resolved','LES_SGS')
    manifest={'suite':suite,'domain_m':[lx,ly,lz],'cells':[nx,ny,nz],
              'dt_s':dt,'release_s':soi,'end_s':end,'averaging_s':averaging,
              'uStar_m_s':ustar,'z0_m':z0,'kPrescribed_m2_s2':kTotal,
              'diameters_m':[20e-6,80e-6],'particleDensity_kg_m3':1000,
              'parcels_per_diameter':count,'release_center_m':[40,0,50],
              'RANS_carrier':'prescribed frozen U/k/epsilon profiles; solvePrimaryRegion false',
              'LES_carrier':'evolving kEqn LES; shared digital-filter inlet and one-way coupling',
              'model':'sgsDispersion', 'statistics':'internal model-owned exact time window',
              'functions_in_controlDict':False,
              'drag_model':'schillerNaumannDrag: paper equation (11), spherical particles',
              'purpose':'pulse dispersion / model integration; not experimental validation',
              'mean_drw_k_over_epsilon_at_release_s':kTotal/eps(50),
              'initial_les_sgs_k_at_release':energy(50)}
    dest=ROOT/suite
    for name in names:
        case=dest/name
        if (case/'system/controlDict').exists():
            raise RuntimeError(f'Refusing to overwrite existing case: {case}')
        les=name.startswith('LES'); corrected=name=='LES_SGS'
        dispersion='sgsDispersion' if corrected else ('none' if les else 'stochasticDispersionRAS')
        write(case,'system/blockMeshDict',f'''scale 1;
vertices ((0 -60 0) (240 -60 0) (240 60 0) (0 60 0)
          (0 -60 120) (240 -60 120) (240 60 120) (0 60 120));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1));
edges ();
boundary
(
 inlet {{ type patch; faces ((0 4 7 3)); }}
 outlet {{ type patch; faces ((1 2 6 5)); }}
 ground {{ type wall; faces ((0 3 2 1)); }}
 top {{ type patch; faces ((4 5 6 7)); }}
 sideA {{ type cyclic; neighbourPatch sideB; faces ((0 1 5 4)); }}
 sideB {{ type cyclic; neighbourPatch sideA; faces ((3 7 6 2)); }}
);
mergePatchPairs ();
''')
        functions=''
        write(case,'system/controlDict',f'''application kinematicParcelFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end};
deltaT {dt}; adjustTimeStep false;
writeControl runTime; writeInterval {1 if smoke else 5};
purgeWrite 0; writeFormat ascii; writePrecision 10; writeCompression off;
timeFormat general; timePrecision 9; runTimeModifiable true;
libs ("libsgsDispersion.so" "libatmosphericModels.so");
{functions}''')
        write(case,'system/fvSchemes','''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
 default none;
 div(phi,U) Gauss linear;
 div(phi,k) Gauss upwind;
 div(phi,epsilon) Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''')
        write(case,'system/fvSolution',f'''solvers
{{
 p {{ solver GAMG; tolerance 1e-7; relTol 0.01; smoother GaussSeidel; }}
 pFinal {{ $p; relTol 0; }}
 "(U|k|epsilon)" {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.01; }}
 "(U|k|epsilon)Final" {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }}
}}
PIMPLE {{ solvePrimaryRegion {'true' if les else 'false'};
 momentumPredictor true; nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0; }}
relaxationFactors {{ equations {{ ".*" 1; }} }}
''')
        write(case,'constant/dynamicMeshDict','dynamicFvMesh staticFvMesh;\n')
        write(case,'constant/transportProperties','transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] 1.5e-5;\nrhoInf [1 -3 0 0 0 0 0] 1.2;\n')
        write(case,'constant/g','dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n','uniformDimensionedVectorField')
        write(case,'constant/turbulenceProperties',
              'simulationType LES;\nLES { LESModel kEqn; turbulence on; printCoeffs on; delta cubeRootVol; kEqnCoeffs { Ce 0.93; Ck 0.094; } }\n' if les else
              'simulationType RAS;\nRAS { RASModel kEpsilon; turbulence on; printCoeffs on; }\n')
        uvalues=[f'({vel(z):.12g} 0 0)' for z in heights for _ in range(nx*ny)]
        uinternal=f'nonuniform List<vector>\n{len(uvalues)}\n(\n'+'\n'.join(uvalues)+'\n)'
        digital=f'''type turbulentDigitalFilterInlet; fsm true; seed 20261011; n ({ny} {nz});
L (5 5 5 20 20 20 10 10 10);
mean {{ type mappedFile; mapMethod nearest; fieldTable UMean; }}
R {{ type mappedFile; mapMethod nearest; fieldTable RResolved; }}
value uniform ({vel(50):.12g} 0 0);'''
        boundaryvals=[f'({vel(z):.12g} 0 0)' for z in heights for _ in range(ny)]
        ub='fixedValue; value nonuniform List<vector>\n'+str(len(boundaryvals))+'\n(\n'+'\n'.join(boundaryvals)+'\n);'
        write(case,'0/U',f'dimensions [0 1 -1 0 0 0 0];\ninternalField {uinternal};\nboundaryField\n{{\n'
              f'inlet {{ {digital if les else "type "+ub} }}\n'
              'outlet { type zeroGradient; }\nground { type noSlip; }\ntop { type slip; }\n'
              'sideA { type cyclic; }\nsideB { type cyclic; }\n}\n','volVectorField')
        write(case,'0/p','dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n'
              'inlet { type zeroGradient; }\noutlet { type fixedValue; value uniform 0; }\n'
              'ground { type zeroGradient; }\ntop { type zeroGradient; }\n'
              'sideA { type cyclic; }\nsideB { type cyclic; }\n}\n','volScalarField')
        ks=[energy(z) if les else kTotal for z in heights for _ in range(nx*ny)]
        kinlet=[energy(z) if les else kTotal for z in heights for _ in range(ny)]
        field(case,'k','[0 2 -2 0 0 0 0]',ks,'type fixedValue; value '+listfield(kinlet)+';',
              'type zeroGradient;' if les else 'type kqRWallFunction; value uniform 1;')
        field(case,'nut','[0 2 -1 0 0 0 0]',0,'type calculated; value uniform 0;',
              'type atmNutkWallFunction; z0 uniform 0.1; value uniform 0;')
        if not les:
            ev=[eps(z) for z in heights for _ in range(nx*ny)]
            ein=[eps(z) for z in heights for _ in range(ny)]
            field(case,'epsilon','[0 2 -3 0 0 0 0]',ev,'type fixedValue; value '+listfield(ein)+';',
                  'type epsilonWallFunction; value uniform 0.01;')
        if les:
            patchpoints=[(0,-ly/2+(j+.5)*dy,z) for z in heights for j in range(ny)]
            pd=case/'constant/boundaryData/inlet'; (pd/'0').mkdir(parents=True,exist_ok=True)
            (pd/'points').write_text(header('points','vectorField')+str(len(patchpoints))+'\n(\n'+
                                   '\n'.join(f'({x} {y} {z})' for x,y,z in patchpoints)+'\n)\n',encoding='utf-8')
            (pd/'0/UMean').write_text(header('UMean','vectorField')+str(len(patchpoints))+'\n(\n'+
                                    '\n'.join(f'({vel(z):.12g} 0 0)' for x,y,z in patchpoints)+'\n)\n',encoding='utf-8')
            stress=[]
            for x,y,z in patchpoints:
                s2=2*energy(z)/3
                xx,yy,zz=1.1*kTotal-s2,.5*kTotal-s2,.4*kTotal-s2
                if min(xx,yy,zz)<=0 or xx*zz<=ustar**4: raise ValueError('Resolved inlet stress is not positive definite')
                stress.append(f'({xx:.12g} 0 {-ustar**2:.12g} {yy:.12g} 0 {zz:.12g})')
            (pd/'0/RResolved').write_text(header('RResolved','symmTensorField')+str(len(stress))+'\n(\n'+
                                        '\n'.join(stress)+'\n)\n',encoding='utf-8')
        write(case,'constant/releasePositions',str(len(points))+'\n(\n'+'\n'.join(f'({x:.12g} {y:.12g} {z:.12g})' for x,y,z in points)+'\n)\n','vectorField')
        injections=''
        for label,diameter in [('fine',20e-6),('coarse',80e-6)]:
            injections+=f'''{label}
{{ type manualInjection; SOI {soi}; massTotal 1e-9; parcelBasisType fixed; nParticle 1;
 positionsFile releasePositions; U0 ({vel(50):.12g} 0 0);
 sizeDistribution {{ type fixedValue; fixedValueDistribution {{ value {diameter}; }} }}
}}
'''
        sgsCoeffs = 'sgsDispersionCoeffs\n {\n     averagingTime 400;\n     C0 4;\n     beta 2;\n     cEpsilon 0.93;\n     verticalDirection (0 0 1);\n     stationaryInitial false;\n }' if corrected else ''
        write(case,'constant/kinematicCloudProperties',f'''solution
{{ active true; coupled false; transient true; cellValueSourceCorrection false; maxCo 0.3;
 interpolationSchemes {{ rho cell; U cellPoint; muc cell; }}
 integrationSchemes {{ U analytical; }}
}}
constantProperties {{ rho0 1000; }}
subModels
{{
 particleForces {{ schillerNaumannDrag; gravity; }}
 injectionModels {{ {injections} }}
 dispersionModel {dispersion};
 {sgsCoeffs}
 patchInteractionModel localInteraction;
 localInteractionCoeffs
 {{ patches (ground {{ type stick; }} "(inlet|outlet|top)" {{ type escape; }}); }}
 stochasticCollisionModel none;
 surfaceFilmModel none;
}}
cloudFunctions {{}}
''')
        (case/(name+'.foam')).touch()
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('suite',choices=['smoke','abl'])
    make(parser.parse_args().suite)
