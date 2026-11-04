#!/usr/bin/env python3
"""Exactly 100 native TEAPOT-body/SC steps over ten metres at R=2 m."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

import numpy as np
import pandas as pd
from orbit.core.bunch import Bunch, BunchTwissAnalysis, SynchPartRedefinitionZdE
from orbit.core import teapot_base
from orbit.core.spacecharge import Grid1D, LSpaceChargeCalc, SpaceChargeCalcSliceBySlice2D
from check_bend_kinematics import body

E, C, MASS = 1.602176634e-19, 299792458., .938272029


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-json', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--zero-sc', action='store_true')
    parser.add_argument('--width-x', type=float, default=.256)
    parser.add_argument('--width-y', type=float, default=.320)
    parser.add_argument('--periodic-check', type=Path)
    parser.add_argument('--reference-path', type=Path, default=Path('/Users/adelmann/git/OPALX-2D5/curved-study/straight-round-paper/curved/opalx/sources/circle-R2-step0.005.dat'))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    assert not (args.out/'manifest.json').exists()
    started = time.monotonic()
    source = json.loads(args.source_json.read_text())
    assert source['geometry'] == 'bend' and source['source_particles'] == 1048576
    assert sha(source['npz']) == source['npz_sha256']
    with np.load(source['npz']) as data:
        initial = data['material_m']
        ids = data['source_particle_id']
        initial_p = data['momenta_betagamma']
    count = len(initial)
    assert np.array_equal(ids, np.arange(count))
    beta = source['beta']
    gamma = 1/math.sqrt(1-beta*beta)
    radius = source['bend_radius_m']
    assert radius == 2.
    path_xyz = np.loadtxt(args.reference_path)[:, [1, 2, 3]]
    chords = np.linalg.norm(np.diff(path_xyz, axis=0), axis=1)
    period = float(chords.sum())
    sagitta = float(radius-np.sqrt(radius*radius-(chords.max()/2)**2))
    total_length, steps, nx, ny, nz, physical_g = 10., 100, 512, 512, 512, 6.36
    width_x, width_y = args.width_x, args.width_y
    ds = total_length/steps
    dt = ds/(beta*C)
    native_g = physical_g*ds/period
    ratio = math.exp((native_g-1)/2)
    macro_size = source['macrocharge_C']/E
    bunch = Bunch()
    bunch.mass(MASS)
    bunch.charge(1.)
    bunch.macroSize(macro_size)
    bunch.addPartAttr('macrosize')
    kinetic = MASS*beta*beta*gamma*gamma/(gamma+1)
    bunch.getSyncParticle().kinEnergy(kinetic)
    assert np.max(np.abs(np.linalg.norm(initial_p, axis=1)-beta*gamma)) < 1e-15
    for x, y, z in initial:
        bunch.addParticle(float(x), 0., float(y), 0., float(z), 0.)
        bunch.partAttrValue('macrosize', bunch.getSize()-1, 0, macro_size)
    half_x, half_y = width_x/2*(1-1/nx), width_y/2*(1-1/ny)
    guards = np.array([[-half_x, 0., -period/2], [half_x, 0., period/2],
                       [0., -half_y, -period/2], [0., half_y, period/2],
                       [0., 0., -period/2], [0., 0., period/2]])
    for x, y, z in guards:
        bunch.addParticle(float(x), 0., float(y), 0., float(z), 0.)
        bunch.partAttrValue('macrosize', bunch.getSize()-1, 0, 0.)
    def reset_guards(reference_s=0.):
        for k, (x, y, z) in enumerate(guards):
            i = count+k
            bunch.x(i, float(x)); bunch.xp(i, 0.)
            bunch.y(i, float(y)); bunch.yp(i, 0.)
            bunch.z(i, float(z-reference_s)); bunch.dE(i, 0.)
    # Exact native deposition/gather test for the guard contribution.
    test = Bunch()
    for x, y, z in guards:
        test.addParticle(float(x), 0., float(y), 0., float(z), 0.)
    grid = Grid1D(nz, period)
    grid.setGridZ(-period/2, period/2)
    grid.setZero()
    grid.binBunchSmoothedByParticle(test)
    checked_z = np.linspace(-period/2+4*period/nz, period/2-4*period/nz, 4097)
    guard_gradient = np.array([grid.calcGradient(float(z)) for z in checked_z])
    assert np.max(np.abs(guard_gradient)) == 0.
    guard_check = dict(passed=True, method='native Grid1D smoothed unit-count deposition and native finite-difference gradient',
                       rolled_coordinate_interval_m=[float(checked_z.min()), float(checked_z.max())], samples=len(checked_z),
                       maximum_native_gradient=0., all_guards_at_lab_ends=True)
    (args.out/'guard-check.json').write_text(json.dumps(guard_check, indent=2)+'\n')
    if not args.zero_sc:
        assert args.periodic_check is not None
        periodic_result = json.loads(args.periodic_check.read_text())
        assert periodic_result['passed'] and periodic_result['selected_strategy'] == 'integer-cell rolling origin'
        solver = SpaceChargeCalcSliceBySlice2D(nx, ny, nz, half_x/half_y)
        solver.longTracking(0)
        longitudinal = LSpaceChargeCalc(ratio, period, 0, 0, nz)
        longitudinal.setUseGrad(1)
        longitudinal.setSmoothGrad(0)
    shift = SynchPartRedefinitionZdE()
    analysis = BunchTwissAnalysis()
    rows, mesh_checks = [], []
    def record(step):
        analysis.analyzeBunch(bunch)
        mean = np.array([analysis.getAverage(i) for i in range(6)])
        covariance = np.array([[analysis.getCorrelation(i, j) for j in range(6)] for i in range(6)])
        charge = analysis.getGlobalMacrosize()*E
        assert math.isclose(charge, source['charge_C'], rel_tol=2e-10)
        assert np.isfinite(mean).all() and np.isfinite(covariance).all()
        local_positions = np.array([[bunch.x(i), bunch.y(i), bunch.z(i)+step*ds] for i in range(count)])
        theta = local_positions[:, 2]/radius
        coordinates = np.c_[local_positions[:, 0]*np.cos(theta)-2*radius*np.sin(theta/2)**2,
                            local_positions[:, 1], (radius+local_positions[:, 0])*np.sin(theta)]
        global_mean = coordinates.mean(axis=0)
        global_cov = np.cov(coordinates, rowvar=False, ddof=0)
        row = dict(step=step, time_s=step*dt, nominal_reference_s_m=step*ds,
                   nominal_design_advance_m=step*ds, particles=count,
                   reference_phase_rad=step*ds/radius,
                   sync_time_s=bunch.getSyncParticle().time(), charge_C=charge,
                   source_particles=count, space_charge_enabled=not args.zero_sc,
                   sigma_global_x_m=math.sqrt(max(0., global_cov[0, 0])),
                   sigma_global_y_m=math.sqrt(max(0., global_cov[1, 1])),
                   sigma_global_z_m=math.sqrt(max(0., global_cov[2, 2])),
                   sigma_u_m=math.sqrt(max(0., covariance[0, 0])),
                   sigma_y_m=math.sqrt(max(0., covariance[2, 2])),
                   sigma_s_m=math.sqrt(max(0., covariance[4, 4])),
                   sigma_frenet_u_m=math.sqrt(max(0., covariance[0, 0])),
                   sigma_frenet_y_m=math.sqrt(max(0., covariance[2, 2])),
                   sigma_frenet_s_m=math.sqrt(max(0., covariance[4, 4])),
                   mean_global_x_m=global_mean[0], mean_global_y_m=global_mean[1], mean_global_z_m=global_mean[2],
                   mean_frenet_u_m=mean[0], mean_frenet_y_m=mean[2], mean_frenet_s_m=mean[4],
                   mean_u_m=mean[0], mean_y_m=mean[2], mean_s_m=mean[4], mean_absolute_s_m=mean[4]+step*ds)
        for i, name in enumerate(['u', 'y', 's']):
            for j, other in enumerate(['u', 'y', 's']):
                row[f'cov_frenet_{name}{other}_m2'] = covariance[[0, 2, 4][i], [0, 2, 4][j]]
                row[f'cov_global_{"xyz"[i]}{"xyz"[j]}_m2'] = global_cov[i, j]
        for i in range(6):
            row[f'native_mean_{i}'] = mean[i]
            for j in range(6):
                row[f'native_cov_{i}{j}'] = covariance[i, j]
        rows.append(row)
        pd.DataFrame(rows).to_csv(args.out/'moments.csv', index=False, float_format='%.17g')
        print(json.dumps({key: row[key] for key in ['step', 'sigma_u_m', 'sigma_y_m', 'sigma_s_m']}), flush=True)
    record(0)
    for step in range(steps):
        body(bunch, ds/2)
        midpoint = (step+.5)*ds
        if not args.zero_sc:
            # Exact integer-cell relabelling of the same periodic LAB mesh.
            # The physical bunch stays interior, leaving numerical markers in
            # vacuum. Selected native LS deposition/gather are already periodic.
            dz = period/nz
            origin_cells = int(np.rint(midpoint/dz))
            origin = origin_cells*dz
            temporary_shift = midpoint-origin
            label_indices = np.unique(np.r_[0, count-1, np.linspace(0, count-1, 256).astype(int)])
            labels_before = np.array([bunch.z(int(i)) for i in label_indices])
            shift.shift_Z(bunch, temporary_shift)
            reset_guards(0.)
            physical = np.array([[bunch.x(i), bunch.y(i), bunch.z(i)] for i in range(count)])
            physical_z = physical[:, 2]
            assert np.max(np.abs(physical[:, 0])) < half_x-2*width_x/nx
            assert np.max(np.abs(physical[:, 1])) < half_y-2*width_y/ny
            assert physical_z.min() > -period/2+4*dz and physical_z.max() < period/2-4*dz
            solver.trackBunch(bunch, ds)
            longitudinal.trackBunch(bunch)
            rho = solver.getRhoGrid()
            actual = [rho.getMinX(), rho.getMaxX(), rho.getMinY(), rho.getMaxY(),
                      rho.getMinZ()/gamma, rho.getMaxZ()/gamma]
            assert np.max(np.abs(np.array(actual)-np.array([-half_x, half_x, -half_y, half_y, -period/2, period/2]))) < 2e-12
            mesh_checks.append(dict(step=step+1, rolled_bounds_m=actual, LAB_origin_m=origin,
                                   origin_cells=origin_cells, temporary_shift_m=temporary_shift,
                                   physical_z_min_m=float(physical_z.min()), physical_z_max_m=float(physical_z.max()),
                                   maximum_physical_abs_u_m=float(np.max(np.abs(physical[:, 0]))),
                                   maximum_physical_abs_y_m=float(np.max(np.abs(physical[:, 1]))),
                                   guard_gradient_zero_from_compact_support=True))
            reset_guards(0.)
            shift.shift_Z(bunch, -temporary_shift)
            label_error = float(np.max(np.abs(np.array([bunch.z(int(i)) for i in label_indices])-labels_before)))
            assert label_error < 2e-14
            mesh_checks[-1]['maximum_unwrapped_label_roundtrip_error_m'] = label_error
        body(bunch, ds/2)
        reset_guards((step+1)*ds)
        record(step+1)
    final = np.array([[bunch.x(i), bunch.xp(i), bunch.y(i), bunch.yp(i), bunch.z(i), bunch.dE(i)]
                      for i in range(count)])
    local_positions = final[:, [0, 2, 4]].copy()
    local_positions[:, 2] += total_length
    theta = local_positions[:, 2]/radius
    positions = np.c_[local_positions[:, 0]*np.cos(theta)-2*radius*np.sin(theta/2)**2,
                      local_positions[:, 1], (radius+local_positions[:, 0])*np.sin(theta)]
    p0 = beta*gamma
    px, py = final[:, 1]*p0, final[:, 3]*p0
    normalized_kinetic = (kinetic+final[:, 5])/MASS
    pz2 = normalized_kinetic*(normalized_kinetic+2)-px*px-py*py
    assert pz2.min() > 0.
    pt = np.sqrt(pz2)
    momentum = np.c_[px*np.cos(theta)-pt*np.sin(theta), py,
                     px*np.sin(theta)+pt*np.cos(theta)]
    np.savez_compressed(args.out/'final-particles.npz', source_particle_id=ids,
                        positions_m=positions, momenta_betagamma=momentum, native_state=final,
                        positions_frenet_m=local_positions)
    if args.zero_sc:
        assert np.max(np.abs(final[:, 5])) < 1e-15
        assert np.max(np.abs(np.linalg.norm(momentum, axis=1)-beta*gamma)) < 1e-15
    assert len(rows) == steps+1 and bunch.getSize() == count+6
    assert abs(bunch.getSyncParticle().time()-steps*dt) < 1e-18
    libraries = {name: dict(path=sys.modules[name].__file__, sha256=sha(sys.modules[name].__file__))
                 for name in ['orbit.core.bunch', 'orbit.core.teapot_base', 'orbit.core.spacecharge']}
    manifest = dict(passed=True, geometry='bend', steps=steps, reference_advance_m=total_length,
                    dt_s=dt, ds_m=ds, source=source, source_json=str(args.source_json),
                    source_json_sha256=sha(args.source_json), native_driver_sha256=sha(__file__),
                    space_charge_enabled=not args.zero_sc, ZeroSC_same_nominal_charge=args.zero_sc,
                    nominal_charge_C=source['charge_C'], mesh=[nx, ny, nz], box_xy_m=[width_x, width_y],
                    longitudinal_period_m=period, physical_g=physical_g, native_encoded_step_g=native_g,
                    native_encoded_b_a=ratio, beta=beta, gamma=gamma, mass_GeV=MASS,
                    integrator='native TEAPOT body half / native transverse+LS step kick / native TEAPOT body half',
                    body_helper_sha256=sha(Path(__file__).with_name('check_bend_kinematics.py')),
                    physical_bend_radius_m=radius, physical_reference_phase_rad=total_length/radius,
                    reference_path=str(args.reference_path), reference_path_sha256=sha(args.reference_path),
                    reference_polyline_period_m=period, maximum_segment_sagitta_m=sagitta,
                    source_coordinates_for_SC='native Frenet coordinates, not temporary OPALX-polyline projections; common period/mesh pinned',
                    periodic_strategy='integer-cell rolling LAB grid origin; source labels restored after each kick, no fractional recentering',
                    periodic_check_file=str(args.periodic_check) if args.periodic_check is not None else None,
                    periodic_check_sha256=sha(args.periodic_check) if args.periodic_check is not None else None,
                    native_periodic_limit='Ideal mesh phase preserved; existing NGP0.999999 near-boundary bias quantified in the periodic check, not corrected.',
                    transport_approximation='native paraxial TEAPOT body maps; global coordinates use circular embedding of native local state',
                    native_steps=steps, mesh_unchanged_from_ten_step=True,
                    timing_boundary='inside main: source/Bunch preparation, all diagnostics and final NPZ; excludes interpreter imports',
                    thread_environment={key: __import__('os').environ.get(key) for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},
                    frenet_mean_s_convention='relative to design reference advance; reference_phase_rad=advance/R',
                    final_positions_frenet_s_convention='absolute physical circular arclength; native_state[:,4] remains relative',
                    fixed_lab_mesh=True, lab_mesh_checks=mesh_checks, guard_check=guard_check,
                    transverse_work_patch=False,
                    limitation='Native transverse SC changes xp/yp without adding transverse work to dE; preserved transparently.',
                    native_libraries=libraries, python=sys.executable, argv=sys.argv,
                    elapsed_seconds=time.monotonic()-started)
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(passed=True, seconds=manifest['elapsed_seconds'], out=str(args.out))), flush=True)


if __name__ == '__main__':
    main()
