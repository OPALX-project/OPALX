#!/usr/bin/env python3
"""Native TEAPOT sign/normalization and body-map controls, without space charge."""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from orbit.core.bunch import Bunch
from orbit.core import teapot_base as TPB

HERE = Path(__file__).resolve().parent
BETA, MASS, RADIUS, LENGTH, STEPS = .001, .938272029, 2., 1., 10


def body(bunch, length):
    angle = length/RADIUS
    TPB.bend1(bunch, length/2, angle/2)
    TPB.bend2(bunch, length/2)
    TPB.bend3(bunch, angle/2)
    TPB.bend4(bunch, angle/2)
    TPB.bend4(bunch, angle/2)
    TPB.bend3(bunch, angle/2)
    TPB.bend2(bunch, length/2)
    TPB.bend1(bunch, length/2, angle/2)


def main():
    initial = np.array([[-.01, 0., -1.], [-.0001, 0., .5], [0., .003, 0.],
                        [.0001, 0., -.5], [.01, 0., 1.]])
    gamma = 1/math.sqrt(1-BETA*BETA)
    outputs = {}
    for name in ['linear_bend1', 'native_body']:
        bunch = Bunch(); bunch.mass(MASS); bunch.charge(1.)
        bunch.getSyncParticle().kinEnergy(MASS*BETA*BETA*gamma*gamma/(gamma+1))
        for u, y, s in initial:
            bunch.addParticle(float(u), 0., float(y), 0., float(s), 0.)
        for _ in range(STEPS):
            for _ in range(2):
                if name == 'linear_bend1':
                    TPB.bend1(bunch, LENGTH/(2*STEPS), LENGTH/(2*STEPS*RADIUS))
                else:
                    body(bunch, LENGTH/(2*STEPS))
        actual = np.array([[bunch.x(i), bunch.xp(i), bunch.y(i), bunch.yp(i), bunch.z(i), bunch.dE(i)]
                           for i in range(len(initial))])
        phase = LENGTH/RADIUS
        expected = np.c_[initial[:, 0]*math.cos(phase), -initial[:, 0]*math.sin(phase)/RADIUS,
                         initial[:, 1], np.zeros(len(initial)),
                         initial[:, 2]-initial[:, 0]*math.sin(phase), np.zeros(len(initial))]
        delta = actual-expected
        if name == 'linear_bend1':
            assert np.max(np.abs(delta)) < 2e-15
        assert actual[-1, 1] < 0 and actual[-1, 4] < initial[-1, 2]
        outputs[name] = dict(state=actual.tolist(), maximum_difference_from_linear_formula=float(np.max(np.abs(delta))),
                             sync_time_s=bunch.getSyncParticle().time())
    lib = Path(TPB.__file__)
    result = dict(passed=True, source='native teapotbase.cc bend1/bend2/bend3/bend4, no field or solver patch',
                  initial_material_uys_m=initial.tolist(), radius_m=RADIUS, advance_m=LENGTH, steps=STEPS,
                  expected_linear_u='u0*cos(phi)', expected_linear_xp='-u0*sin(phi)/R',
                  expected_linear_phase_shift_m='-u0*sin(phi)',
                  xp_definition='p_radial/p_sync; confirmed by native linac drift and SynchPartRedefinition energy transform',
                  outputs=outputs, native_library=str(lib), native_library_sha256=hashlib.sha256(lib.read_bytes()).hexdigest())
    (HERE/'bend-kinematics.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
