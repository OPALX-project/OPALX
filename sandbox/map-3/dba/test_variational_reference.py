"""Independent checks of coordinates, analytic sensitivities and face events."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from analytic_fields import load_parameters
from validate_opalx_fields import global_field
from variational_reference import (
    Face, PlanarField, VariationalDynamics, integrate, make_reference, particle_constants,
)

spec = importlib.util.spec_from_file_location("map2_exact", Path(__file__).resolve().parents[2] / "map-2/check_maps.py")
exact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exact)


class VariationalTests(unittest.TestCase):
    def test_drift_full_matrix(self):
        field = lambda r, mode: (np.zeros(3), np.zeros((3, 3)))
        matrix, _ = integrate(VariationalDynamics(field, exact.P, 1), 1.7)
        np.testing.assert_allclose(matrix, exact.drift(1.7), rtol=0, atol=2e-14)

    def test_quadrupole_full_matrix(self):
        gradient = 1.3
        jac = np.array([[0., gradient, 0], [gradient, 0, 0], [0, 0, 0]])
        field = lambda r, mode: (jac @ r, jac)
        matrix, _ = integrate(VariationalDynamics(field, exact.P, -1), 0.4)
        np.testing.assert_allclose(matrix, exact.quadrupole(gradient, 0.4), rtol=0, atol=3e-14)

    def test_uniform_bend_full_relativistic_matrix(self):
        rho, angle = 2.0, 0.4
        field = lambda r, mode: (np.array([0., 1 / rho, 0.]), np.zeros((3, 3)))
        matrix, _ = integrate(VariationalDynamics(field, exact.P, 1), rho * angle)
        np.testing.assert_allclose(matrix, exact.sector_bend(rho, angle), rtol=0, atol=3e-14)

    def test_field_jacobian_including_vertical_fringe(self):
        p = load_parameters()
        field = PlanarField(p)
        for row, origin, rotation in field.poses:
            for z in (-0.15, 0.23, p.body_length + 0.19) if row.kind == "bend" else (0.1,):
                point = origin + rotation @ np.array([0.008, 0, z])
                mode = row.kind == "quadrupole"
                b, jac = field(point, mode)
                np.testing.assert_allclose(b, global_field(point, p)[0], rtol=0, atol=3e-14)
                numeric = np.column_stack([
                    (global_field(point + step, p)[0] - global_field(point - step, p)[0]) / 2e-6
                    for step in np.eye(3) * 1e-6])
                np.testing.assert_allclose(jac, numeric, rtol=2e-7, atol=1e-9)

    def test_lorentz_jacobian(self):
        # Affine field lets all seven state variables be independently perturbed.
        gradient = np.array([[.2, .3, -.1], [.3, -.2, .4], [-.1, .4, 0]])
        field = lambda r, mode: (np.array([.2, -.4, .3]) + gradient @ r, gradient)
        dynamics = VariationalDynamics(field, 3.0, -.7)
        state = np.array([.01, .03, .02, .07, -.02, 1.03, .2])
        _, jac = dynamics.values(state)
        numeric = np.column_stack([(dynamics.values(state + d)[0] - dynamics.values(state - d)[0]) / 2e-6
                                   for d in np.eye(7) * 1e-6])
        np.testing.assert_allclose(jac, numeric, rtol=2e-7, atol=2e-10)

    def test_off_axis_quad_face_sensitivity(self):
        # The displaced quad has nonzero B on the reference at both hard faces.
        # Compare against independently integrated, event-located nonlinear rays.
        gradient = np.array([[0., 1.3, 0], [1.3, 0, 0], [0, 0, 0]])
        def field(r, mode):
            return (gradient @ (r - [.02, 0, 0]), gradient) if mode else (np.zeros(3), np.zeros((3, 3)))
        dynamics = VariationalDynamics(field, 3.0, -1)
        faces = [Face("in", np.array([0., 0, .3]), np.array([0., 0, 1]), True),
                 Face("out", np.array([0., 0, .8]), np.array([0., 0, 1]), False)]
        matrix, info = integrate(dynamics, 1.3, faces, rtol=3e-14, atol=3e-16, max_step=.01)
        ref = np.array(info["final_state"])
        norm = np.linalg.norm(ref[3:6])
        n = ref[3:6] / norm
        ex = np.array([n[2], 0, -n[0]])
        beta_ratio = dynamics.u0 * norm / np.sqrt(1 + (dynamics.u0 * norm)**2) / dynamics.beta0

        def ray(offset):
            x, xp, y, yp, zeta, delta = offset
            k = (1 + delta) * np.array([xp, yp, 1.]) / np.sqrt(1 + xp*xp + yp*yp)
            state = np.r_[x, y, 0., k, -zeta]
            s, mode = 0., False
            end = Face("measurement", ref[:3], n, False)
            for face in [*faces, end]:
                def event(t, z):
                    return (z[:3] - face.origin) @ face.normal
                event.terminal, event.direction = True, 1
                result = solve_ivp(lambda t, z: dynamics.values(z, mode)[0], (s, 1.5), state,
                                   method="DOP853", events=event, rtol=3e-14, atol=3e-16, max_step=.01)
                self.assertEqual(result.status, 1)
                s, state = result.t[-1], result.y[:, -1]
                mode = face.next_quad_on
            return np.array([ex @ (state[:3] - ref[:3]), ex @ state[3:6] / (n @ state[3:6]),
                             state[1], state[4] / (n @ state[3:6]), -beta_ratio * (state[6] - ref[6]),
                             np.linalg.norm(state[3:6]) / norm - 1])
        numeric = np.column_stack([(ray(d) - ray(-d)) / 2e-5 for d in np.eye(6) * 1e-5])
        np.testing.assert_allclose(matrix, numeric, rtol=0, atol=2e-9)

    def test_dba_cross_method_refinement(self):
        p = load_parameters()
        fine, information = make_reference(p, method="DOP853", rtol=3e-14, atol=3e-16, max_step=.005)
        check, _ = make_reference(p, method="RK45", rtol=3e-14, atol=3e-16, max_step=.005)
        np.testing.assert_allclose(fine, check, rtol=0, atol=1e-9)
        self.assertLess(information["max_relative_momentum_drift"], 1e-12)

    def test_particle_momentum_not_silently_replaced_by_field_P0(self):
        u0, mass = particle_constants()
        self.assertAlmostEqual(u0 * mass, .2505104777748827, places=15)
        self.assertGreater(abs(u0 * mass - load_parameters().momentum_GeV_c), 3e-10)


if __name__ == "__main__":
    unittest.main()
