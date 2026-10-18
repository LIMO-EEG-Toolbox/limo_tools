"""Independent arithmetic checks for MATLAB fixtures; NOT execution of MATLAB."""
import unittest
import numpy as np


class NumericalReferences(unittest.TestCase):
    def test_error_metrics(self):
        difference = np.array([3.0, 4.0])
        self.assertEqual(np.abs(difference).mean(), 3.5)
        self.assertAlmostEqual(np.sqrt(np.square(difference).mean()), np.sqrt(12.5))

    def test_fitted_response(self):
        design = np.array([[1, 0], [0, 2], [1, -1]])
        beta_delta = np.array([[1, 2], [3, 4]])
        np.testing.assert_array_equal(design @ beta_delta, [[1, 2], [6, 8], [-2, -2]])
        self.assertEqual(np.abs(design @ beta_delta).mean(), 3.5)
        self.assertEqual(np.abs(design[:, :1] @ beta_delta[:1]).mean(), 1)

    def test_bisquare_curve(self):
        z = np.array([0, 3, 5.3425, 7.685, 10])
        excess = np.maximum(z - 3, 0)
        q = np.square(1 - np.square(np.minimum(excess / 4.685, 1)))
        np.testing.assert_allclose(q, [1, 1, 0.5625, 0, 0], atol=1e-14)

    def test_robust_scores_normalization_and_ess(self):
        errors = np.array([1., 2., 3., 4., 20.])
        location = np.median(errors)
        scale = 1.4826 * np.median(np.abs(errors - location))
        z = (errors - location) / scale
        q = np.square(1 - np.square(np.minimum(np.maximum(z - 3, 0) / 4.685, 1)))
        np.testing.assert_array_equal(q, [1, 1, 1, 1, 0])
        np.testing.assert_allclose(q / q.mean(), [1.25, 1.25, 1.25, 1.25, 0])
        self.assertEqual(q.sum() ** 2 / np.square(q).sum(), 4)

    def test_degenerate_scale_fixtures(self):
        for values, expected_iqr in [([0, 0, 0, 1, 2], 1), ([0, 0, 0, 0, 2], 0)]:
            values = np.array(values)
            self.assertEqual(np.median(np.abs(values - np.median(values))), 0)
            self.assertEqual(np.diff(np.quantile(values, [.25, .75], method='linear'))[0], expected_iqr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
