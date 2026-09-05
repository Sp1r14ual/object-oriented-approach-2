using System;
using System.Collections.Generic;
using System.Linq;
using LinearAlgebra;

namespace Functions
{
    /// <summary>
    /// Нелинейный кубический сплайн в одномерном пространстве.
    /// Задаётся узлами сетки (x_0 < x_1 < ... < x_{k-1}).
    /// Параметрами являются значения функции в узлах: (y_0, y_1, ..., y_{k-1}).
    /// По требованию задания: сплайн НЕ реализует IDifferentiableFunction (только IFunction).
    /// </summary>
    public class CubicSplineFunction : IParametricFunction
    {
        private readonly List<double> _knots;

        public IReadOnlyList<double> Knots => _knots;
        public int KnotCount => _knots.Count;

        public CubicSplineFunction(IEnumerable<double> knots)
        {
            if (knots == null) throw new ArgumentNullException(nameof(knots));
            _knots = knots.OrderBy(x => x).ToList();
            if (_knots.Count < 3)
                throw new ArgumentException("At least 3 knots are required for a cubic spline.", nameof(knots));

            for (int i = 0; i < _knots.Count - 1; i++)
            {
                if (Math.Abs(_knots[i + 1] - _knots[i]) < 1e-14)
                    throw new ArgumentException("Knots must be strictly distinct.");
            }
        }

        public CubicSplineFunction(int knotCount, double xMin = 0.0, double xMax = 1.0)
        {
            if (knotCount < 3)
                throw new ArgumentException("At least 3 knots are required for a cubic spline.", nameof(knotCount));
            if (xMax <= xMin)
                throw new ArgumentException("xMax must be greater than xMin.");

            _knots = new List<double>(knotCount);
            double step = (xMax - xMin) / (knotCount - 1);
            for (int i = 0; i < knotCount; i++)
            {
                _knots.Add(xMin + i * step);
            }
        }

        public IFunction Bind(IVector parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != _knots.Count)
            {
                throw new ArgumentException($"Expected {_knots.Count} parameters, got {parameters.Count}.", nameof(parameters));
            }

            return new BoundCubicSplineFunction(_knots, parameters);
        }

        /// <summary>
        /// Реализует только IFunction (не реализует IDifferentiableFunction).
        /// </summary>
        private class BoundCubicSplineFunction : IFunction
        {
            private readonly List<double> _knots;
            private readonly double[] _y;
            private readonly double[] _m; // Вторые производные в узлах

            public BoundCubicSplineFunction(List<double> knots, IVector parameters)
            {
                _knots = knots;
                int n = knots.Count;
                _y = new double[n];
                for (int i = 0; i < n; i++) _y[i] = parameters[i];

                _m = ComputeSecondDerivatives(knots, _y);
            }

            private static double[] ComputeSecondDerivatives(List<double> x, double[] y)
            {
                int n = x.Count;
                double[] h = new double[n - 1];
                for (int i = 0; i < n - 1; i++)
                {
                    h[i] = x[i + 1] - x[i];
                }

                // Естественный кубический сплайн: M_0 = M_{n-1} = 0
                // Для внутренних узлов решаем трёхдиагональную систему:
                // h_{i-1} * M_{i-1} + 2*(h_{i-1} + h_i) * M_i + h_i * M_{i+1} = 6 * ((y_{i+1}-y_i)/h_i - (y_i-y_{i-1})/h_{i-1})
                int mCount = n - 2;
                if (mCount <= 0)
                {
                    return new double[n];
                }

                double[] diag = new double[mCount];
                double[] sub = new double[mCount - 1];
                double[] sup = new double[mCount - 1];
                double[] rhs = new double[mCount];

                for (int i = 1; i < n - 1; i++)
                {
                    int row = i - 1;
                    diag[row] = 2.0 * (h[i - 1] + h[i]);
                    if (row > 0) sub[row - 1] = h[i - 1];
                    if (row < mCount - 1) sup[row] = h[i];

                    double d1 = (y[i + 1] - y[i]) / h[i];
                    double d0 = (y[i] - y[i - 1]) / h[i - 1];
                    rhs[row] = 6.0 * (d1 - d0);
                }

                // Метод прогонки (алгоритм Томаса)
                double[] cPrime = new double[mCount];
                double[] dPrime = new double[mCount];

                cPrime[0] = (mCount > 1) ? sup[0] / diag[0] : 0.0;
                dPrime[0] = rhs[0] / diag[0];

                for (int i = 1; i < mCount; i++)
                {
                    double denom = diag[i] - sub[i - 1] * cPrime[i - 1];
                    if (Math.Abs(denom) < 1e-15) denom = 1e-12;

                    if (i < mCount - 1)
                    {
                        cPrime[i] = sup[i] / denom;
                    }
                    dPrime[i] = (rhs[i] - sub[i - 1] * dPrime[i - 1]) / denom;
                }

                double[] solution = new double[mCount];
                solution[mCount - 1] = dPrime[mCount - 1];
                for (int i = mCount - 2; i >= 0; i--)
                {
                    solution[i] = dPrime[i] - cPrime[i] * solution[i + 1];
                }

                double[] M = new double[n];
                M[0] = 0.0;
                M[n - 1] = 0.0;
                for (int i = 1; i < n - 1; i++)
                {
                    M[i] = solution[i - 1];
                }
                return M;
            }

            public double Value(IVector point)
            {
                if (point == null) throw new ArgumentNullException(nameof(point));
                if (point.Count == 0)
                    throw new ArgumentException("Point must contain at least 1 coordinate.");

                double xVal = point[0];
                int n = _knots.Count;

                int i;
                if (xVal <= _knots[0])
                {
                    i = 0;
                }
                else if (xVal >= _knots[n - 1])
                {
                    i = n - 2;
                }
                else
                {
                    int left = 0, right = n - 1;
                    while (right - left > 1)
                    {
                        int mid = (left + right) / 2;
                        if (_knots[mid] <= xVal) left = mid;
                        else right = mid;
                    }
                    i = left;
                }

                double h = _knots[i + 1] - _knots[i];
                double dxR = _knots[i + 1] - xVal;
                double dxL = xVal - _knots[i];

                double s = (_m[i] * dxR * dxR * dxR + _m[i + 1] * dxL * dxL * dxL) / (6.0 * h)
                         + (_y[i] * dxR + _y[i + 1] * dxL) / h
                         - (h / 6.0) * (_m[i] * dxR + _m[i + 1] * dxL);

                return s;
            }
        }
    }
}
