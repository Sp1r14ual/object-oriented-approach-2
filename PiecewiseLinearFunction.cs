using System;
using System.Collections.Generic;
using System.Linq;
using LinearAlgebra;

namespace Functions
{
    /// <summary>
    /// Кусочно-линейная функция в одномерном пространстве.
    /// Задаётся узлами сетки (x_0 < x_1 < ... < x_{k-1}).
    /// Параметрами являются значения функции в узлах: (y_0, y_1, ..., y_{k-1}).
    /// Реализует IDifferentiableFunction.
    /// </summary>
    public class PiecewiseLinearFunction : IParametricFunction
    {
        private readonly List<double> _knots;

        public IReadOnlyList<double> Knots => _knots;
        public int KnotCount => _knots.Count;

        public PiecewiseLinearFunction(IEnumerable<double> knots)
        {
            if (knots == null) 
                throw new ArgumentNullException(nameof(knots));
            _knots = knots.OrderBy(x => x).ToList();
            if (_knots.Count < 2)
                throw new ArgumentException("At least 2 knots are required for a piecewise linear function.", nameof(knots));

            // Проверка на совпадение узлов
            for (int i = 0; i < _knots.Count - 1; i++)
            {
                if (Math.Abs(_knots[i + 1] - _knots[i]) < 1e-14)
                    throw new ArgumentException("Knots must be strictly distinct.");
            }
        }

        public PiecewiseLinearFunction(int knotCount, double xMin = 0.0, double xMax = 1.0)
        {
            if (knotCount < 2)
                throw new ArgumentException("At least 2 knots are required.", nameof(knotCount));
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
            if (parameters == null) 
                throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != _knots.Count)
                throw new ArgumentException($"Expected {_knots.Count} parameters (values at knots), but got {parameters.Count}.", nameof(parameters));

            return new BoundPiecewiseLinearFunction(_knots, parameters);
        }

        private class BoundPiecewiseLinearFunction : IDifferentiableFunction
        {
            private readonly List<double> _knots;
            private readonly Vector _values;

            public BoundPiecewiseLinearFunction(List<double> knots, IVector values)
            {
                _knots = knots;
                _values = new Vector(values);
            }

            private (int index, double t) GetInterval(double x)
            {
                int k = _knots.Count;
                if (x <= _knots[0])
                {
                    double t = (x - _knots[0]) / (_knots[1] - _knots[0]);
                    return (0, t);
                }
                if (x >= _knots[k - 1])
                {
                    double t = (x - _knots[k - 2]) / (_knots[k - 1] - _knots[k - 2]);
                    return (k - 2, t);
                }

                // Бинарный поиск интервала
                int left = 0;
                int right = k - 1;
                while (right - left > 1)
                {
                    int mid = (left + right) / 2;
                    if (_knots[mid] <= x) left = mid;
                    else right = mid;
                }

                double frac = (x - _knots[left]) / (_knots[left + 1] - _knots[left]);
                return (left, frac);
            }

            public double Value(IVector point)
            {
                if (point == null) throw new ArgumentNullException(nameof(point));
                if (point.Count == 0)
                    throw new ArgumentException("Point must contain at least 1 coordinate.");

                double x = point[0];
                var (idx, t) = GetInterval(x);
                return (1.0 - t) * _values[idx] + t * _values[idx + 1];
            }

            /// <summary>
            /// Градиент функции по параметрам (значениям в узлах y_0, ..., y_{k-1}).
            /// df/dy_idx = 1 - t, df/dy_{idx+1} = t, остальные 0.
            /// </summary>
            public IVector Gradient(IVector point)
            {
                if (point == null) throw new ArgumentNullException(nameof(point));
                if (point.Count == 0)
                    throw new ArgumentException("Point must contain at least 1 coordinate.");

                double x = point[0];
                var (idx, t) = GetInterval(x);

                var grad = Vector.CreateZeros(_knots.Count);
                grad[idx] = 1.0 - t;
                grad[idx + 1] = t;
                return grad;
            }
        }
    }
}
