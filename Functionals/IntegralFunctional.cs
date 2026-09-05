using System;
using System.Collections.Generic;
using Functions;
using LinearAlgebra;

namespace Functionals
{
    /// <summary>
    /// Численный интеграл функции по некоторой n-мерной параллелепипедной области:
    /// Omega = [a_0, b_0] x [a_1, b_1] x ... x [a_{n-1}, b_{n-1}]
    /// F(f) = integral_Omega f(x) dx
    /// По требованию задания: реализует IFunctional.
    /// </summary>
    public class IntegralFunctional : IFunctional
    {
        private readonly List<(double Min, double Max)> _bounds;
        private readonly int _stepsPerDimension;

        public IReadOnlyList<(double Min, double Max)> Bounds => _bounds;
        public int StepsPerDimension => _stepsPerDimension;

        /// <summary>
        /// Интеграл по 1D-интервалу [a, b].
        /// </summary>
        public IntegralFunctional(double a, double b, int steps = 100)
        {
            if (b <= a) throw new ArgumentException("Upper bound b must be greater than lower bound a.");
            if (steps < 2) throw new ArgumentException("Steps must be at least 2.");
            if (steps % 2 != 0) steps++; // Для формулы Симпсона требуется чётное число интервалов

            _bounds = new List<(double Min, double Max)> { (a, b) };
            _stepsPerDimension = steps;
        }

        /// <summary>
        /// Интеграл по n-мерному брусу [a_0, b_0] x ... x [a_{n-1}, b_{n-1}].
        /// </summary>
        public IntegralFunctional(IEnumerable<(double Min, double Max)> bounds, int stepsPerDimension = 20)
        {
            if (bounds == null) throw new ArgumentNullException(nameof(bounds));
            _bounds = new List<(double Min, double Max)>(bounds);
            if (_bounds.Count == 0)
                throw new ArgumentException("At least one dimension is required.", nameof(bounds));

            foreach (var (min, max) in _bounds)
            {
                if (max <= min)
                    throw new ArgumentException($"Invalid interval [{min}, {max}]: max must be greater than min.");
            }

            if (stepsPerDimension < 2) stepsPerDimension = 2;
            if (stepsPerDimension % 2 != 0) stepsPerDimension++;
            _stepsPerDimension = stepsPerDimension;
        }

        public double Value(IFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));

            if (_bounds.Count == 1)
            {
                // Составная формула Симпсона (1/3) в 1D
                var (a, b) = _bounds[0];
                int n = _stepsPerDimension;
                double h = (b - a) / n;

                double sum = function.Value(new Vector(a)) + function.Value(new Vector(b));

                for (int i = 1; i < n; i++)
                {
                    double x = a + i * h;
                    double weight = (i % 2 == 1) ? 4.0 : 2.0;
                    sum += weight * function.Value(new Vector(x));
                }

                return sum * (h / 3.0);
            }
            else
            {
                // Многомерное интегрирование методом прямоугольников по регулярной сетке
                return IntegrateRecursive(function, 0, new double[_bounds.Count]);
            }
        }

        private double IntegrateRecursive(IFunction function, int dim, double[] currentPoint)
        {
            if (dim == _bounds.Count)
            {
                return function.Value(new Vector(currentPoint));
            }

            var (a, b) = _bounds[dim];
            int n = _stepsPerDimension;
            double h = (b - a) / n;
            double dimSum = 0.0;

            for (int i = 0; i < n; i++)
            {
                // Средняя точка отрезка (метод средних прямоугольников)
                currentPoint[dim] = a + (i + 0.5) * h;
                dimSum += IntegrateRecursive(function, dim + 1, currentPoint);
            }

            return dimSum * h;
        }
    }
}
