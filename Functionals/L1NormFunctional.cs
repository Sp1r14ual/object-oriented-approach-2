using System;
using System.Collections.Generic;
using Functions;
using LinearAlgebra;

namespace Functionals
{
    /// <summary>
    /// l1 норма разности с требуемыми значениями в наборе точек.
    /// F(f) = sum |f(x_i) - y_i|
    /// По требованию задания: реализует IDifferentiableFunctional, НЕ реализует ILeastSquaresFunctional.
    /// </summary>
    public class L1NormFunctional : IDifferentiableFunctional
    {
        private readonly List<(IVector Point, double Value)> _dataPoints;

        public IReadOnlyList<(IVector Point, double Value)> DataPoints => _dataPoints;

        public L1NormFunctional(IEnumerable<(IVector Point, double Value)> dataPoints)
        {
            if (dataPoints == null) throw new ArgumentNullException(nameof(dataPoints));
            _dataPoints = new List<(IVector Point, double Value)>(dataPoints);
            if (_dataPoints.Count == 0)
                throw new ArgumentException("At least one data point is required.", nameof(dataPoints));
        }

        public L1NormFunctional(IEnumerable<(double X, double Y)> points1D)
        {
            if (points1D == null) throw new ArgumentNullException(nameof(points1D));
            _dataPoints = new List<(IVector Point, double Value)>();
            foreach (var (x, y) in points1D)
            {
                _dataPoints.Add((new Vector(x), y));
            }
            if (_dataPoints.Count == 0)
                throw new ArgumentException("At least one data point is required.", nameof(points1D));
        }

        public double Value(IFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));

            double sum = 0.0;
            foreach (var (pt, y) in _dataPoints)
            {
                double diff = function.Value(pt) - y;
                sum += Math.Abs(diff);
            }
            return sum;
        }

        /// <summary>
        /// Вычисляет градиент l1-нормы по параметрам функции:
        /// grad F = sum sgn(f(x_i) - y_i) * grad f(x_i).
        /// Требует, чтобы переданная функция реализовывала IDifferentiableFunction.
        /// </summary>
        public IVector Gradient(IFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (function is not IDifferentiableFunction diffFunc)
            {
                throw new ArgumentException(
                    $"Function of type '{function.GetType().Name}' does not implement IDifferentiableFunction. " +
                    "Analytical gradient cannot be computed.", nameof(function));
            }

            Vector totalGrad = null;
            foreach (var (pt, y) in _dataPoints)
            {
                double diff = diffFunc.Value(pt) - y;
                double sgn = diff > 1e-12 ? 1.0 : (diff < -1e-12 ? -1.0 : 0.0);

                IVector ptGrad = diffFunc.Gradient(pt);
                if (totalGrad == null)
                {
                    totalGrad = Vector.CreateZeros(ptGrad.Count);
                }

                for (int j = 0; j < ptGrad.Count; j++)
                {
                    totalGrad[j] += sgn * ptGrad[j];
                }
            }

            return totalGrad ?? new Vector();
        }
    }
}
