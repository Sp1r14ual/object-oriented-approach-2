using System;
using System.Collections.Generic;
using Functions;
using LinearAlgebra;

namespace Functionals
{
    /// <summary>
    /// linf норма (чебышёвская норма) разности с требуемыми значениями в наборе точек.
    /// F(f) = max |f(x_i) - y_i|
    /// По требованию задания: НЕ реализует IDifferentiableFunctional, НЕ реализует ILeastSquaresFunctional.
    /// </summary>
    public class LInfNormFunctional : IFunctional
    {
        private readonly List<(IVector Point, double Value)> _dataPoints;

        public IReadOnlyList<(IVector Point, double Value)> DataPoints => _dataPoints;

        public LInfNormFunctional(IEnumerable<(IVector Point, double Value)> dataPoints)
        {
            if (dataPoints == null) throw new ArgumentNullException(nameof(dataPoints));
            _dataPoints = new List<(IVector Point, double Value)>(dataPoints);
            if (_dataPoints.Count == 0)
                throw new ArgumentException("At least one data point is required.", nameof(dataPoints));
        }

        public LInfNormFunctional(IEnumerable<(double X, double Y)> points1D)
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

            double maxDiff = 0.0;
            foreach (var (pt, y) in _dataPoints)
            {
                double diff = Math.Abs(function.Value(pt) - y);
                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
            return maxDiff;
        }
    }
}
