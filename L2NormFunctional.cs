using System;
using System.Collections.Generic;
using Functions;
using LinearAlgebra;

namespace Functionals
{
    /// <summary>
    /// l2 норма разности с требуемыми значениями в наборе точек.
    /// F(f) = sqrt(sum (f(x_i) - y_i)^2)  (или sum (f(x_i) - y_i)^2 при Squared = true)
    /// По требованию задания: реализует IDifferentiableFunctional И ILeastSquaresFunctional.
    /// </summary>
    public class L2NormFunctional : IDifferentiableFunctional, ILeastSquaresFunctional
    {
        private readonly List<(IVector Point, double Value)> _dataPoints;

        public IReadOnlyList<(IVector Point, double Value)> DataPoints => _dataPoints;

        /// <summary>
        /// Если true, Value возвращает сумму квадратов невязок sum r_i^2,
        /// если false (по умолчанию), возвращает евклидову l2-норму sqrt(sum r_i^2).
        /// </summary>
        public bool Squared { get; set; } = false;

        public L2NormFunctional(IEnumerable<(IVector Point, double Value)> dataPoints, bool squared = false)
        {
            if (dataPoints == null) 
                throw new ArgumentNullException(nameof(dataPoints));
            _dataPoints = new List<(IVector Point, double Value)>(dataPoints);
            if (_dataPoints.Count == 0)
                throw new ArgumentException("At least one data point is required.", nameof(dataPoints));
            Squared = squared;
        }

        public L2NormFunctional(IEnumerable<(double X, double Y)> points1D, bool squared = false)
        {
            if (points1D == null) 
                throw new ArgumentNullException(nameof(points1D));
            _dataPoints = new List<(IVector Point, double Value)>();
            foreach (var (x, y) in points1D)
            {
                _dataPoints.Add((new Vector(x), y));
            }
            if (_dataPoints.Count == 0)
                throw new ArgumentException("At least one data point is required.", nameof(points1D));
            Squared = squared;
        }

        public double Value(IFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));

            double sumSq = 0.0;
            foreach (var (pt, y) in _dataPoints)
            {
                double diff = function.Value(pt) - y;
                sumSq += diff * diff;
            }

            return Squared ? sumSq : Math.Sqrt(sumSq);
        }

        /// <summary>
        /// Вычисляет вектор невязок r_i = f(x_i) - y_i.
        /// </summary>
        public IVector Residual(IFunction function)
        {
            if (function == null) 
                throw new ArgumentNullException(nameof(function));

            var res = new Vector(_dataPoints.Count);
            foreach (var (pt, y) in _dataPoints)
            {
                res.Add(function.Value(pt) - y);
            }
            return res;
        }

        /// <summary>
        /// Вычисляет матрицу Якоби невязок по параметрам:
        /// J_{i, j} = d r_i / d theta_j = d f(x_i) / d theta_j.
        /// Требует, чтобы функция реализовывала IDifferentiableFunction.
        /// </summary>
        public IMatrix Jacobian(IFunction function)
        {
            if (function == null) 
                throw new ArgumentNullException(nameof(function));
            if (function is not IDifferentiableFunction diffFunc)
            {
                throw new ArgumentException(
                    $"Function of type '{function.GetType().Name}' does not implement IDifferentiableFunction. " +
                    "Jacobian matrix cannot be computed.", nameof(function));
            }

            var jacobian = new Matrix();
            foreach (var (pt, _) in _dataPoints)
            {
                IVector grad = diffFunc.Gradient(pt);
                jacobian.Add(new Vector(grad));
            }
            return jacobian;
        }

        /// <summary>
        /// Вычисляет градиент функционала по параметрам:
        /// Если Squared = false: grad F = (J^T * r) / ||r||
        /// Если Squared = true:  grad F = 2 * (J^T * r)
        /// </summary>
        public IVector Gradient(IFunction function)
        {
            if (function == null) 
                throw new ArgumentNullException(nameof(function));
            if (function is not IDifferentiableFunction diffFunc)
                throw new ArgumentException(
                    $"Function of type '{function.GetType().Name}' does not implement IDifferentiableFunction. " +
                    "Analytical gradient cannot be computed.", nameof(function));
            

            var r = Residual(function);
            var J = Jacobian(diffFunc);

            int m = _dataPoints.Count;
            int p = J.Count > 0 ? J[0].Count : 0;

            var grad = Vector.CreateZeros(p);
            for (int i = 0; i < m; i++)
            {
                double ri = r[i];
                for (int j = 0; j < p; j++)
                {
                    grad[j] += ri * J[i][j];
                }
            }

            if (Squared)
                return grad * 2.0;
            else
            {
                double norm = ((Vector)r).Norm();
                if (norm < 1e-15)
                    return Vector.CreateZeros(p);
                return grad / norm;
            }
        }
    }
}
