using System;
using Functions;
using Functionals;
using LinearAlgebra;

namespace Optimizers
{
    /// <summary>
    /// Алгоритм Гаусса-Ньютона для решения задач наименьших квадратов.
    /// По требованию задания: требует реализации ILeastSquaresFunctional.
    /// Решает систему нормальных уравнений: (J^T * J + lambda * I) * delta = -J^T * r.
    /// </summary>
    public class GaussNewtonOptimizer : IOptimizator
    {
        public int MaxIterations { get; set; } = 100;
        public double Tolerance { get; set; } = 1e-8;
        public double Damping { get; set; } = 1e-6;

        public GaussNewtonOptimizer(int maxIterations = 100, double tolerance = 1e-8, double damping = 1e-6)
        {
            MaxIterations = maxIterations;
            Tolerance = tolerance;
            Damping = damping;
        }

        public IVector Minimize(IFunctional objective,
                                 IParametricFunction function,
                                 IVector initialParameters,
                                 IVector minimumParameters = default,
                                 IVector maximumParameters = default)
        {
            if (objective == null) 
                throw new ArgumentNullException(nameof(objective));
            if (function == null) 
                throw new ArgumentNullException(nameof(function));
            if (initialParameters == null) 
                throw new ArgumentNullException(nameof(initialParameters));

            // Проверка строгого требования: требуется ILeastSquaresFunctional
            if (objective is not ILeastSquaresFunctional lsqObjective)
                throw new ArgumentException(
                    $"Objective of type '{objective.GetType().Name}' does not implement ILeastSquaresFunctional. " +
                    "GaussNewtonOptimizer requires an objective that implements ILeastSquaresFunctional.",
                    nameof(objective));
            

            var current = new Vector(initialParameters);
            Clamp(current, minimumParameters, maximumParameters);

            for (int iter = 0; iter < MaxIterations; iter++)
            {
                IFunction boundFunc = function.Bind(current);

                IVector r = lsqObjective.Residual(boundFunc);
                IMatrix J = lsqObjective.Jacobian(boundFunc);

                // Решение нормальных уравнений: (J^T * J + lambda * I) * delta = -J^T * r
                Vector delta = Matrix.SolveNormalEquations(J, r, Damping);

                double deltaNorm = delta.Norm();
                if (deltaNorm < Tolerance)
                    break;

                // Линейный поиск шага с делением пополам (backtracking)
                double currentCost = objective.Value(boundFunc);
                double stepSize = 1.0;
                Vector candidate = null;
                bool improved = false;

                for (int ls = 0; ls < 10; ls++)
                {
                    candidate = current + delta * stepSize;
                    Clamp(candidate, minimumParameters, maximumParameters);

                    double candidateCost = objective.Value(function.Bind(candidate));
                    if (candidateCost < currentCost)
                    {
                        improved = true;
                        break;
                    }
                    stepSize *= 0.5;
                }

                if (!improved)
                    break;

                current = candidate;
            }

            return current;
        }

        private static void Clamp(Vector v, IVector min, IVector max)
        {
            for (int i = 0; i < v.Count; i++)
            {
                if (min != null && min.Count > i && v[i] < min[i]) v[i] = min[i];
                if (max != null && max.Count > i && v[i] > max[i]) v[i] = max[i];
            }
        }
    }
}
