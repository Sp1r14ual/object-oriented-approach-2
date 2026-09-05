using System;
using Functions;
using Functionals;
using LinearAlgebra;

namespace Optimizers
{
    /// <summary>
    /// Метод градиентного спуска с дроблением шага (Backtracking Armijo line search).
    /// По требованию задания: требует реализации IDifferentiableFunctional.
    /// </summary>
    public class GradientDescentOptimizer : IOptimizator
    {
        public int MaxIterations { get; set; } = 2000;
        public double Tolerance { get; set; } = 1e-6;
        public double InitialLearningRate { get; set; } = 1.0;

        public GradientDescentOptimizer(int maxIterations = 2000, double tolerance = 1e-6, double initialLearningRate = 1.0)
        {
            MaxIterations = maxIterations;
            Tolerance = tolerance;
            InitialLearningRate = initialLearningRate;
        }

        public IVector Minimize(IFunctional objective,
                                 IParametricFunction function,
                                 IVector initialParameters,
                                 IVector minimumParameters = default,
                                 IVector maximumParameters = default)
        {
            if (objective == null) throw new ArgumentNullException(nameof(objective));
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (initialParameters == null) throw new ArgumentNullException(nameof(initialParameters));

            if (objective is not IDifferentiableFunctional diffObjective)
            {
                throw new ArgumentException(
                    $"Objective of type '{objective.GetType().Name}' does not implement IDifferentiableFunctional. " +
                    "GradientDescentOptimizer requires an objective that implements IDifferentiableFunctional.",
                    nameof(objective));
            }

            var current = new Vector(initialParameters);
            Clamp(current, minimumParameters, maximumParameters);

            for (int k = 0; k < MaxIterations; k++)
            {
                IFunction boundFunc = function.Bind(current);
                double currentVal = diffObjective.Value(boundFunc);
                var grad = new Vector(diffObjective.Gradient(boundFunc));

                double gradNorm = grad.Norm();
                if (gradNorm < Tolerance)
                {
                    break;
                }

                // Направление антиградиента
                var dir = grad * -1.0;
                double slope = -gradNorm * gradNorm;

                // Backtracking line search
                double alpha = InitialLearningRate;
                const double c1 = 1e-4;
                const double rho = 0.5;
                Vector next = null;
                bool stepFound = false;

                for (int ls = 0; ls < 30; ls++)
                {
                    next = current + dir * alpha;
                    Clamp(next, minimumParameters, maximumParameters);

                    double nextVal = diffObjective.Value(function.Bind(next));
                    if (nextVal <= currentVal + c1 * alpha * slope)
                    {
                        stepFound = true;
                        break;
                    }
                    alpha *= rho;
                }

                if (!stepFound || (next - current).Norm() < 1e-12)
                {
                    break;
                }

                current = next;
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
