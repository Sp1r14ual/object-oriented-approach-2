using System;
using Functions;
using Functionals;
using LinearAlgebra;

namespace Optimizers
{
    /// <summary>
    /// Универсальный метод оптимизации: классический метод Монте-Карло (случайный поиск).
    /// </summary>
    public class MonteCarloOptimizer : IOptimizator
    {
        public int MaxIterations { get; set; } = 50000;
        public double StepScale { get; set; } = 2.0;
        public int? Seed { get; set; }

        public MonteCarloOptimizer(int maxIterations = 50000, double stepScale = 2.0, int? seed = 0)
        {
            MaxIterations = maxIterations;
            StepScale = stepScale;
            Seed = seed;
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

            var rand = Seed.HasValue ? new Random(Seed.Value) : new Random();
            int pCount = initialParameters.Count;

            var best = new Vector(initialParameters);
            Clamp(best, minimumParameters, maximumParameters);
            double bestVal = objective.Value(function.Bind(best));

            for (int k = 0; k < MaxIterations; k++)
            {
                var candidate = new Vector(pCount);
                for (int i = 0; i < pCount; i++)
                {
                    double val;
                    if (minimumParameters != null && maximumParameters != null &&
                        minimumParameters.Count > i && maximumParameters.Count > i)
                    {
                        val = minimumParameters[i] + rand.NextDouble() * (maximumParameters[i] - minimumParameters[i]);
                    }
                    else
                    {
                        val = initialParameters[i] + (rand.NextDouble() * 2.0 - 1.0) * StepScale;
                    }
                    candidate.Add(val);
                }

                Clamp(candidate, minimumParameters, maximumParameters);

                double currentVal = objective.Value(function.Bind(candidate));
                if (currentVal < bestVal)
                {
                    bestVal = currentVal;
                    best = candidate;
                }
            }

            return best;
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
