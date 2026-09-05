using System;
using Functions;
using Functionals;
using LinearAlgebra;

namespace Optimizers
{
    /// <summary>
    /// Универсальный метод оптимизации: алгоритм имитации отжига (Simulated Annealing).
    /// Не требует дифференцируемости и структуры наименьших квадратов,
    /// работает с любым IFunctional и IParametricFunction через стандартные интерфейсы.
    /// </summary>
    public class SimulatedAnnealingOptimizer : IOptimizator
    {
        public int MaxIterations { get; set; } = 10000;
        public double InitialTemperature { get; set; } = 100.0;
        public double CoolingFactor { get; set; } = 0.995;
        public double StepScale { get; set; } = 1.0;
        public int? Seed { get; set; }

        public SimulatedAnnealingOptimizer(int maxIterations = 10000, double initialTemperature = 100.0, double coolingFactor = 0.995, int? seed = null)
        {
            MaxIterations = maxIterations;
            InitialTemperature = initialTemperature;
            CoolingFactor = coolingFactor;
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

            var current = new Vector(initialParameters);
            Clamp(current, minimumParameters, maximumParameters);

            double currentEnergy = objective.Value(function.Bind(current));

            var best = current.Clone();
            double bestEnergy = currentEnergy;

            double temperature = InitialTemperature;

            for (int k = 0; k < MaxIterations; k++)
            {
                if (temperature < 1e-12) break;

                // Генерация соседней точки
                var candidate = new Vector(pCount);
                for (int i = 0; i < pCount; i++)
                {
                    double delta;
                    if (minimumParameters != null && maximumParameters != null &&
                        minimumParameters.Count > i && maximumParameters.Count > i)
                    {
                        double span = maximumParameters[i] - minimumParameters[i];
                        delta = (rand.NextDouble() * 2.0 - 1.0) * span * (temperature / InitialTemperature + 0.05);
                    }
                    else
                    {
                        delta = (rand.NextDouble() * 2.0 - 1.0) * StepScale * (temperature / InitialTemperature + 0.05);
                    }

                    candidate.Add(current[i] + delta);
                }

                Clamp(candidate, minimumParameters, maximumParameters);

                double candidateEnergy = objective.Value(function.Bind(candidate));
                double deltaEnergy = candidateEnergy - currentEnergy;

                // Критерий Метрополиса
                if (deltaEnergy < 0 || rand.NextDouble() < Math.Exp(-deltaEnergy / temperature))
                {
                    current = candidate;
                    currentEnergy = candidateEnergy;

                    if (currentEnergy < bestEnergy)
                    {
                        best = current.Clone();
                        bestEnergy = currentEnergy;
                    }
                }

                temperature *= CoolingFactor;
            }

            return best;
        }

        private static void Clamp(Vector v, IVector min, IVector max)
        {
            for (int i = 0; i < v.Count; i++)
            {
                if (min != null && min.Count > i && v[i] < min[i])
                {
                    v[i] = min[i];
                }
                if (max != null && max.Count > i && v[i] > max[i])
                {
                    v[i] = max[i];
                }
            }
        }
    }
}
