using System;
using Functions;
using Functionals;
using LinearAlgebra;

namespace Optimizers
{
    /// <summary>
    /// Метод нелинейных сопряжённых градиентов (Polak-Ribière с рестартом).
    /// По требованию задания: требует реализации IDifferentiableFunctional.
    /// </summary>
    public class ConjugateGradientOptimizer : IOptimizator
    {
        public int MaxIterations { get; set; } = 1000;
        public double Tolerance { get; set; } = 1e-6;
        public double InitialStep { get; set; } = 1.0;

        public ConjugateGradientOptimizer(int maxIterations = 1000, double tolerance = 1e-6, double initialStep = 1.0)
        {
            MaxIterations = maxIterations;
            Tolerance = tolerance;
            InitialStep = initialStep;
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

            // Проверка строгого требования: требуется IDifferentiableFunctional
            if (objective is not IDifferentiableFunctional diffObjective)
            {
                throw new ArgumentException(
                    $"Objective of type '{objective.GetType().Name}' does not implement IDifferentiableFunctional. " +
                    "ConjugateGradientOptimizer requires an objective that implements IDifferentiableFunctional.",
                    nameof(objective));
            }

            int p = initialParameters.Count;
            var current = new Vector(initialParameters);
            Clamp(current, minimumParameters, maximumParameters);

            IFunction boundFunc = function.Bind(current);
            double currentVal = diffObjective.Value(boundFunc);
            var g = new Vector(diffObjective.Gradient(boundFunc));
            var d = g * -1.0; // Начальное направление спуска

            for (int k = 0; k < MaxIterations; k++)
            {
                double gNorm = g.Norm();
                if (gNorm < Tolerance)
                {
                    break;
                }

                // Поиск шага одномерной минимизацией (backtracking Armijo line search)
                double alpha = InitialStep;
                double slope = g.Dot(d);
                if (slope > 0) // Если направление не является спусковым, сбрасываем в антиградиент
                {
                    d = g * -1.0;
                    slope = g.Dot(d);
                }

                const double c1 = 1e-4;
                const double rho = 0.5;
                Vector candidate = null;
                double candidateVal = currentVal;
                bool stepAccepted = false;

                for (int ls = 0; ls < 30; ls++)
                {
                    candidate = current + d * alpha;
                    Clamp(candidate, minimumParameters, maximumParameters);

                    candidateVal = diffObjective.Value(function.Bind(candidate));
                    if (candidateVal <= currentVal + c1 * alpha * slope)
                    {
                        stepAccepted = true;
                        break;
                    }
                    alpha *= rho;
                }

                if (!stepAccepted || (candidate - current).Norm() < 1e-12)
                {
                    // Рестарт в антиградиент
                    d = g * -1.0;
                    alpha = InitialStep;
                    candidate = current + d * alpha;
                    Clamp(candidate, minimumParameters, maximumParameters);
                    candidateVal = diffObjective.Value(function.Bind(candidate));
                    if (candidateVal >= currentVal)
                    {
                        break; // Дальнейшее улучшение невозможно
                    }
                }

                current = candidate;
                currentVal = candidateVal;

                boundFunc = function.Bind(current);
                var gNew = new Vector(diffObjective.Gradient(boundFunc));

                // Коэффициент Полака-Рибьера с рестартом: beta = max(0, gNew^T * (gNew - g) / ||g||^2)
                var y = gNew - g;
                double beta = gNorm > 1e-15 ? gNew.Dot(y) / (gNorm * gNorm) : 0.0;
                if (beta < 0 || (k + 1) % p == 0)
                {
                    beta = 0.0; // Рестарт
                }

                var dNew = (gNew * -1.0) + (d * beta);

                // Проверка условия спуска
                if (dNew.Dot(gNew) >= 0)
                {
                    dNew = gNew * -1.0;
                }

                g = gNew;
                d = dNew;
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
