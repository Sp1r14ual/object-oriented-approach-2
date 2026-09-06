using System;
using System.Collections.Generic;
using Functions;
using Functionals;
using LinearAlgebra;
using Optimizers;

namespace OptimizationApp
{
    class Program
    {
        static void Main(string[] args)
        {
            var knots = new List<double> { 0.0, 1.0, 2.0, 3.0, 4.0 };
            IParametricFunction parametricFunc = new PiecewiseLinearFunction(knots);

            IVector initParams = new Vector(0.0, 0.0, 0.0, 0.0, 0.0);
            IFunction boundFunc = parametricFunc.Bind(initParams);

            var targetPoints = new List<(IVector Point, double Value)>
            {
                (new Vector(0.0), 1.0),
                (new Vector(0.5), 2.0),
                (new Vector(1.0), 3.0),
                (new Vector(1.5), 2.5),
                (new Vector(2.0), 2.0),
                (new Vector(2.5), 3.5),
                (new Vector(3.0), 5.0),
                (new Vector(3.5), 4.5),
                (new Vector(4.0), 4.0)
            };

            IFunctional functional = new L2NormFunctional(targetPoints);
            IOptimizator optimizer = new GaussNewtonOptimizer(maxIterations: 50, tolerance: 1e-10);

            Console.WriteLine("Решение задачи минимизации L2-норма:");

            double initialError = functional.Value(boundFunc);
            Console.WriteLine($"Начальная L2-невязка с параметрами {initParams}: {initialError:F6}");

            IVector optimizedParams = optimizer.Minimize(
                objective: functional,
                function: parametricFunc,
                initialParameters: initParams
            );

            IFunction optimizedFunc = parametricFunc.Bind(optimizedParams);
            double finalError = functional.Value(optimizedFunc);

            Console.WriteLine("Результаты:");
            Console.WriteLine($"Найденные значения в узлах сетки x = [0, 1, 2, 3, 4]:");
            for (int i = 0; i < knots.Count; i++)
            {
                Console.WriteLine($"  x = {knots[i]:F1} -> y = {optimizedParams[i]:F4}");
            }

            Console.WriteLine($"Итоговая L2-невязка: {finalError:E6}");

            Console.WriteLine("Проверка значений в контрольных точках:");
            foreach (var (pt, expected) in targetPoints)
            {
                double actual = optimizedFunc.Value(pt);
                Console.WriteLine($"  f({pt[0]:F1}) = {actual:F4} (ожидалось: {expected:F4}, ошибка: {Math.Abs(actual - expected):E2})");
            }

            Console.WriteLine("Тест завершён");
        }
    }
}
