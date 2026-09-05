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
            Console.OutputEncoding = System.Text.Encoding.UTF8;
            Console.WriteLine("================================================================================");
            Console.WriteLine("  СИСТЕМА КЛАССОВ ОПТИМИЗАЦИИ: ВЫБРАННЫЙ ВАРИАНТ РЕАЛИЗАЦИЙ");
            Console.WriteLine("================================================================================\n");

            Console.WriteLine("В Example.cs уже реализованы:");
            Console.WriteLine("  - IParametricFunction: LineFunction (линейная функция)");
            Console.WriteLine("  - IFunctional:         MyFunctional (базовый функционал без дифференцирования)");
            Console.WriteLine("  - IOptimizator:        MinimizerMonteCarlo (метод Монте-Карло)\n");

            Console.WriteLine("Выбранные реализации из задания (которых не было в Example.cs):");
            Console.WriteLine("  1. IParametricFunction -> PiecewiseLinearFunction (кусочно-линейная, реализует IDifferentiableFunction)");
            Console.WriteLine("  2. IFunctional         -> L2NormFunctional (l2 норма, реализует IDifferentiableFunctional и ILeastSquaresFunctional)");
            Console.WriteLine("  3. IOptimizator        -> GaussNewtonOptimizer (алгоритм Гаусса-Ньютона, требующий ILeastSquaresFunctional)\n");

            Console.WriteLine("--------------------------------------------------------------------------------");
            Console.WriteLine("ПРОВЕРКА КОНТРАКТОВ ИНТЕРФЕЙСОВ:");
            Console.WriteLine("--------------------------------------------------------------------------------");

            var knots = new List<double> { 0.0, 1.0, 2.0, 3.0, 4.0 };
            IParametricFunction parametricFunc = new PiecewiseLinearFunction(knots);

            IVector initParams = new Vector(0.0, 0.0, 0.0, 0.0, 0.0);
            IFunction boundFunc = parametricFunc.Bind(initParams);

            Console.WriteLine($"• PiecewiseLinearFunction.Bind возвращает IFunction:                 {boundFunc is IFunction} (True)");
            Console.WriteLine($"• PiecewiseLinearFunction.Bind возвращает IDifferentiableFunction:  {boundFunc is IDifferentiableFunction} (True)");

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
            Console.WriteLine($"• L2NormFunctional реализует IFunctional:                          {functional is IFunctional} (True)");
            Console.WriteLine($"• L2NormFunctional реализует IDifferentiableFunctional:            {functional is IDifferentiableFunctional} (True)");
            Console.WriteLine($"• L2NormFunctional реализует ILeastSquaresFunctional:              {functional is ILeastSquaresFunctional} (True)");

            IOptimizator optimizer = new GaussNewtonOptimizer(maxIterations: 50, tolerance: 1e-10);
            Console.WriteLine($"• GaussNewtonOptimizer реализует IOptimizator:                     {optimizer is IOptimizator} (True)");

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("РЕШЕНИЕ ЗАДАЧИ ОПТИМИЗАЦИИ (ВЗАИМОДЕЙСТВИЕ ТОЛЬКО ЧЕРЕЗ ИНТЕРФЕЙСЫ):");
            Console.WriteLine("--------------------------------------------------------------------------------");

            double initialError = functional.Value(boundFunc);
            Console.WriteLine($"Начальная L2-невязка с параметрами {initParams}: {initialError:F6}");

            // Оптимизатор взаимодействует только через IFunctional, IParametricFunction, IVector
            IVector optimizedParams = optimizer.Minimize(
                objective: functional,
                function: parametricFunc,
                initialParameters: initParams
            );

            IFunction optimizedFunc = parametricFunc.Bind(optimizedParams);
            double finalError = functional.Value(optimizedFunc);

            Console.WriteLine($"\nОптимизация завершена!");
            Console.WriteLine($"Найденные значения в узлах сетки x = [0, 1, 2, 3, 4]:");
            for (int i = 0; i < knots.Count; i++)
            {
                Console.WriteLine($"  x = {knots[i]:F1} -> y = {optimizedParams[i]:F4}");
            }

            Console.WriteLine($"\nИтоговая L2-невязка: {finalError:E6}");

            Console.WriteLine("\nПроверка значений в контрольных точках:");
            foreach (var (pt, expected) in targetPoints)
            {
                double actual = optimizedFunc.Value(pt);
                Console.WriteLine($"  f({pt[0]:F1}) = {actual:F4} (ожидалось: {expected:F4}, ошибка: {Math.Abs(actual - expected):E2})");
            }

            Console.WriteLine("\n================================================================================");
            Console.WriteLine("                      ТЕСТ УСПЕШНО ЗАВЕРШЁН!");
            Console.WriteLine("================================================================================");
        }
    }
}
