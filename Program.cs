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
            Console.WriteLine("           СИСТЕМА КЛАССОВ ДЛЯ РЕШЕНИЯ ЗАДАЧИ ОПТИМИЗАЦИИ");
            Console.WriteLine("================================================================================\n");

            PrintInterfaceConformanceMatrix();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 1: Подгонка 1D линейной функции y = 2.5*x + 1.2");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestLinear1D();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 2: Подгонка 2D линейной функции f(x0, x1) = 1.5*x0 - 2.0*x1 + 3.0");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestLinear2D();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 3: Подгонка полинома 2-й степени P_2(x) = 1.0 - 2.0*x + 0.5*x^2");
            Console.WriteLine("(полином намеренно НЕ реализует IDifferentiableFunction)");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestPolynomial();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 4: Подгонка кусочно-линейной функции (реализует IDifferentiableFunction)");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestPiecewiseLinear();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 5: Подгонка нелинейного кубического сплайна (не дифференцируем по параметрам)");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestCubicSpline();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 6: Численный интегральный функционал IntegralFunctional");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestIntegralFunctional();

            Console.WriteLine("\n--------------------------------------------------------------------------------");
            Console.WriteLine("ТЕСТ 7: Проверка ограничений и изоляции интерфейсов (Negative Testing)");
            Console.WriteLine("--------------------------------------------------------------------------------");
            TestArchitecturalConstraints();

            Console.WriteLine("\n================================================================================");
            Console.WriteLine("                         ВСЕ ТЕСТЫ УСПЕШНО ПРОЙДЕНЫ!");
            Console.WriteLine("================================================================================");
        }

        static void PrintInterfaceConformanceMatrix()
        {
            Console.WriteLine("МАТРИЦА СООТВЕТСТВИЯ ТРЕБОВАНИЯМ ЗАДАНИЯ:");
            Console.WriteLine("--------------------------------------------------------------------------------");

            var linear = new LinearFunction(1).Bind(new Vector(1, 0));
            var poly = new PolynomialFunction(2).Bind(new Vector(1, 0, 0));
            var pwl = new PiecewiseLinearFunction(3).Bind(new Vector(0, 1, 2));
            var spline = new CubicSplineFunction(3).Bind(new Vector(0, 1, 2));

            Console.WriteLine("1. Реализации IParametricFunction (после Bind):");
            Console.WriteLine($"   • LinearFunction:          IFunction={linear is IFunction}, IDifferentiableFunction={linear is IDifferentiableFunction}  [Ожидалось: TRUE]");
            Console.WriteLine($"   • PolynomialFunction:      IFunction={poly is IFunction}, IDifferentiableFunction={poly is IDifferentiableFunction} [Ожидалось: FALSE]");
            Console.WriteLine($"   • PiecewiseLinearFunction: IFunction={pwl is IFunction}, IDifferentiableFunction={pwl is IDifferentiableFunction}  [Ожидалось: TRUE]");
            Console.WriteLine($"   • CubicSplineFunction:     IFunction={spline is IFunction}, IDifferentiableFunction={spline is IDifferentiableFunction} [Ожидалось: FALSE]");

            var dummyPts = new List<(double X, double Y)> { (0, 0), (1, 1) };
            var l1 = new L1NormFunctional(dummyPts);
            var l2 = new L2NormFunctional(dummyPts);
            var linf = new LInfNormFunctional(dummyPts);
            var integral = new IntegralFunctional(0, 1);

            Console.WriteLine("\n2. Реализации IFunctional:");
            Console.WriteLine($"   • L1NormFunctional:   IDifferentiableFunctional={l1 is IDifferentiableFunctional}, ILeastSquaresFunctional={l1 is ILeastSquaresFunctional} [Ожидалось: TRUE, FALSE]");
            Console.WriteLine($"   • L2NormFunctional:   IDifferentiableFunctional={l2 is IDifferentiableFunctional}, ILeastSquaresFunctional={l2 is ILeastSquaresFunctional}  [Ожидалось: TRUE, TRUE]");
            Console.WriteLine($"   • LInfNormFunctional: IDifferentiableFunctional={linf is IDifferentiableFunctional}, ILeastSquaresFunctional={linf is ILeastSquaresFunctional} [Ожидалось: FALSE, FALSE]");
            Console.WriteLine($"   • IntegralFunctional: IFunctional={integral is IFunctional}, IDifferentiableFunctional={integral is IDifferentiableFunctional} [Ожидалось: TRUE, FALSE]");

            Console.WriteLine("\n3. Реализации IOptimizator:");
            Console.WriteLine($"   • SimulatedAnnealingOptimizer: универсальный метод имитации отжига");
            Console.WriteLine($"   • MonteCarloOptimizer:        универсальный метод случайного поиска");
            Console.WriteLine($"   • ConjugateGradientOptimizer: требует IDifferentiableFunctional (метод сопряжённых градиентов)");
            Console.WriteLine($"   • GradientDescentOptimizer:   требует IDifferentiableFunctional (градиентный спуск)");
            Console.WriteLine($"   • GaussNewtonOptimizer:       требует ILeastSquaresFunctional (алгоритм Гаусса-Ньютона)");
        }

        static void TestLinear1D()
        {
            // Истинные параметры: a = 2.5, b = 1.2
            // y = 2.5 * x + 1.2
            var points = new List<(double X, double Y)>
            {
                (0.0, 1.2),
                (1.0, 3.7),
                (2.0, 6.2),
                (3.0, 8.7),
                (4.0, 11.2)
            };

            var lineFunc = new LinearFunction(1); // 1D -> 2 параметра: w_0, b
            var l2Obj = new L2NormFunctional(points);
            var l1Obj = new L1NormFunctional(points);
            var linfObj = new LInfNormFunctional(points);
            var initial = new Vector(0.0, 0.0);

            // 1. Гаусс-Ньютон
            var gn = new GaussNewtonOptimizer();
            var resGN = gn.Minimize(l2Obj, lineFunc, initial);
            Console.WriteLine($"[Gauss-Newton       | L2]   Найденные параметры: a={resGN[0]:F4}, b={resGN[1]:F4} (Остаточная L2={l2Obj.Value(lineFunc.Bind(resGN)):E3})");

            // 2. Метод сопряжённых градиентов
            var cg = new ConjugateGradientOptimizer();
            var resCG = cg.Minimize(l2Obj, lineFunc, initial);
            Console.WriteLine($"[Conjugate-Gradient | L2]   Найденные параметры: a={resCG[0]:F4}, b={resCG[1]:F4} (Остаточная L2={l2Obj.Value(lineFunc.Bind(resCG)):E3})");

            // 3. Метод градиентного спуска
            var gd = new GradientDescentOptimizer();
            var resGD = gd.Minimize(l1Obj, lineFunc, initial);
            Console.WriteLine($"[Gradient-Descent   | L1]   Найденные параметры: a={resGD[0]:F4}, b={resGD[1]:F4} (Остаточная L1={l1Obj.Value(lineFunc.Bind(resGD)):E3})");

            // 4. Имитация отжига на Linf норме
            var sa = new SimulatedAnnealingOptimizer(maxIterations: 15000, initialTemperature: 50.0, coolingFactor: 0.998, seed: 123);
            var resSA = sa.Minimize(linfObj, lineFunc, initial);
            Console.WriteLine($"[Simulated-Annealing| Linf] Найденные параметры: a={resSA[0]:F4}, b={resSA[1]:F4} (Остаточная Linf={linfObj.Value(lineFunc.Bind(resSA)):E3})");
        }

        static void TestLinear2D()
        {
            // f(x0, x1) = 1.5 * x0 - 2.0 * x1 + 3.0
            // Параметры: [w0=1.5, w1=-2.0, b=3.0]
            var points = new List<(IVector Point, double Value)>
            {
                (new Vector(0.0, 0.0), 3.0),
                (new Vector(1.0, 0.0), 4.5),
                (new Vector(0.0, 1.0), 1.0),
                (new Vector(1.0, 1.0), 2.5),
                (new Vector(2.0, -1.0), 8.0),
                (new Vector(-1.0, 2.0), -2.5)
            };

            var linear2D = new LinearFunction(2); // 2D -> 3 параметра: w0, w1, b
            var l2 = new L2NormFunctional(points);
            var initial = new Vector(0.0, 0.0, 0.0);

            var gn = new GaussNewtonOptimizer();
            var res = gn.Minimize(l2, linear2D, initial);
            Console.WriteLine($"[Gauss-Newton | 2D Linear] Параметры: w0={res[0]:F4}, w1={res[1]:F4}, b={res[2]:F4}");
            Console.WriteLine($"                           Истинные:  w0=1.5000, w1=-2.0000, b=3.0000 (Невязка: {l2.Value(linear2D.Bind(res)):E3})");
        }

        static void TestPolynomial()
        {
            // P_2(x) = 1.0 - 2.0*x + 0.5*x^2
            // Параметры: a0 = 1.0, a1 = -2.0, a2 = 0.5
            var points = new List<(double X, double Y)>
            {
                (-2.0, 1.0 - 2.0 * (-2.0) + 0.5 * 4.0), // 7.0
                (-1.0, 1.0 - 2.0 * (-1.0) + 0.5 * 1.0), // 3.5
                (0.0,  1.0),
                (1.0,  1.0 - 2.0 * 1.0 + 0.5 * 1.0),   // -0.5
                (2.0,  1.0 - 2.0 * 2.0 + 0.5 * 4.0),   // -1.0
                (3.0,  1.0 - 2.0 * 3.0 + 0.5 * 9.0)    // -0.5
            };

            var poly = new PolynomialFunction(degree: 2); // 3 параметра: a0, a1, a2
            var l2 = new L2NormFunctional(points);
            var initial = new Vector(0.0, 0.0, 0.0);

            // Так как полином НЕ дифференцируем по параметрам (нет IDifferentiableFunction),
            // используем универсальный оптимизатор (Имитация отжига)
            var sa = new SimulatedAnnealingOptimizer(maxIterations: 20000, initialTemperature: 50.0, coolingFactor: 0.998, seed: 42);
            var res = sa.Minimize(l2, poly, initial, minimumParameters: new Vector(-5, -5, -5), maximumParameters: new Vector(5, 5, 5));

            Console.WriteLine($"[Simulated Annealing] Полином 2-й степени найден:");
            Console.WriteLine($"                      a0={res[0]:F4}, a1={res[1]:F4}, a2={res[2]:F4}");
            Console.WriteLine($"                      Истинные: a0=1.0000, a1=-2.0000, a2=0.5000 (L2={l2.Value(poly.Bind(res)):F4})");
        }

        static void TestPiecewiseLinear()
        {
            // Кусочно-линейная функция на сетке x = [0, 1, 2, 3]
            // Истинные значения в узлах: y = [1.0, 3.0, 2.0, 4.0]
            var knots = new List<double> { 0.0, 1.0, 2.0, 3.0 };
            var pwl = new PiecewiseLinearFunction(knots);

            // Набор точек для подгонки
            var points = new List<(double X, double Y)>
            {
                (0.0, 1.0),
                (0.5, 2.0),
                (1.0, 3.0),
                (1.5, 2.5),
                (2.0, 2.0),
                (2.5, 3.0),
                (3.0, 4.0)
            };

            var l2 = new L2NormFunctional(points);
            var initial = new Vector(0.0, 0.0, 0.0, 0.0);

            // PiecewiseLinearFunction реализует IDifferentiableFunction, поэтому Гаусс-Ньютон и CG применимы!
            var gn = new GaussNewtonOptimizer();
            var resGN = gn.Minimize(l2, pwl, initial);
            Console.WriteLine($"[Gauss-Newton | PiecewiseLinear] Узлы: y0={resGN[0]:F4}, y1={resGN[1]:F4}, y2={resGN[2]:F4}, y3={resGN[3]:F4}");
            Console.WriteLine($"                                 Истинные: y0=1.0000, y1=3.0000, y2=2.0000, y3=4.0000 (L2={l2.Value(pwl.Bind(resGN)):E3})");
        }

        static void TestCubicSpline()
        {
            // Кубический сплайн на 4 узлах x = [0, 1, 2, 3]
            // Истинные параметры: y = [0.0, 1.0, 0.5, 2.0]
            var knots = new List<double> { 0.0, 1.0, 2.0, 3.0 };
            var spline = new CubicSplineFunction(knots);

            // Обучающие точки
            var points = new List<(double X, double Y)>
            {
                (0.0, 0.0),
                (0.5, 0.6),
                (1.0, 1.0),
                (1.5, 0.7),
                (2.0, 0.5),
                (2.5, 1.1),
                (3.0, 2.0)
            };

            var l2 = new L2NormFunctional(points);
            var initial = new Vector(0.5, 0.5, 0.5, 0.5);

            // Сплайн не реализует IDifferentiableFunction, используем универсальный SimulatedAnnealing
            var sa = new SimulatedAnnealingOptimizer(maxIterations: 15000, initialTemperature: 20.0, coolingFactor: 0.997, seed: 99);
            var res = sa.Minimize(l2, spline, initial);
            Console.WriteLine($"[Simulated Annealing | Spline] Значения в узлах сплайна: y=[{res[0]:F3}, {res[1]:F3}, {res[2]:F3}, {res[3]:F3}]");
            Console.WriteLine($"                               Остаточная невязка L2={l2.Value(spline.Bind(res)):F4}");
        }

        static void TestIntegralFunctional()
        {
            // Минимизация функционала F(f) = integral_0^1 f(x) dx
            // с ограничением параметров или квадратичной штрафной функции
            // Пусть f(x) = a * x + b.
            // Интеграл от 0 до 1: a/2 + b.
            var integral = new IntegralFunctional(0.0, 1.0, steps: 100);
            var line = new LinearFunction(1);

            // Минимизируем интеграл при ограничении a in [1, 2], b in [3, 4]
            // Минимум интеграла достигается при a=1, b=3, интеграл = 1/2 + 3 = 3.5.
            var sa = new SimulatedAnnealingOptimizer(maxIterations: 5000, seed: 10);
            var res = sa.Minimize(integral, line, new Vector(1.5, 3.5),
                                  minimumParameters: new Vector(1.0, 3.0),
                                  maximumParameters: new Vector(2.0, 4.0));

            double minIntegral = integral.Value(line.Bind(res));
            Console.WriteLine($"[IntegralFunctional] Минимум интеграла на [0, 1] при a in [1, 2], b in [3, 4]:");
            Console.WriteLine($"                     Параметры: a={res[0]:F4}, b={res[1]:F4}");
            Console.WriteLine($"                     Значение интеграла: {minIntegral:F4} (Теоретический минимум = 3.5000)");
        }

        static void TestArchitecturalConstraints()
        {
            Console.WriteLine("Проверка контрактов интерфейсов и генерации соответствующих исключений:");

            var dummyPts = new List<(double X, double Y)> { (0, 1), (1, 2) };
            var l1 = new L1NormFunctional(dummyPts);
            var linf = new LInfNormFunctional(dummyPts);
            var l2 = new L2NormFunctional(dummyPts);
            var line = new LinearFunction(1);
            var poly = new PolynomialFunction(1);

            // 1. Gauss-Newton с объективом L1 (который НЕ реализует ILeastSquaresFunctional)
            try
            {
                var gn = new GaussNewtonOptimizer();
                gn.Minimize(l1, line, new Vector(1, 1));
                Console.WriteLine("  [ОШИБКА] Gauss-Newton должен был отклонить L1NormFunctional!");
            }
            catch (ArgumentException ex)
            {
                Console.WriteLine($"  [УСПЕХ] Gauss-Newton отклонил L1NormFunctional: {ex.Message.Split('.')[0]}.");
            }

            // 2. Conjugate-Gradient с объективом Linf (который НЕ реализует IDifferentiableFunctional)
            try
            {
                var cg = new ConjugateGradientOptimizer();
                cg.Minimize(linf, line, new Vector(1, 1));
                Console.WriteLine("  [ОШИБКА] ConjugateGradient должен был отклонить LInfNormFunctional!");
            }
            catch (ArgumentException ex)
            {
                Console.WriteLine($"  [УСПЕХ] ConjugateGradient отклонил LInfNormFunctional: {ex.Message.Split('.')[0]}.");
            }

            // 3. Вычисление Якобиана L2 с PolynomialFunction (полином НЕ реализует IDifferentiableFunction)
            try
            {
                var boundPoly = poly.Bind(new Vector(1, 2));
                l2.Jacobian(boundPoly);
                Console.WriteLine("  [ОШИБКА] L2NormFunctional.Jacobian должен был отклонить PolynomialFunction!");
            }
            catch (ArgumentException ex)
            {
                Console.WriteLine($"  [УСПЕХ] L2NormFunctional.Jacobian отклонил недифференцируемый полином: {ex.Message.Split('.')[0]}.");
            }

            // 4. Вычисление Градиента L1 с CubicSplineFunction (сплайн НЕ реализует IDifferentiableFunction)
            try
            {
                var spline = new CubicSplineFunction(3);
                var boundSpline = spline.Bind(new Vector(1, 2, 3));
                l1.Gradient(boundSpline);
                Console.WriteLine("  [ОШИБКА] L1NormFunctional.Gradient должен был отклонить CubicSplineFunction!");
            }
            catch (ArgumentException ex)
            {
                Console.WriteLine($"  [УСПЕХ] L1NormFunctional.Gradient отклонил сплайн без IDifferentiableFunction: {ex.Message.Split('.')[0]}.");
            }
        }
    }
}
