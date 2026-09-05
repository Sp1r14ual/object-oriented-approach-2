# Система классов для решения задачи оптимизации

Проект реализует гибкую архитектуру для решения задач оптимизации параметров функций на основе системы фиксированных интерфейсов в соответствии с требованиями задания.

---

## 1. Архитектура и обязательные требования

### Принцип взаимодействия через интерфейсы
Классы оптимизаторов, функционалов и функций взаимодействуют **исключительно** через фиксированные интерфейсы, не привязываясь к конкретным типам.
- Оптимизаторы работают через интерфейсы `IFunctional`, `IDifferentiableFunctional`, `ILeastSquaresFunctional`, `IParametricFunction`.
- Функционалы работают через `IFunction` и `IDifferentiableFunction`.
- Векторы и матрицы абстрагированы через `IVector` и `IMatrix`.

---

## 2. Фиксированные интерфейсы

```csharp
public interface IVector : IList<double> { }
public interface IMatrix : IList<IList<double>> { }

namespace Functions
{
    public interface IParametricFunction
    {
        IFunction Bind(IVector parameters);
    }

    public interface IFunction
    {
        double Value(IVector point);
    }

    public interface IDifferentiableFunction : IFunction
    {
        // Вычисляет градиент по параметрам исходной IParametricFunction
        IVector Gradient(IVector point);
    }
}

namespace Functionals
{
    public interface IFunctional
    {
        double Value(IFunction function);
    }

    public interface IDifferentiableFunctional : IFunctional
    {
        IVector Gradient(IFunction function);
    }

    public interface ILeastSquaresFunctional : IFunctional
    {
        IVector Residual(IFunction function);
        IMatrix Jacobian(IFunction function);
    }
}

public interface IOptimizator
{
    IVector Minimize(IFunctional objective,
                     IParametricFunction function,
                     IVector initialParameters,
                     IVector minimumParameters = default,
                     IVector maximumParameters = default);
}
```

---

## 3. Матрица реализаций

| Компонент | Класс | Реализуемые интерфейсы | Описание / Особенности |
|---|---|---|---|
| **Функция 1** | `LinearFunction` | `IParametricFunction` $\to$ `IDifferentiableFunction` | Линейная $n$-мерная функция $f(x) = \sum w_i x_i + b$. Число параметров: $n+1$. Дифференцируема по параметрам. |
| **Функция 2** | `PolynomialFunction` | `IParametricFunction` $\to$ `IFunction` | Полином $n$-й степени в 1D $P_n(x) = \sum a_k x^k$. Число параметров: $n+1$. **Не реализует** `IDifferentiableFunction`. |
| **Функция 3** | `PiecewiseLinearFunction` | `IParametricFunction` $\to$ `IDifferentiableFunction` | Кусочно-линейная функция на сетке узлов. Параметры — значения в узлах. **Реализует** `IDifferentiableFunction`. |
| **Функция 4** | `CubicSplineFunction` | `IParametricFunction` $\to$ `IFunction` | Нелинейный кубический сплайн. **Не реализует** `IDifferentiableFunction`. |
| **Функционал 1** | `L1NormFunctional` | `IDifferentiableFunctional` | $L_1$-норма разности $\sum \|f(x_i) - y_i\|$. **Не реализует** `ILeastSquaresFunctional`. |
| **Функционал 2** | `L2NormFunctional` | `IDifferentiableFunctional`, `ILeastSquaresFunctional` | $L_2$-норма разности $\sqrt{\sum (f(x_i) - y_i)^2}$. Вектор невязок и матрица Якоби. |
| **Функционал 3** | `LInfNormFunctional` | `IFunctional` | $L_\infty$-норма $\max \|f(x_i) - y_i\|$. **Не реализует** `IDifferentiableFunctional` и `ILeastSquaresFunctional`. |
| **Функционал 4** | `IntegralFunctional` | `IFunctional` | Численный интеграл по области $\int_\Omega f(x) dx$ (формула Симпсона 1/3 в 1D, кратный интеграл в $n$D). |
| **Оптимизатор 1** | `SimulatedAnnealingOptimizer` | `IOptimizator` | **Универсальный**: алгоритм имитации отжига (работает с любыми функционалами и функциями). |
| **Оптимизатор 1'** | `MonteCarloOptimizer` | `IOptimizator` | **Универсальный**: метод случайного поиска Монте-Карло. |
| **Оптимизатор 2** | `ConjugateGradientOptimizer` | `IOptimizator` | **Требует `IDifferentiableFunctional`**: метод нелинейных сопряжённых градиентов Полака-Рибьера с одномерным поиском Армихо. |
| **Оптимизатор 2'** | `GradientDescentOptimizer` | `IOptimizator` | **Требует `IDifferentiableFunctional`**: градиентный спуск с дроблением шага. |
| **Оптимизатор 3** | `GaussNewtonOptimizer` | `IOptimizator` | **Требует `ILeastSquaresFunctional`**: алгоритм Гаусса-Ньютона с решением нормальных уравнений $(J^T J + \lambda I)\Delta = -J^T r$. |

---

## 4. Запуск и тестирование

Сборка и запуск программы:
```bash
dotnet run
```

Программа выполняет комплексный набор тестов:
1. Вывод матрицы соответствия контрактам интерфейсов;
2. Подгонка 1D линейной функции всеми тремя группами оптимизаторов;
3. Подгонка 2D линейной функции методом Гаусса-Ньютона;
4. Подгонка полинома методом имитации отжига;
5. Подгонка кусочно-линейной функции методом Гаусса-Ньютона;
6. Подгонка нелинейного кубического сплайна;
7. Минимизация интегрального функционала;
8. Негативные тесты: проверка того, что оптимизаторы и функционалы корректно отклоняют неподдерживаемые интерфейсы с понятным `ArgumentException`.