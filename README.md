# Система классов для решения задачи оптимизации

Проект реализует систему классов для решения задачи параметрической оптимизации на базе фиксированных интерфейсов.

---

## 1. Архитектура и обязательные требования

- **Интерфейсы неизменяемы**: Все интерфейсы объявлены в точном соответствии с заданием.
- **Взаимодействие строго через интерфейсы**:
  - Оптимизатор `IOptimizator` принимает только `IFunctional`, `IParametricFunction`, `IVector`.
  - Функционал `IFunctional` принимает только `IFunction`.
  - При необходимости вычисления градиента или матрицы Якоби функционалы и оптимизаторы запрашивают специализированные интерфейсы: `IDifferentiableFunction`, `IDifferentiableFunctional`, `ILeastSquaresFunctional`.

---

## 2. Реализации в проекте

### Из Example.cs (исходные базовые реализации):
1. **`LineFunction : IParametricFunction`** — одномерная линейная функция $f(x) = ax + b$ (реализует только `IFunction`).
2. **`MyFunctional : IFunctional`** — функционал суммы квадратов ошибок.
3. **`MinimizerMonteCarlo : IOptimizator`** — универсальный метод случайного поиска Монте-Карло.

### Выбранные реализации из задания (которых не было в Example.cs):
1. **`PiecewiseLinearFunction : IParametricFunction`** (Кусочно-линейная функция):
   - Задаётся сеткой узлов $x_0 < x_1 < \dots < x_{k-1}$.
   - Параметрами являются значения функции в узлах: $(y_0, y_1, \dots, y_{k-1})$.
   - Метод `Bind` возвращает реализацию **`IDifferentiableFunction`**.
   - Аналитический градиент по параметрам: $\frac{\partial f}{\partial y_i} = 1 - t$, $\frac{\partial f}{\partial y_{i+1}} = t$.

2. **`L2NormFunctional : IDifferentiableFunctional, ILeastSquaresFunctional`** ($l_2$ норма разности):
   - Реализует евклидову норму $\|r\|_2 = \sqrt{\sum (f(x_i) - y_i)^2}$.
   - Реализует **`ILeastSquaresFunctional`**:
     - `Residual(IFunction function)`: вектор невязок $r_i = f(x_i) - y_i$.
     - `Jacobian(IFunction function)`: матрица Якоби $J_{i, j} = \frac{\partial f(x_i)}{\partial \theta_j}$, вычисляемая через `IDifferentiableFunction.Gradient`.
   - Реализует **`IDifferentiableFunctional`**: градиент $\nabla F = \frac{J^T r}{\|r\|_2}$.

3. **`GaussNewtonOptimizer : IOptimizator`** (Алгоритм Гаусса-Ньютона):
   - Требует реализации **`ILeastSquaresFunctional`**.
   - На каждой итерации решает систему нормальных уравнений:
     $$(J^T J + \lambda I)\Delta \theta = -J^T r$$
   - Выполняет одномерный поиск шага и обновляет вектор параметров с квадратичной скоростью сходимости.

---

## 3. Структура файлов

```
object-oriented-approach-2/
├── Interfaces/
│   └── Interfaces.cs             # Фиксированные интерфейсы (IVector, IMatrix, Functions, Functionals, IOptimizator)
├── LinearAlgebra/
│   ├── Vector.cs                 # Реализация IVector (операции, нормы, скалярное произведение)
│   └── Matrix.cs                 # Реализация IMatrix (СЛАУ по Гауссу, нормальные уравнения)
├── Functions/
│   └── PiecewiseLinearFunction.cs# Кусочно-линейная функция (IDifferentiableFunction)
├── Functionals/
│   └── L2NormFunctional.cs       # L2-норма (IDifferentiableFunctional, ILeastSquaresFunctional)
├── Optimizers/
│   └── GaussNewtonOptimizer.cs   # Алгоритм Гаусса-Ньютона (ILeastSquaresFunctional)
├── Example.cs                    # Исходный пример из задания
├── Program.cs                    # Демонстрационный запуск и проверка контрактов
├── README.md                     # Документация проекта
└── OptimizationApp.csproj        # Файл проекта .NET
```

---

## 4. Запуск

```bash
dotnet run
```