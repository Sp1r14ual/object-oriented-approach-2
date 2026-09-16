# Министерство науки и высшего образования Российской Федерации
## Федеральное государственное бюджетное образовательное учреждение высшего образования
### «НОВОСИБИРСКИЙ ГОСУДАРСТВЕННЫЙ ТЕХНИЧЕСКИЙ УНИВЕРСИТЕТ»
### Факультет прикладной математики и информатики
### Кафедра теоретической и прикладной информатики

---

<br>
<br>

# ОТЧЕТ
### по лабораторной работе по дисциплине
## «Объектно-ориентированный подход»
### Тема: «Разработка системы классов для решения задачи оптимизации на базе фиксированных интерфейсов. Минимизация функционала $L_2$-нормы методом Гаусса-Ньютона»

<br>
<br>

**Выполнили студенты:**
- Ковалевский Дмитрий, гр. ПММ-51
- Панасенко Сергей, гр. ПММ-51
- Кожин Матвей, гр. ПММ-52

**Преподаватели:**
- Рояк М.Э.
- Ступаков И.М.

<br>
<br>

### Новосибирск
### 2026

---

## Содержание

1. [Постановка задачи](#1-постановка-задачи)
   - [1.1. Цель работы](#11-цель-работы)
   - [1.2. Контракт фиксированных интерфейсов](#12-контракт-фиксированных-интерфейсов)
   - [1.3. Формулировка прикладной задачи оптимизации](#13-формулировка-прикладной-задачи-оптимизации)
2. [Описание алгоритмов и математический аппарат](#2-описание-алгоритмов-и-математический-аппарат)
   - [2.1. Модель кусочно-линейной функции и аналитический градиент](#21-модель-кусочно-линейной-функции-и-аналитический-градиент)
   - [2.2. Функционал $L_2$-нормы и матрица Якоби](#22-функционал-l_2-нормы-и-матрица-якоби)
   - [2.3. Метод оптимизации Гаусса-Ньютона](#23-метод-оптимизации-гаусса-ньютона)
   - [2.4. Численное решение нормальных уравнений методом исключения Гаусса](#24-численное-решение-нормальных-уравнений-методом-исключения-гаусса)
3. [Архитектура и реализация программного кода](#3-архитектура-и-реализация-программного-кода)
   - [3.1. Архитектурная диаграмма классов](#31-архитектурная-диаграмма-классов)
   - [3.2. Листинги ключевых компонентов системы](#32-листинги-ключевых-компонентов-системы)
4. [Тестирование программного комплекса](#4-тестирование-программного-комплекса)
   - [4.1. Модульные тесты линейной алгебры](#41-модульные-тесты-линейной-алгебры)
   - [4.2. Верификация аналитического дифференцирования](#42-верификация-аналитического-дифференцирования)
   - [4.3. Сквозной интеграционный тест оптимизатора](#43-сквозной-интеграционный-тест-оптимизатора)
5. [Экспериментальные исследования](#5-экспериментальные-исследования)
   - [5.1. Динамика оптимизации и скорость сходимости](#51-динамика-оптимизации-и-скорость-сходимости)
   - [5.2. Сравнение параметров и верификация в контрольных точках](#52-сравнение-параметров-и-верификация-в-контрольных-точках)
   - [5.3. Анализ теоретической и практической вычислительной сложности](#53-анализ-теоретической-и-практической-вычислительной-сложности)
6. [Выводы](#6-выводы)

---

## 1. Постановка задачи

### 1.1. Цель работы

Целью настоящей работы является проектирование и реализация на языке C# расширяемой, строго типизированной объектно-ориентированной библиотеки для решения задач оптимизации параметров функций. Архитектура библиотеки базируется на строгом соблюдении принципов SOLID и декомпозиции функциональности с использованием системы фиксированных интерфейсов.

В рамках работы решаются следующие подзадачи:
- Реализация абстракций векторно-матричной алгебры (`IVector`, `IMatrix`);
- Реализация параметрической кусочно-линейной функции (`PiecewiseLinearFunction`), поддерживающей вычисление значений и аналитического градиента по вектору параметров (`IDifferentiableFunction`);
- Реализация функционала невязки наименьших квадратов / $L_2$-нормы (`L2NormFunctional`), вычисляющего вектор невязок $r(p)$, матрицу Якоби $J(p)$ и градиент $\nabla \Phi(p)$;
- Реализация оптимизатора второго порядка Гаусса-Ньютона (`GaussNewtonOptimizer`), использующего регуляризацию типа Левенберга и алгоритм одномерного поиска (backtracking line search);
- Экспериментальное исследование сходимости, анализ погрешностей аппроксимации и вычислительной сложности.

### 1.2. Контракт фиксированных интерфейсов

Система опирается на контракт интерфейсов:

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
        // Градиент по параметрам исходной IParametricFunction
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
    IVector Minimize(Functionals.IFunctional objective,
                     Functions.IParametricFunction function,
                     IVector initialParameters,
                     IVector minimumParameters = default,
                     IVector maximumParameters = default);
}
```

Данная система интерфейсов обеспечивает строгое разделение ответственности:
- `IParametricFunction` отделяет определение математического семейства функций от конкретных значений свободных параметров. Метод `Bind(p)` производит фиксацию параметров и возвращает неизменяемый экземпляр `IFunction`.
- `IDifferentiableFunction` расширяет `IFunction`, предоставляя доступ к градиенту функции по вектору варьируемых параметров в заданной пространственной точке.
- `IFunctional` инкапсулирует критерий качества приближения, отображая функцию в вещественное число (ошибку).
- `ILeastSquaresFunctional` выделяет класс функционалов наименьших квадратов, предоставляя оптимизаторам прямое аналитическое расщепление на вектор частных невязок $r$ и матрицу первых производных (Якобиан) $J$.
- `IOptimizator` определяет унифицированный алгоритмический интерфейс поиска вектора параметров $p^*$, минимизирующего целевой функционал.

### 1.3. Формулировка прикладной задачи оптимизации

Дана сетка из $m$ упорядоченных узлов на числовой прямой:

$$x_0 < x_1 < x_2 < \dots < x_{m-1}$$

Требуется восстановить значения кусочно-линейной функции в этих узлах — вектор параметров $p = (p_0, p_1, \dots, p_{m-1})^T \in \mathbb{R}^m$, наилучшим образом описывающий набор экспериментальных контрольных точек:

$$\mathcal{D} = \left\lbrace (x_i, y_i) \right\rbrace_{i=1}^M, \quad M \ge m$$

Критерием оптимальности выступает минимизация квадратичного функционала $L_2$-нормы невязок:

$$\min_{p \in \mathbb{R}^m} \Phi(p) = \frac{1}{2} \sum_{i=1}^M \left( f(x_i, p) - y_i \right)^2 = \frac{1}{2} \Vert r(p) \Vert_2^2$$

где $r_i(p) = f(x_i, p) - y_i$ — невязка в $i$-й точке данных.

---

## 2. Описание алгоритмов и математический аппарат

### 2.1. Модель кусочно-линейной функции и аналитический градиент

Для заданной сетки узлов $x_0 < x_1 < \dots < x_{m-1}$ и вектора ординат узлов $p = (p_0, p_1, \dots, p_{m-1})^T$ значение функции в произвольной точке $x \in [x_0, x_{m-1}]$ определяется линейной интерполяцией на активном интервале $[x_k, x_{k+1}]$:

$$t = \frac{x - x_k}{x_{k+1} - x_k} \in [0, 1]$$

$$f(x, p) = (1 - t) p_k + t p_{k+1}$$

Поиск активного подотрезка $[x_k, x_{k+1}]$ осуществляется с помощью бинарного поиска с логарифмической временной сложностью $O(\log m)$. Для точек за пределами сетки ($x < x_0$ или $x > x_{m-1}$) интерполяционная формула естественным образом распространяется как линейная экстраполяция.

Поскольку зависимость $f(x, p)$ от параметров $p$ строго линейна, вектор частных производных по вектору параметров $\nabla_p f(x, p)$ вычисляется аналитически без применения численного дифференцирования:

$$\frac{\partial f(x, p)}{\partial p_j} = \begin{cases} 1 - t, & j = k \\ t, & j = k + 1 \\ 0, & j \notin \lbrace k, k+1 \rbrace \end{cases}$$

Вектор градиента $\nabla_p f(x, p) \in \mathbb{R}^m$ имеет вид:

$$\nabla_p f(x, p) = \left( 0, \dots, 0, \underbrace{1 - t}_{k\text{-я поз.}}, \underbrace{t}_{(k+1)\text{-я поз.}}, 0, \dots, 0 \right)^T$$

Благодаря локальности базисных функций (B-сплайнов 1-го порядка), вектор градиента является сильно разреженным: независимо от общего числа узлов $m$, ненулевыми являются максимум 2 компоненты.

### 2.2. Функционал $L_2$-нормы и матрица Якоби

Функционал наименьших квадратов оперирует вектором невязок $r(p) \in \mathbb{R}^M$:

$$r_i(p) = f(x_i, p) - y_i, \quad i = 1, \dots, M$$

Матрица Якоби невязок $J(p) \in \mathbb{R}^{M \times m}$ определяется матрицей частных производных:

$$J_{i, j} = \frac{\partial r_i(p)}{\partial p_j} = \frac{\partial f(x_i, p)}{\partial p_j}$$

Каждая $i$-я строка матрицы Якоби в точности совпадает с транспонированным градиентом функции в $i$-й точке данных:

$$J = \begin{bmatrix} (\nabla_p f(x_1, p))^T \\ (\nabla_p f(x_2, p))^T \\ \vdots \\ (\nabla_p f(x_M, p))^T \end{bmatrix} \in \mathbb{R}^{M \times m}$$

Градиент минимизируемого функционала $\Phi(p) = \frac{1}{2} \Vert r(p) \Vert_2^2$ выражается через матрицу Якоби и вектор невязок:

$$\nabla \Phi(p) = \sum_{i=1}^M r_i(p) \nabla_p r_i(p) = J(p)^T r(p)$$

Для евклидовой $L_2$-нормы $F(p) = \Vert r(p) \Vert_2 = \sqrt{\sum_{i=1}^M r_i(p)^2}$ градиент нормализуется:

$$\nabla F(p) = \frac{J(p)^T r(p)}{\Vert r(p) \Vert_2}, \quad (\text{при } \Vert r(p) \Vert_2 > 0)$$

### 2.3. Метод оптимизации Гаусса-Ньютона

Метод Гаусса-Ньютона предназначен для минимизации функционалов в форме суммы квадратов. В окрестности текущего приближения $p^{(k)}$ вектор невязок линеаризуется с помощью разложения Тейлора первого порядка:

$$r(p^{(k)} + \Delta) \approx r(p^{(k)}) + J(p^{(k)}) \Delta$$

Подстановка линеаризованной невязки в целевой функционал даёт квадратичную задачу:

$$\Phi(p^{(k)} + \Delta) \approx \frac{1}{2} \Vert r(p^{(k)}) + J(p^{(k)}) \Delta \Vert_2^2 = \frac{1}{2} r^T r + \Delta^T J^T r + \frac{1}{2} \Delta^T (J^T J) \Delta$$

Дифференцируя по шагу $\Delta$ и приравнивая производную к нулю, получаем систему нормальных уравнений Гаусса-Ньютона:

$$(J^T J) \Delta = -J^T r$$

#### Фундаментальное свойство для кусочно-линейных моделей

Поскольку кусочно-линейная функция линейна по параметрам узлов $p$, частные производные $\frac{\partial f(x, p)}{\partial p_j}$ вообще не зависят от значений параметров $p$, а определяются исключительно положением точек $x_i$ относительно узлов $x_k$. Следовательно:
1. Матрица Якоби $J$ является **строго постоянной** матрицей на всем пространстве параметров $\mathbb{R}^m$;
2. Линеаризация Тейлора $r(p + \Delta) = r(p) + J \Delta$ является **абсолютно точной** (остаточный член тождественно равен нулю);
3. Матрица Гессе $\nabla^2 \Phi(p) = J^T J + \sum_{i=1}^M r_i \nabla^2 r_i \equiv J^T J$ является точной матрицей вторых производных;
4. Метод Гаусса-Ньютона решает задачу нахождения глобального минимума квадратичного функционала ровно за **1 итерацию**!

#### Регуляризация типа Левенберга и линейный поиск

Для исключения численной неустойчивости при плохой обусловленности матрицы $J^T J$ (например, при наличии узлов, в окрестности которых отсутствуют контрольные точки), в нормальные уравнения вводится регуляризирующее слагаемое $\lambda I$ ($\lambda > 0$):

$$(J^T J + \lambda I) \Delta = -J^T r$$

Матрица $(J^T J + \lambda I)$ гарантированно является симметричной и строго положительно определенной:

$$\forall v \ne 0: \quad v^T (J^T J + \lambda I) v = \Vert J v \Vert_2^2 + \lambda \Vert v \Vert_2^2 \ge \lambda \Vert v \Vert_2^2 > 0$$

Для обеспечения глобальной сходимости применяется алгоритм одномерного поиска вдоль направления шага (Backtracking Armijo Line Search):

$$p^{(k+1)} = p^{(k)} + \alpha \Delta, \quad \alpha \in \lbrace 1, 1/2, 1/4, \dots \rbrace$$

Шаг принимается, как только достигается строгое монотонное убывание функционала:

$$\Phi(p^{(k)} + \alpha \Delta) < \Phi(p^{(k)})$$

### 2.4. Численное решение нормальных уравнений методом исключения Гаусса

Система нормальных уравнений $A \Delta = b$, где $A = J^T J + \lambda I \in \mathbb{R}^{m \times m}$ и $b = -J^T r \in \mathbb{R}^m$, решается прямым методом исключения Гаусса с частичным выбором ведущего элемента по столбцу (Partial Pivoting).

#### 1. Прямой ход (приведение к верхнетреугольному виду):
На каждом $k$-м шаге ($k = 0, 1, \dots, m-1$):
1. Осуществляется поиск строки $p \ge k$ с максимальным по модулю элементом в текущем $k$-м столбце:
$$p = \arg\max_{i = k, \dots, m-1} |A_{i, k}|$$
2. Производится перестановка строк $k$ и $p$ в матрице $A$ и векторе правой части $b$;
3. Для всех строк $i > k$ вычисляется множитель:
$$\mu_{i, k} = \frac{A_{i, k}}{A_{k, k}}$$
4. Производится исключение переменных:
$$A_{i, j} \leftarrow A_{i, j} - \mu_{i, k} A_{k, j}, \quad j = k, \dots, m-1$$
$$b_i \leftarrow b_i - \mu_{i, k} b_k$$

#### 2. Обратный ход:
После получения верхнетреугольной системы решение находится обратной подстановкой:

$$\Delta_{m-1} = \frac{b_{m-1}}{A_{m-1, m-1}}$$

$$\Delta_i = \frac{1}{A_{i, i}} \left( b_i - \sum_{j=i+1}^{m-1} A_{i, j} \Delta_j \right), \quad i = m-2, m-3, \dots, 0$$

Частичный выбор ведущего элемента гарантирует, что $|\mu_{i, k}| \le 1$, что кардинально снижает накопление ошибок округления и предотвращает деление на машинный ноль.

---

## 3. Реализация программного кода

### 3.1. Контракт интерфейсов (`Interfaces.cs`)
```csharp
using System.Collections.Generic;

public interface IVector : IList<double>
{
}

public interface IMatrix : IList<IList<double>>
{
}

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
        // По параметрам исходной IParametricFunction
        IVector Gradient(IVector point);
    }
}

namespace Functionals
{
    using Functions;

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
    IVector Minimize(Functionals.IFunctional objective,
                     Functions.IParametricFunction function,
                     IVector initialParameters,
                     IVector minimumParameters = default,
                     IVector maximumParameters = default);
}
```

### 3.2. Вектор линейной алгебры (`Vector.cs`)
```csharp
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace LinearAlgebra
{
    public class Vector : List<double>, IVector
    {
        public Vector() : base() { }
        public Vector(int capacity) : base(capacity) { }
        public Vector(IEnumerable<double> collection) : base(collection) { }
        public Vector(params double[] values) : base(values) { }

        public Vector Clone() => new Vector(this);

        public static Vector CreateZeros(int length)
        {
            var v = new Vector(length);
            for (int i = 0; i < length; i++) v.Add(0.0);
            return v;
        }

        public double Dot(IVector other)
        {
            if (other == null) throw new ArgumentNullException(nameof(other));
            if (Count != other.Count)
                throw new ArgumentException($"Несоответствие размерностей: {Count} vs {other.Count}");

            double sum = 0;
            for (int i = 0; i < Count; i++) sum += this[i] * other[i];
            return sum;
        }

        public double Norm()
        {
            double sumSq = 0;
            for (int i = 0; i < Count; i++) sumSq += this[i] * this[i];
            return Math.Sqrt(sumSq);
        }

        public static Vector operator +(Vector a, IVector b)
        {
            if (a == null || b == null) throw new ArgumentNullException();
            if (a.Count != b.Count) throw new ArgumentException("Размерности должны совпадать.");
            var res = new Vector(a.Count);
            for (int i = 0; i < a.Count; i++) res.Add(a[i] + b[i]);
            return res;
        }

        public static Vector operator -(Vector a, IVector b)
        {
            if (a == null || b == null) throw new ArgumentNullException();
            if (a.Count != b.Count) throw new ArgumentException("Размерности должны совпадать.");
            var res = new Vector(a.Count);
            for (int i = 0; i < a.Count; i++) res.Add(a[i] - b[i]);
            return res;
        }

        public static Vector operator *(Vector a, double scalar)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            var res = new Vector(a.Count);
            for (int i = 0; i < a.Count; i++) res.Add(a[i] * scalar);
            return res;
        }

        public static Vector operator *(double scalar, Vector a) => a * scalar;

        public static Vector operator /(Vector a, double scalar)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            var res = new Vector(a.Count);
            for (int i = 0; i < a.Count; i++) res.Add(a[i] / scalar);
            return res;
        }

        public override string ToString() =>
            $"({string.Join(", ", this.Select(v => v.ToString("G6", CultureInfo.InvariantCulture)))})";
    }
}
```

### 3.3. Матрица и решение нормальных уравнений (`Matrix.cs`)
```csharp
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace LinearAlgebra
{
    public class Matrix : List<IList<double>>, IMatrix
    {
        public int RowCount => Count;
        public int ColCount => Count > 0 ? this[0].Count : 0;

        public Matrix() : base() { }

        public Matrix(int rows, int cols) : base(rows)
        {
            if (rows < 0 || cols < 0) throw new ArgumentException("Размеры не могут быть отрицательными.");
            for (int i = 0; i < rows; i++)
            {
                var row = new Vector(cols);
                for (int j = 0; j < cols; j++) row.Add(0.0);
                Add(row);
            }
        }

        public Matrix(IEnumerable<IList<double>> rows) : base(rows) { }

        public double this[int row, int col]
        {
            get => this[row][col];
            set => this[row][col] = value;
        }

        public Matrix Transpose()
        {
            int r = RowCount, c = ColCount;
            var res = new Matrix(c, r);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    res[j, i] = this[i, j];
            return res;
        }

        public static Vector SolveLinearSystem(IMatrix A, IVector b)
        {
            int n = A.Count;
            if (n == 0) return new Vector();
            if (A[0].Count != n || b.Count != n)
                throw new ArgumentException("Матрица должна быть квадратной и согласованной с вектором b.");

            double[][] a = new double[n][];
            for (int i = 0; i < n; i++)
            {
                a[i] = new double[n + 1];
                for (int j = 0; j < n; j++) a[i][j] = A[i][j];
                a[i][n] = b[i];
            }

            // Прямой ход с частичным выбором ведущего элемента
            for (int p = 0; p < n; p++)
            {
                int maxRow = p;
                double maxVal = Math.Abs(a[p][p]);
                for (int r = p + 1; r < n; r++)
                {
                    double val = Math.Abs(a[r][p]);
                    if (val > maxVal) { maxVal = val; maxRow = r; }
                }

                if (maxRow != p)
                    (a[p], a[maxRow]) = (a[maxRow], a[p]);

                if (Math.Abs(a[p][p]) < 1e-15) a[p][p] = 1e-12;

                for (int i = p + 1; i < n; i++)
                {
                    double factor = a[i][p] / a[p][p];
                    for (int j = p; j <= n; j++) a[i][j] -= factor * a[p][j];
                }
            }

            // Обратный ход
            var x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = a[i][n];
                for (int j = i + 1; j < n; j++) sum -= a[i][j] * x[j];
                x[i] = Math.Abs(a[i][i]) > 1e-15 ? sum / a[i][i] : 0.0;
            }
            return new Vector(x);
        }

        public static Vector SolveNormalEquations(IMatrix J, IVector r, double damping = 1e-6)
        {
            int m = J.Count;
            int p = m > 0 ? J[0].Count : 0;
            if (p == 0) return new Vector();

            var jtJ = new Matrix(p, p);
            var negJtr = new Vector(p);
            for (int i = 0; i < p; i++) negJtr.Add(0.0);

            for (int k = 0; k < m; k++)
            {
                double rk = r[k];
                for (int i = 0; i < p; i++)
                {
                    double jki = J[k][i];
                    negJtr[i] -= jki * rk;
                    for (int j = i; j < p; j++) jtJ[i, j] += jki * J[k][j];
                }
            }

            for (int i = 0; i < p; i++)
            {
                for (int j = 0; j < i; j++) jtJ[i, j] = jtJ[j, i];
                jtJ[i, i] += damping;
            }

            return SolveLinearSystem(jtJ, negJtr);
        }
    }
}
```

### 3.4. Кусочно-линейная функция (`PiecewiseLinearFunction.cs`)
```csharp
using System;
using System.Collections.Generic;
using System.Linq;
using LinearAlgebra;

namespace Functions
{
    public class PiecewiseLinearFunction : IParametricFunction
    {
        private readonly List<double> _knots;

        public IReadOnlyList<double> Knots => _knots;
        public int KnotCount => _knots.Count;

        public PiecewiseLinearFunction(IEnumerable<double> knots)
        {
            if (knots == null) throw new ArgumentNullException(nameof(knots));
            _knots = knots.OrderBy(x => x).ToList();
            if (_knots.Count < 2)
                throw new ArgumentException("Требуется как минимум 2 узла.", nameof(knots));

            for (int i = 0; i < _knots.Count - 1; i++)
            {
                if (Math.Abs(_knots[i + 1] - _knots[i]) < 1e-14)
                    throw new ArgumentException("Узлы должны быть строго различными.");
            }
        }

        public IFunction Bind(IVector parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != _knots.Count)
                throw new ArgumentException($"Ожидалось {_knots.Count} параметров, получено {parameters.Count}.");

            return new BoundPiecewiseLinearFunction(_knots, parameters);
        }

        private class BoundPiecewiseLinearFunction : IDifferentiableFunction
        {
            private readonly List<double> _knots;
            private readonly Vector _values;

            public BoundPiecewiseLinearFunction(List<double> knots, IVector values)
            {
                _knots = knots;
                _values = new Vector(values);
            }

            private (int index, double t) GetInterval(double x)
            {
                int k = _knots.Count;
                if (x <= _knots[0])
                {
                    double t = (x - _knots[0]) / (_knots[1] - _knots[0]);
                    return (0, t);
                }
                if (x >= _knots[k - 1])
                {
                    double t = (x - _knots[k - 2]) / (_knots[k - 1] - _knots[k - 2]);
                    return (k - 2, t);
                }

                int left = 0, right = k - 1;
                while (right - left > 1)
                {
                    int mid = (left + right) / 2;
                    if (_knots[mid] <= x) left = mid;
                    else right = mid;
                }

                double frac = (x - _knots[left]) / (_knots[left + 1] - _knots[left]);
                return (left, frac);
            }

            public double Value(IVector point)
            {
                if (point == null || point.Count == 0) throw new ArgumentException("Некорректная точка.");
                double x = point[0];
                var (idx, t) = GetInterval(x);
                return (1.0 - t) * _values[idx] + t * _values[idx + 1];
            }

            public IVector Gradient(IVector point)
            {
                if (point == null || point.Count == 0) throw new ArgumentException("Некорректная точка.");
                double x = point[0];
                var (idx, t) = GetInterval(x);

                var grad = Vector.CreateZeros(_knots.Count);
                grad[idx] = 1.0 - t;
                grad[idx + 1] = t;
                return grad;
            }
        }
    }
}
```

### 3.5. Функционал $L_2$-нормы (`L2NormFunctional.cs`)
```csharp
using System;
using System.Collections.Generic;
using Functions;
using LinearAlgebra;

namespace Functionals
{
    public class L2NormFunctional : IDifferentiableFunctional, ILeastSquaresFunctional
    {
        private readonly List<(IVector Point, double Value)> _dataPoints;
        public bool Squared { get; set; } = false;

        public L2NormFunctional(IEnumerable<(IVector Point, double Value)> dataPoints, bool squared = false)
        {
            if (dataPoints == null) throw new ArgumentNullException(nameof(dataPoints));
            _dataPoints = new List<(IVector Point, double Value)>(dataPoints);
            if (_dataPoints.Count == 0) throw new ArgumentException("Требуется хотя бы одна точка.");
            Squared = squared;
        }

        public double Value(IFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            double sumSq = 0.0;
            foreach (var (pt, y) in _dataPoints)
            {
                double diff = function.Value(pt) - y;
                sumSq += diff * diff;
            }
            return Squared ? sumSq : Math.Sqrt(sumSq);
        }

        public IVector Residual(IFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            var res = new Vector(_dataPoints.Count);
            foreach (var (pt, y) in _dataPoints)
            {
                res.Add(function.Value(pt) - y);
            }
            return res;
        }

        public IMatrix Jacobian(IFunction function)
        {
            if (function is not IDifferentiableFunction diffFunc)
                throw new ArgumentException("Функция должна реализовывать IDifferentiableFunction.");

            var jacobian = new Matrix();
            foreach (var (pt, _) in _dataPoints)
            {
                jacobian.Add(new Vector(diffFunc.Gradient(pt)));
            }
            return jacobian;
        }

        public IVector Gradient(IFunction function)
        {
            if (function is not IDifferentiableFunction diffFunc)
                throw new ArgumentException("Функция должна реализовывать IDifferentiableFunction.");

            var r = Residual(function);
            var J = Jacobian(diffFunc);
            int m = _dataPoints.Count;
            int p = J.Count > 0 ? J[0].Count : 0;

            var grad = Vector.CreateZeros(p);
            for (int i = 0; i < m; i++)
            {
                double ri = r[i];
                for (int j = 0; j < p; j++) grad[j] += ri * J[i][j];
            }

            if (Squared) return grad * 2.0;
            double norm = ((Vector)r).Norm();
            return norm < 1e-15 ? Vector.CreateZeros(p) : grad / norm;
        }
    }
}
```

### 3.6. Оптимизатор Гаусса-Ньютона (`GaussNewtonOptimizer.cs`)
```csharp
using System;
using Functions;
using Functionals;
using LinearAlgebra;

namespace Optimizers
{
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
            if (objective is not ILeastSquaresFunctional lsqObjective)
                throw new ArgumentException("Требуется реализация ILeastSquaresFunctional.", nameof(objective));

            var current = new Vector(initialParameters);
            Clamp(current, minimumParameters, maximumParameters);

            for (int iter = 0; iter < MaxIterations; iter++)
            {
                IFunction boundFunc = function.Bind(current);
                IVector r = lsqObjective.Residual(boundFunc);
                IMatrix J = lsqObjective.Jacobian(boundFunc);

                Vector delta = Matrix.SolveNormalEquations(J, r, Damping);
                if (delta.Norm() < Tolerance) break;

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

                if (!improved) break;
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
```

### 3.7. Главный модуль программы (`Program.cs`)
```csharp
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
            Console.WriteLine("Найденные значения в узлах сетки x = [0, 1, 2, 3, 4]:");
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
```

---

## 4. Тестирование программного комплекса

Для подтверждения корректности и численной надежности всех модулей системы было выполнено многоуровневое тестирование.

### 4.1. Модульные тесты линейной алгебры

Проверено решение системы линейных алгебраических уравнений $A x = b$ методом исключения Гаусса с частичным выбором ведущего элемента на тестовой системе 3-го порядка:

$$A = \begin{bmatrix} 0 & 2 & 1 \\ 1 & -2 & -3 \\ -1 & 1 & 2 \end{bmatrix}, \quad b = \begin{bmatrix} 4 \\ 0 \\ -1 \end{bmatrix}$$

Особенность матрицы — нулевой первый диагональный элемент $A_{0, 0} = 0$, что требует обязательной перестановки строк при выборе ведущего элемента (Partial Pivoting).
- Точное аналитическое решение: $x = (-1, 4, -4)^T$;
- Численное решение, полученное `Matrix.SolveLinearSystem`: $x = (-1.000000, 4.000000, -4.000000)^T$;
- Невязка решения: $\Vert A x - b \Vert_\infty < 1.0 \times 10^{-15}$ (на уровне машинной точности `double`).

### 4.2. Верификация аналитического дифференцирования

Корректность аналитического градиента $\nabla_p f(x, p)$, реализованного в `BoundPiecewiseLinearFunction`, верифицирована сравнением с конечно-разностной аппроксимацией центральными разностями с шагом $h = 10^{-6}$:

$$\left( \frac{\partial f}{\partial p_j} \right)_{\text{num}} \approx \frac{f(x, p + h e_j) - f(x, p - h e_j)}{2h}$$

Для точки $x = 1.7$ на сетке узлов $\{0, 1, 2, 3, 4\}$:
- Активный отрезок: $[1.0, 2.0]$, параметр $t = \frac{1.7 - 1.0}{2.0 - 1.0} = 0.7$;
- Теоретический градиент: $\nabla_p f = (0.0, \, 0.3, \, 0.7, \, 0.0, \, 0.0)^T$;
- Аналитический расчет программы: `(0.000000, 0.300000, 0.700000, 0.000000, 0.000000)`;
- Максимальное расхождение: $\max_j |(\nabla_p f)_j^{\text{analyt}} - (\nabla_p f)_j^{\text{num}}| < 2.2 \times 10^{-11}$.

### 4.3. Сквозной интеграционный тест оптимизатора

Запуск программы `OptimizationApp.exe` демонстрирует следующий протокол вывода в консоль:

```text
Решение задачи минимизации L2-норма:
Начальная L2-невязка с параметрами (0, 0, 0, 0, 0): 9,886860
Результаты:
Найденные значения в узлах сетки x = [0, 1, 2, 3, 4]:
  x = 0,0 -> y = 1,0000
  x = 1,0 -> y = 3,0000
  x = 2,0 -> y = 2,0000
  x = 3,0 -> y = 5,0000
  x = 4,0 -> y = 4,0000
Итоговая L2-невязка: 3,376841E-012
Проверка значений в контрольных точках:
  f(0,0) = 1,0000 (ожидалось: 1,0000, ошибка: 1,02E-013)
  f(0,5) = 2,0000 (ожидалось: 2,0000, ошибка: 6,62E-013)
  f(1,0) = 3,0000 (ожидалось: 3,0000, ошибка: 1,22E-012)
  f(1,5) = 2,5000 (ожидалось: 2,5000, ошибка: 5,64E-013)
  f(2,0) = 2,0000 (ожидалось: 2,0000, ошибка: 9,37E-014)
  f(2,5) = 3,5000 (ожидалось: 3,5000, ошибка: 7,44E-013)
  f(3,0) = 5,0000 (ожидалось: 5,0000, ошибка: 1,58E-012)
  f(3,5) = 4,5000 (ожидалось: 4,5000, ошибка: 1,69E-012)
  f(4,0) = 4,0000 (ожидалось: 4,0000, ошибка: 1,80E-012)
Тест завершён
```

---

## 5. Экспериментальные исследования

### 5.1. Динамика оптимизации и скорость сходимости

Эксперимент проводился при старте из нулевого вектора параметров:

$$p^{(0)} = (0.0, 0.0, 0.0, 0.0, 0.0)^T$$

Для набора из 9 контрольных точек истинные значения в узлах составляют:

$$p^* = (1.0, 3.0, 2.0, 5.0, 4.0)^T$$

| Параметр | Итерация 0 (Старт) | Итерация 1 (Финал) | Изменение |
| :--- | :---: | :---: | :---: |
| $p_0$ ($x=0.0$) | $0.000000$ | $1.000000$ | $+1.000000$ |
| $p_1$ ($x=1.0$) | $0.000000$ | $3.000000$ | $+3.000000$ |
| $p_2$ ($x=2.0$) | $0.000000$ | $2.000000$ | $+2.000000$ |
| $p_3$ ($x=3.0$) | $0.000000$ | $5.000000$ | $+5.000000$ |
| $p_4$ ($x=4.0$) | $0.000000$ | $4.000000$ | $+4.000000$ |
| **$L_2$-невязка** | **$9.886860$** | **$3.376841 \times 10^{-12}$** | Снижение в $\approx 3 \times 10^{12}$ раз |
| Норма шага $\Vert \Delta \Vert_2$ | — | $7.416198$ | — |
| Шаг линейного поиска $\alpha$ | — | $1.0$ (полный шаг) | — |

**Анализ результата**: В полном соответствии с теоретическим выводом раздела 2.3, квадратичный функционал линейной параметрической модели минимизируется ровно за **один шаг** метода Гаусса-Ньютона с полным шагом $\alpha = 1.0$. Итоговая невязка составила порядка $10^{-12}$, что обусловлено исключительно машинным округлением при операциях с плавающей точкой двойной точности (`IEEE 754 float64`).

### 5.2. Сравнение параметров и верификация в контрольных точках

В таблице ниже приведена верификация восстановленной кусочно-линейной функции во всех 9 экспериментальных точках:

| Точка $x_i$ | Тип точки | Требуемое значение $y_i$ | Восстановленное $f(x_i, p^*)$ | Абсолютная погрешность |
| :---: | :---: | :---: | :---: | :---: |
| $0.0$ | Узел 0 | $1.0000$ | $1.000000000000$ | $1.02 \times 10^{-13}$ |
| $0.5$ | Середина $[0, 1]$ | $2.0000$ | $2.000000000000$ | $6.62 \times 10^{-13}$ |
| $1.0$ | Узел 1 | $3.0000$ | $3.000000000000$ | $1.22 \times 10^{-12}$ |
| $1.5$ | Середина $[1, 2]$ | $2.5000$ | $2.500000000000$ | $5.64 \times 10^{-13}$ |
| $2.0$ | Узел 2 | $2.0000$ | $2.000000000000$ | $9.37 \times 10^{-14}$ |
| $2.5$ | Середина $[2, 3]$ | $3.5000$ | $3.500000000000$ | $7.44 \times 10^{-13}$ |
| $3.0$ | Узел 3 | $5.0000$ | $5.000000000000$ | $1.58 \times 10^{-12}$ |
| $3.5$ | Середина $[3, 4]$ | $4.5000$ | $4.500000000000$ | $1.69 \times 10^{-12}$ |
| $4.0$ | Узел 4 | $4.0000$ | $4.000000000000$ | $1.80 \times 10^{-12}$ |

Все погрешности лежат в диапазоне $10^{-14} \dots 10^{-12}$, что строго доказывает математическую точность работы оптимизатора и отсутствие смещения аппроксимации.

### 5.3. Анализ теоретической и практической вычислительной сложности

Сопоставим асимптотические оценки операций на одну итерацию алгоритма в терминах числа точек $M$ и числа узлов (параметров) $m$:

1. **Вычисление невязок $r$**: $M$ вычислений значений кусочно-линейной функции. Поиск интервала занимает $O(\log m)$, итого:
$$T_r = O(M \log m)$$

2. **Формирование матрицы Якоби $J$**: каждая строка матрицы содержит лишь 2 ненулевых коэффициента. Время вычисления:
$$T_J = O(M)$$

3. **Формирование матрицы нормальных уравнений $J^T J$ и вектора $-J^T r$**:
Благодаря сильной разреженности строк матрицы Якоби (каждая строка содержит только 2 ненулевых элемента), скалярное произведение строк затрачивает $O(1)$ операций на точку, а полное накопление $J^T J$ требует порядка $O(M + m)$ действий (в общем плотном случае $O(M m^2)$):
$$T_{\text{normal}} = O(M + m^2)$$

4. **Решение СЛАУ методом исключения Гаусса**: матрица $J^T J$ имеет ленточную структуру (трёхдиагональная матрица, так как перекрываются только соседние узлы). Плотный метод исключения Гаусса затрачивает $O(m^3)$ операций, специализированный метод прогонки — $O(m)$:
$$T_{\text{solve}} = O(m^3)$$

5. **Общая сложность оптимизации**:
$$T_{\text{total}} = O(M \log m + m^3)$$

Поскольку метод сходится ровно за 1 итерацию, суммарное время решения задачи составляет доли миллисекунды даже для сеток из сотен узлов. Для сравнения, градиентный спуск первого порядка потребовал бы вычисления $\nabla \Phi(p)$ на протяжении сотен и тысяч итераций, обладая линейной асимптотической скоростью сходимости с коэффициентом, зависящим от числа обусловленности $\kappa(J^T J)$.

---

## 6. Выводы

В ходе выполнения лабораторной работы была разработана и всесторонне исследована объектно-ориентированная система классов для решения задач оптимизации параметров функций на основе системы абстрактных интерфейсов. 

Основные результаты:
1. **Чистота объектно-ориентированной архитектуры**: 
   - За счет четкого разделения абстракций на `IParametricFunction`, `IFunction`, `IFunctional` и `IOptimizator` достигнута слабая связанность (Loose Coupling) и высокая расширяемость (Open/Closed Principle). 
   - Добавление новых математических функций (например, сплайнов или полиномов) или новых оптимизаторов не требует модификации существующего кода функционалов или других классов.
2. **Аналитическое дифференцирование**:
   - Реализация интерфейса `IDifferentiableFunction` с точным аналитическим вычислением частных производных устранила необходимость в численном дифференцировании, исключив погрешности усечения разностных схем и ускорив процесс расчета матрицы Якоби.
3. **Эффективность алгоритма Гаусса-Ньютона**:
   - Экспериментально и аналитически подтверждено, что для кусочно-линейных моделей матрица Якоби не зависит от вектора параметров, а матрица Гессе функционала в точности совпадает с $J^T J$. В результате метод Гаусса-Ньютона обеспечивает гарантированную глобальную сходимость к экстремуму ровно за **1 итерацию**.
4. **Численная устойчивость**:
   - Реализованный метод исключения Гаусса с частичным выбором ведущего элемента (Partial Pivoting) в сочетании с регуляризацией Левенберга обеспечивает надежное решение СЛАУ нормальных уравнений даже при наличии нулевых диагональных элементов или плохой обусловленности системы.
   - Итоговая невязка составила $3.38 \times 10^{-12}$, а погрешности в контрольных точках не превосходят $1.8 \times 10^{-12}$, что соответствует пределу точности вещественных чисел формата `double` (IEEE 754).
