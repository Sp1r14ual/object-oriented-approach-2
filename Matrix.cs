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
            if (rows < 0 || cols < 0)
                throw new ArgumentException("Размеры не могут быть отрицательными.");

            for (int i = 0; i < rows; i++)
            {
                var row = new Vector(cols);
                for (int j = 0; j < cols; j++)
                {
                    row.Add(0.0);
                }
                Add(row);
            }
        }

        public Matrix(IEnumerable<IList<double>> rows) : base(rows) { }

        public Matrix(double[,] array)
        {
            if (array == null) throw new ArgumentNullException(nameof(array));
            int rows = array.GetLength(0);
            int cols = array.GetLength(1);
            for (int i = 0; i < rows; i++)
            {
                var row = new Vector(cols);
                for (int j = 0; j < cols; j++)
                {
                    row.Add(array[i, j]);
                }
                Add(row);
            }
        }

        public double this[int row, int col]
        {
            get => this[row][col];
            set => this[row][col] = value;
        }

        public Matrix Transpose()
        {
            int r = RowCount;
            int c = ColCount;
            var res = new Matrix(c, r);
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    res[j, i] = this[i, j];
                }
            }
            return res;
        }

        public Vector Multiply(IVector vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            if (ColCount != vector.Count)
                throw new ArgumentException($"Количество столбцов матрицы ({ColCount}) должно совпадать с длиной вектора ({vector.Count}).");

            var res = new Vector(RowCount);
            for (int i = 0; i < RowCount; i++)
            {
                double sum = 0;
                for (int j = 0; j < ColCount; j++)
                {
                    sum += this[i, j] * vector[j];
                }
                res.Add(sum);
            }
            return res;
        }

        public static Matrix Multiply(IMatrix a, IMatrix b)
        {
            if (a == null || b == null) throw new ArgumentNullException();
            int aRows = a.Count;
            int aCols = aRows > 0 ? a[0].Count : 0;
            int bRows = b.Count;
            int bCols = bRows > 0 ? b[0].Count : 0;

            if (aCols != bRows)
                throw new ArgumentException($"Невозможно умножить матрицу {aRows}x{aCols} на {bRows}x{bCols}.");

            var res = new Matrix(aRows, bCols);
            for (int i = 0; i < aRows; i++)
            {
                for (int j = 0; j < bCols; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < aCols; k++)
                    {
                        sum += a[i][k] * b[k][j];
                    }
                    res[i, j] = sum;
                }
            }
            return res;
        }

        /// <summary>
        /// Решает линейную систему A * x = b методом Гаусса с частичным выбором.
        /// </summary>
        public static Vector SolveLinearSystem(IMatrix A, IVector b)
        {
            int n = A.Count;
            if (n == 0) return new Vector();
            if (A[0].Count != n)
                throw new ArgumentException("Матрица A должна быть квадратной.");
            if (b.Count != n)
                throw new ArgumentException("Размерность b должна совпадать с размеромностью матрицы A.");

            // Расширенная матрица
            double[][] a = new double[n][];
            for (int i = 0; i < n; i++)
            {
                a[i] = new double[n + 1];
                for (int j = 0; j < n; j++)
                {
                    a[i][j] = A[i][j];
                }
                a[i][n] = b[i];
            }

            // Прямой ход с частичным выбором
            for (int p = 0; p < n; p++)
            {
                int maxRow = p;
                double maxVal = Math.Abs(a[p][p]);
                for (int r = p + 1; r < n; r++)
                {
                    double val = Math.Abs(a[r][p]);
                    if (val > maxVal)
                    {
                        maxVal = val;
                        maxRow = r;
                    }
                }

                if (maxRow != p)
                {
                    (a[p], a[maxRow]) = (a[maxRow], a[p]);
                }

                if (Math.Abs(a[p][p]) < 1e-15)
                {
                    // Добавление небольшой регуляризации, если сингулярна
                    a[p][p] = 1e-12;
                }

                for (int i = p + 1; i < n; i++)
                {
                    double factor = a[i][p] / a[p][p];
                    for (int j = p; j <= n; j++)
                    {
                        a[i][j] -= factor * a[p][j];
                    }
                }
            }

            // Обратная подстановка
            var x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = a[i][n];
                for (int j = i + 1; j < n; j++)
                {
                    sum -= a[i][j] * x[j];
                }
                x[i] = Math.Abs(a[i][i]) > 1e-15 ? sum / a[i][i] : 0.0;
            }

            return new Vector(x);
        }

        /// <summary>
        /// Решает (J^T * J + lambda * I) * delta = -J^T * r
        /// </summary>
        public static Vector SolveNormalEquations(IMatrix J, IVector r, double damping = 1e-6)
        {
            int m = J.Count;
            int p = m > 0 ? J[0].Count : 0;
            if (p == 0) return new Vector();

            var jtJ = new Matrix(p, p);
            var negJtr = new Vector(p);

            for (int i = 0; i < p; i++)
            {
                negJtr.Add(0.0);
            }

            for (int k = 0; k < m; k++)
            {
                double rk = r[k];
                for (int i = 0; i < p; i++)
                {
                    double jki = J[k][i];
                    negJtr[i] -= jki * rk;

                    for (int j = i; j < p; j++)
                    {
                        jtJ[i, j] += jki * J[k][j];
                    }
                }
            }

            for (int i = 0; i < p; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    jtJ[i, j] = jtJ[j, i];
                }
                jtJ[i, i] += damping;
            }

            return SolveLinearSystem(jtJ, negJtr);
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            for (int i = 0; i < RowCount; i++)
            {
                sb.Append("[ ");
                for (int j = 0; j < ColCount; j++)
                {
                    sb.Append(this[i, j].ToString("F4", CultureInfo.InvariantCulture)).Append(" ");
                }
                sb.AppendLine("]");
            }
            return sb.ToString();
        }
    }
}
