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

        public Vector Clone()
        {
            return new Vector(this);
        }

        public static Vector CreateZeros(int length)
        {
            var v = new Vector(length);
            for (int i = 0; i < length; i++)
            {
                v.Add(0.0);
            }
            return v;
        }

        public double Dot(IVector other)
        {
            if (other == null) throw new ArgumentNullException(nameof(other));
            if (Count != other.Count)
                throw new ArgumentException($"Несоответствие размерностей: {Count} vs {other.Count}");

            double sum = 0;
            for (int i = 0; i < Count; i++)
            {
                sum += this[i] * other[i];
            }
            return sum;
        }

        public double Norm()
        {
            double sumSq = 0;
            for (int i = 0; i < Count; i++)
            {
                sumSq += this[i] * this[i];
            }
            return Math.Sqrt(sumSq);
        }

        public static Vector operator +(Vector a, IVector b)
        {
            if (a == null || b == null) throw new ArgumentNullException();
            if (a.Count != b.Count) throw new ArgumentException("Размерности векторов должны совпадать.");
            var res = new Vector(a.Count);
            for (int i = 0; i < a.Count; i++) res.Add(a[i] + b[i]);
            return res;
        }

        public static Vector operator -(Vector a, IVector b)
        {
            if (a == null || b == null) throw new ArgumentNullException();
            if (a.Count != b.Count) throw new ArgumentException("Размерности векторов должны совпадать.");
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

        public override string ToString()
        {
            return $"({string.Join(", ", this.Select(v => v.ToString("G6", CultureInfo.InvariantCulture)))})";
        }
    }
}
