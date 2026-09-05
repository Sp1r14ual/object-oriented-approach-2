using System;
using LinearAlgebra;

namespace Functions
{
    /// <summary>
    /// Полином n-й степени в одномерном пространстве (число параметров n+1).
    /// P_n(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_n * x^n
    /// По требованию задания: НЕ реализует IDifferentiableFunction (только IFunction).
    /// </summary>
    public class PolynomialFunction : IParametricFunction
    {
        private readonly int _degree;

        public int Degree => _degree;

        public PolynomialFunction(int degree = 1)
        {
            if (degree < 0)
                throw new ArgumentException("Degree must be non-negative.", nameof(degree));
            _degree = degree;
        }

        public IFunction Bind(IVector parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != _degree + 1)
            {
                throw new ArgumentException($"Expected {_degree + 1} parameters for degree {_degree} polynomial, but got {parameters.Count}.", nameof(parameters));
            }

            return new BoundPolynomialFunction(parameters);
        }

        /// <summary>
        /// Реализует ТОЛЬКО IFunction, намеренно НЕ реализует IDifferentiableFunction.
        /// </summary>
        private class BoundPolynomialFunction : IFunction
        {
            private readonly Vector _coefficients;

            public BoundPolynomialFunction(IVector parameters)
            {
                _coefficients = new Vector(parameters);
            }

            public double Value(IVector point)
            {
                if (point == null) throw new ArgumentNullException(nameof(point));
                if (point.Count == 0)
                    throw new ArgumentException("Point must contain at least 1 coordinate.");

                double x = point[0];

                // Вычисление полинома по схеме Горнера: a_0 + x * (a_1 + x * (a_2 + ...))
                double result = 0.0;
                for (int i = _coefficients.Count - 1; i >= 0; i--)
                {
                    result = result * x + _coefficients[i];
                }
                return result;
            }
        }
    }
}
