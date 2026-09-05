using System;
using LinearAlgebra;

namespace Functions
{
    /// <summary>
    /// Линейная функция в n-мерном пространстве (число параметров n+1).
    /// f(x) = w_0 * x_0 + w_1 * x_1 + ... + w_{n-1} * x_{n-1} + b
    /// Реализует IDifferentiableFunction.
    /// </summary>
    public class LinearFunction : IParametricFunction
    {
        private readonly int _dimension;

        public int Dimension => _dimension;

        public LinearFunction(int dimension = 1)
        {
            if (dimension <= 0)
                throw new ArgumentException("Dimension must be positive.", nameof(dimension));
            _dimension = dimension;
        }

        public IFunction Bind(IVector parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != _dimension + 1)
            {
                throw new ArgumentException($"Expected {_dimension + 1} parameters for a {_dimension}-dimensional linear function, but got {parameters.Count}.", nameof(parameters));
            }

            return new BoundLinearFunction(parameters, _dimension);
        }

        private class BoundLinearFunction : IDifferentiableFunction
        {
            private readonly Vector _parameters;
            private readonly int _dimension;

            public BoundLinearFunction(IVector parameters, int dimension)
            {
                _dimension = dimension;
                _parameters = new Vector(parameters);
            }

            public double Value(IVector point)
            {
                if (point == null) throw new ArgumentNullException(nameof(point));
                if (point.Count < _dimension)
                    throw new ArgumentException($"Point must have at least {_dimension} dimensions, got {point.Count}.");

                double result = _parameters[_dimension]; // bias term b
                for (int i = 0; i < _dimension; i++)
                {
                    result += _parameters[i] * point[i];
                }
                return result;
            }

            /// <summary>
            /// Градиент функции по параметрам исходной IParametricFunction.
            /// df/dw_i = x_i, df/db = 1
            /// </summary>
            public IVector Gradient(IVector point)
            {
                if (point == null) throw new ArgumentNullException(nameof(point));
                if (point.Count < _dimension)
                    throw new ArgumentException($"Point must have at least {_dimension} dimensions, got {point.Count}.");

                var grad = new Vector(_dimension + 1);
                for (int i = 0; i < _dimension; i++)
                {
                    grad.Add(point[i]);
                }
                grad.Add(1.0); // df/db = 1
                return grad;
            }
        }
    }
}
