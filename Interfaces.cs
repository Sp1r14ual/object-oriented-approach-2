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
