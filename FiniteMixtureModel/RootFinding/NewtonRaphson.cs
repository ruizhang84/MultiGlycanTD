using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.RootFinding
{
    public class NewtonRaphson
    {
        public delegate double Func(double x);
        public delegate double DerivativeFunc(double x);
        public static double Find(Func f, DerivativeFunc df, double x, double eps=0.001)
        {
            double h = f(x) / df(x);
            while (Math.Abs(h) >= eps)
            {
                h = f(x) / df(x);
                x = x - h;
            }
            return x;
        }
    }
}
