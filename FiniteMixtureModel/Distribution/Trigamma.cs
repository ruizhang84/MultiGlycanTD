using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.Distribution
{
    public class Trigamma
    {
        public static double Value(double x, int order = 15)
        {
            if (x > 6)
            {
                return Expansion(x, order);
            }
            int step = (int)Math.Ceiling(6 - x);
            double shift = 0;
            for (int i = 0; i < step; i++)
            {
                shift += 1 / (x * x);
                x++;
            }
            return Expansion(x, order) + shift;
        }

        public static double Expansion(double x, int order)
        {
            double sum = 0;
            for (int i = 0; i <= order; i++)
            {
                sum += BernoulliNumbers.Compute(i) / (Math.Pow(x, i + 1));
            }
            return sum;
        }
    }
}
