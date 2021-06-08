using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.Distribution
{
    public class Digamma
    {
        public static double Value(double x, int order = 14)
        {
            if (x > 6)
            {
                return Expansion(x, order);
            }
            int step = (int) Math.Ceiling(6 - x);
            double shift = 0;
            for (int i = 0; i < step; i++)
            {
                shift += 1 / x;
                x++;
            }
            return Expansion(x, order) - shift;
        }

        public static double Expansion(double x, int order)
        {
            double sum = 0;
            for(int i = 1; i <= order; i++)
            {
                sum += BernoulliNumbers.Compute(2 * i) / (2 * i * Math.Pow(x, 2 * i));
            }
            return Math.Log(x) - 0.5 / x - sum;
        }
    }
}
