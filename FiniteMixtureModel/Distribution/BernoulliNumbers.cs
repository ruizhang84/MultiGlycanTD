using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.Distribution
{
    // Akiyama–Tanigawa algorithm
    // cite: https://rosettacode.org/wiki/Bernoulli_numbers
    public class BernoulliNumbers
    {
        // returns nth Bernoulli number
        public static double Compute(int n)
        {
            BigInteger f;
            BigInteger[] nu = new BigInteger[n + 1],
                         de = new BigInteger[n + 1];
            for (int m = 0; m <= n; m++)
            {
                nu[m] = 1; de[m] = m + 1;
                for (int j = m; j > 0; j--)
                    if ((f = BigInteger.GreatestCommonDivisor(
                        nu[j - 1] = j * (de[j] * nu[j - 1] - de[j - 1] * nu[j]),
                        de[j - 1] *= de[j])) != BigInteger.One)
                    { nu[j - 1] /= f; de[j - 1] /= f; }
            }
            double numer = (double) nu[0];
            double denomin =  (double) de[0];
            return numer / denomin;
        }
    }
}
