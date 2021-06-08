using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.LanczosApproximation
{
    public class Godfrey
    {
        public int g { get; set; }
        public int n { get; set; }
        public ChebyshevPolynomials polynomials;
        
        public Godfrey(int g = 5, int n = 7)
        {
            this.g = g;
            this.n = n;
            polynomials = new ChebyshevPolynomials(n);
        }

        private int DoubleFactorial(int n)
        {
            if (n < 2)
                return 1;

            return n * DoubleFactorial(n - 2);
        }

        public double Coefficient(int k)
        {
            //                          (2 * a - 1)!! e ^ (a + g + 1 / 2)
            //Fg(a) = sqrt(2 / pi) *    ------------------------
            //                          2 ^ a(a + g + 1 / 2) ^ (a + 1 / 2)
            double sum = 0;
            for(int l = 0; l <= k; l++)
            {
                double coeff = polynomials.Get(2 * k + 1, 2 * l + 1);
                double fac = DoubleFactorial(2 * l - 1);
                double power = Math.Pow(2, l) * Math.Pow(l + g + 0.5, (l + 0.5));
                double exp = Math.Exp(l + g + 0.5);
                sum += coeff * fac * exp / power;
            }
            return Math.Sqrt(2 / Math.PI)  * sum;
        }

        public double SeriesA(double z)
        {
            // http://www.numericana.com/answer/info/godfrey.htm
            throw new NotImplementedException();
        }

        public double Gamma(double z)
        {
            if (z < 0.5)
                return Math.PI / (Math.Sin(Math.PI * z) * Gamma(1 - z));
            double t = z + g + 0.5;
            return Math.Sqrt(Math.PI * 2) * Math.Pow(t, z + 0.5) *
                Math.Exp(-t) * SeriesA(z);

        }

    }
}
