using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.LanczosApproximation
{
    public class ChebyshevPolynomials
    {
        public int[,] Coefficiet { get; set; }

        public ChebyshevPolynomials(int n = 7)
        {
            Coefficiet = new int[2 * n + 2, 2 * n + 2];
        }

        public int Get(int n, int m)
        { 
            // base case
            if (n == 1 && m == 1)
            {
                Coefficiet[1, 1] = 1;
            }
            else if (n == 2 && m == 2)
            {
                Coefficiet[2, 2] = 1;
            }
            else if (m == 1)
            {
                Coefficiet[n, m] = -Get(n - 2, 1);
            }
            else if(n == m)
            {
                Coefficiet[n, m] = 2 * Get(n - 1, n - 1);
            }
            else if(n > m)
            {
                Coefficiet[n, m] = 
                    2 * Get(n - 1, m - 1) - Get(n - 2, m);
            }
            return Coefficiet[n, m];
        }


    }
}
