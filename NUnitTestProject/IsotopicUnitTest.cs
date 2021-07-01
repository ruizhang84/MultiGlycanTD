using SpectrumProcess.brain;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class IsotopicUnitTest
    {
        [Test]
        public void Test()
        {
            Dictionary<ElementType, int> composition = new Dictionary<ElementType, int>()
                { { ElementType.C, 50}, { ElementType.H, 71}, { ElementType.N, 13}, { ElementType.O, 12} };
            //{ { new C(), 8}, { new H(), 15}, { new N(), 1}, { new O(), 6} };


            Compound compound = new Compound(composition);

            int order = 10;
            List<double> coeff = Brain.Run.Distribute(compound, order);
            List<double> massList = Brain.Run.CenterMass(compound, order);
            for(int i = 0; i < order; i++)
            {
                string output = "Coeff: " + coeff[i].ToString() + 
                    " Mass: " + massList[i].ToString();
                 Console.WriteLine(output);
            }

            Assert.Pass();
        }
    }
}
