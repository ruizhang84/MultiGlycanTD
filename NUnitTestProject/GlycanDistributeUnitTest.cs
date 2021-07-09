using MultiGlycanTDLibrary.engine.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class GlycanDistributeUnitTest
    {
        [Test]
        public void DistributeTest()
        {
            GlycanBuilder glycanBuilder =
               new GlycanBuilder(7, 7, 5, 6, 0,
               true, false, false,
               10, true, true);
            glycanBuilder.Build();

            // distribution maps
            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();

            foreach (double mass in distr_map["GlcNAc-4 Man-3 Gal-2 Fuc-1 NeuAc-1 "])
            {
                Console.WriteLine(mass);
            }

            Console.WriteLine();

            var sorted = distr_map["GlcNAc-4 Man-3 Gal-2 Fuc-1 NeuAc-1 "]
                    .Select((x, i) => new KeyValuePair<double, int>(x, i))
                    .OrderByDescending(x => x.Key)
                    .Take(3).Where(x => x.Key > 0.05)
                    .ToList();
            List<int> idx = sorted.Select(x => x.Value).ToList();



            foreach (double mass in mass_map["GlcNAc-4 Man-3 Gal-2 Fuc-1 NeuAc-1 "])
            {
                double mz = MultiGlycanTDLibrary.util.mass.Spectrum.To.ComputeMZ(mass, 1.0078, 3);
                Console.WriteLine(mz);
            }

        }
    }
}
