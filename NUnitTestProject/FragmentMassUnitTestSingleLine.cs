using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class FragmentMassUnitTestSingleLine
    {
        [Test]
        public void FragmentsTest()
        {
            GlycanBuilder glycanBuilder =
               new GlycanBuilder(7, 7, 5, 6, 0,
               true, false, false,
               10, true, true);
            glycanBuilder.Build();

            // distribution maps
            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();

            // fragmentation maps
            object obj = new object();
            Dictionary<double, List<string>> fragments = new Dictionary<double, List<string>>();
            var map = glycanBuilder.GlycanMaps();

            GlycanIonsBuilder.Build.Permethylated = false;
            GlycanIonsBuilder.Build.Reduced = false;
            MultiGlycanClassLibrary.util.mass.Glycan.To.permethylation = false;
            MultiGlycanClassLibrary.util.mass.Glycan.To.reduced = false;

            GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.kABEE;
            Glycan.To.Derivatization = Glycan.kABEE;

            GlycanIonsBuilder.Build.Types = new List<FragmentType>()
            { FragmentType.YZ };
            //{
            //    FragmentType.B, FragmentType.C, FragmentType.Y, FragmentType.Z,
            //    FragmentType.BY, FragmentType.BZ, FragmentType.CY, FragmentType.YY,
            //    FragmentType.YZ, FragmentType.ZZ
            //};
            foreach(var pair in map)
            {
                var id = pair.Key;
                //if (!id.StartsWith("2 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 0"))
                //    continue;
                if (!id.StartsWith("2 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"))
                    continue;
                var glycan = pair.Value;
                //foreach(IGlycan g in glycan.FragmentMap())
                //{
                //    Console.WriteLine(g.ID());
                //}
                Console.WriteLine(Glycan.To.Compute(glycan));
                if (glycan.IsValid())
                {
                    List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                        .OrderBy(m => m).Select(m => Math.Round(m, 4)).ToList();
                    Console.WriteLine(string.Join("\n", massList.Select(m => m.ToString())));
                }
            }
        }

    }
}




