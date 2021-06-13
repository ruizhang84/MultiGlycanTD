using MultiGlycanTDLibrary.engine.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class FragmentMassUnitTest
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

            GlycanIonsBuilder.Build.Permethylated = true;
            GlycanIonsBuilder.Build.Reduced = true;
            GlycanIonsBuilder.Build.Types = new List<FragmentTypes>()
            {
                FragmentTypes.B, FragmentTypes.C, FragmentTypes.Y, FragmentTypes.Z,
                FragmentTypes.BY, FragmentTypes.BZ, FragmentTypes.CY, FragmentTypes.YY,
                FragmentTypes.YZ, FragmentTypes.ZZ
            };

            string path = @"C:\Users\iruiz\Downloads\MSMS\builds.csv";
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("glycan_id,name,mass");
                    foreach(var pair in map)
                    {
                        var id = pair.Key;
                        var glycan = pair.Value;
                        if (glycan.IsValid())
                        {
                            List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                                .OrderBy(m => m).Select(m => Math.Round(m, 4)).ToList();
                            string output = glycan.ID() + "," + glycan.Name() + ","
                                + string.Join(" ", massList.Select(m => m.ToString()));
                            writer.WriteLine(output);
                        }
                    }

                    writer.Flush();
                }
            }
        }

    }
}


