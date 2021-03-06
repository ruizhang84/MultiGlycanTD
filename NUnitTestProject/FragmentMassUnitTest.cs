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

            GlycanIonsBuilder.Build.Permethylated = false;
            GlycanIonsBuilder.Build.Reduced = false;
            MultiGlycanClassLibrary.util.mass.Glycan.To.permethylation = false;
            MultiGlycanClassLibrary.util.mass.Glycan.To.reduced = false;

            GlycanIonsBuilder.Build.Types = new List<FragmentType>()
            { FragmentType.YYY };
            //2 1 0 0 1 1 2 0 2 0 2 0 2 0 1 0 1 0 1 0 1 0 0 0 0 0
            string path = @"C:\Users\iruiz\Downloads\MSMS\test_build.csv";
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("glycan_id,name,mass");
                    foreach(var pair in map)
                    {
                        var id = pair.Key;
                        var glycan = pair.Value;
                        if (id.StartsWith("2 1 1 0 1 1 2 0 2 0 2 0 2 0 1 0 1 0 1 0 1 0 0 0 0 0"))
                        {
                            List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                                .OrderBy(m => m).Select(m => Math.Round(m, 4)).Distinct().ToList();
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


