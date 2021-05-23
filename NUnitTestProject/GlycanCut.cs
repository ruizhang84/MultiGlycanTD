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
    public class GlycanCut
    {
        [Test]
        public void Test()
        {
            GlycanBuilder glycanBuilder =
                new GlycanBuilder(12, 12, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            var map = glycanBuilder.GlycanMaps();
            Console.WriteLine(map.Count);

            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();

            string path = @"C:\Users\Rui Zhang\Downloads\cuts.csv";
            MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("glycan_id,detail,name,cut");
                    string id = "2 1 0 0 1 1 3 1 1 0 2 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0";
                    if (map.ContainsKey(id))
                    {
                        var glycan = map[id];
                        foreach (var g in glycan.Fragments())
                        {
                            int diff = GlycanFragmentBuilderHelper.CountYCut(g, glycan, 10);
                            int[] table = g.Table();
                            string output = g.ID() + "," 
                                + table[6].ToString() + "-"
                                + table[7].ToString() + "-"
                                + table[8].ToString() + "-"
                                + table[9].ToString() + "|-"
                                + table[10].ToString() + "-"
                                + table[11].ToString() + "-"
                                + table[12].ToString() + "-"
                                + table[13].ToString() + "|-"
                                + table[14].ToString() + "-"
                                + table[15].ToString() + "-"
                                + table[16].ToString() + "-"
                                + table[17].ToString() + ","
                                + g.Name() + "," + diff.ToString();
                            writer.WriteLine(output);
                        }
                    }
                    writer.Flush();
                }
            }
            watch.Stop();

            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

            Assert.Pass();
            
        }
    }
}
