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
    public class BuildTest
    {

        [Test]
        public void TestBuild()
        {
            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();

            GlycanBuilder glycanBuilder =
               new GlycanBuilder(12, 12, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            var map = glycanBuilder.GlycanMaps();
            Console.WriteLine(map.Count);

            //string path = @"C:\Users\Rui Zhang\Downloads\builds.csv";
            //MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            //using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            //{
            //    using (StreamWriter writer = new StreamWriter(ostrm))
            //    {
            //        writer.WriteLine("glycan_id,name");
            //        foreach(var name in map.Keys)
            //        {
            //            IGlycan g = map[name];
            //            string output = g.ID() + "," + g.Name();
            //            writer.WriteLine(output);
            //        }
            //        //string id = "2 1 0 0 1 1 3 1 0 0 3 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0";
            //        //if (map.ContainsKey(id))
            //        //{
            //        //    var glycan = map[id];
            //        //    foreach (var g in glycan.FragmentMap())
            //        //    {
            //        //        int diff = GlycanFragmentBuilderHelper.CountYCut(g, glycan, 10);
            //        //        //2 3 0 0 3 1 0 0 3 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
            //        //        string output = g.ID() + "," + g.Name() + "," + diff.ToString();
            //        //        writer.WriteLine(output);
            //        //    }
            //        //}
            //        writer.Flush();
            //    }
            //}
            watch.Stop();

            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

            Assert.Pass();
        }
    }
}
