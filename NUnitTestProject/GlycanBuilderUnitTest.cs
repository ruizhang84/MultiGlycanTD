using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class GlycanBuilderTest
    {
        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void Test1()
        {
            GlycanBuilder glycanBuilder =
                new GlycanBuilder(12, 12, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            var map = glycanBuilder.GlycanMaps();
            Console.WriteLine(map.Count);

            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();
            string output = "";
            Object obj = new Object();
            //foreach (var id in map.Keys)
            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                //2 3 0 0 3 1 0 0 3 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<string> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                        .OrderBy(m => m).Select(m => Math.Round(m, 4).ToString()).ToList();
                    lock(obj)
                    {
                        output += id + "," + string.Join(" ", massList) + "\n";
                    }
                }
            });
            watch.Stop();
            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

            string path = @"C:\Users\Rui Zhang\Downloads\fragments.csv";
            MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("glycan_id,fragments");
                    writer.WriteLine(output);
                    writer.Flush();
                }
            }
            Assert.Pass();
        }
    }
}