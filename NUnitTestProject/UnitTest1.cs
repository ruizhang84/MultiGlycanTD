using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace NUnitTestProject
{
    public class Tests
    {
        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void Test1()
        {
            GlycanBuilder glycanBuilder =
                new GlycanBuilder(6, 6, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            var map = glycanBuilder.GlycanMaps();
            Console.WriteLine(map.Count);

            string path = @"C:\Users\Rui Zhang\Downloads\fragments.csv";
            MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("glycan_id,fragments");
                    foreach (var id in map.Keys)
                    {
                        //2 3 0 0 3 1 0 0 3 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
                        var glycan = map[id];
                        if (!glycan.IsValid())
                            continue;
                        List<string> massList = GlycanIonsBuilder.Build.Yions(glycan)
                                                .OrderBy(m => m).Select(m => Math.Round(m, 4).ToString()).ToList();
                        string output = id + "," + string.Join(" ", massList);
                        writer.WriteLine(output);
                    }
                    writer.Flush();
                }
            }
            Assert.Pass();
        }
    }
}