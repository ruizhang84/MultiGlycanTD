using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SerializationJasonTest
    {
        [Test]
        public void JasonTest()
        {

            GlycanBuilder glycanBuilder =
                new GlycanBuilder(6, 6, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            var map = glycanBuilder.GlycanMaps();
            Console.WriteLine(map.Count);

            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();

            Object obj = new Object();

            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();
            CompdJson compdJson = new CompdJson()
            {
                DistrMap = distr_map,
                MassMap = mass_map
            };

            Dictionary<string, List<double>> fragments = new Dictionary<string, List<double>>();
            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                        .OrderBy(m => m).Select(m => Math.Round(m, 4)).ToList();
                    lock (obj)
                    {
                        fragments[id] = massList;
                    }
                }
            });
            var composition_map = glycanBuilder.GlycanCompositionMaps();
            Dictionary<string, List<string>> id_map = new Dictionary<string, List<string>>();
            foreach (var pair in composition_map)
            {
                id_map[pair.Key] = pair.Value.Select(p => p.ID()).ToList();
            }

            GlycanJson glycanJson = new GlycanJson()
            {
                Compound = compdJson,
                IDMap = id_map,
                Fragments = fragments
            };

            watch.Stop();
            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

            string fileName = @"C:\Users\iruiz\Downloads\massList.json";
            string jsonString = JsonSerializer.Serialize(glycanJson);
            File.WriteAllText(fileName, jsonString);

            string jsonStringRead = File.ReadAllText(fileName);
            GlycanJson glycanJsonRead = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            Assert.AreEqual(map.Count, glycanJsonRead.Fragments.Count);


        }
    }
}
