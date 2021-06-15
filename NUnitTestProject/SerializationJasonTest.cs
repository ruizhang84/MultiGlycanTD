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
                new GlycanBuilder(10, 10, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();

            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();
            CompdJson compdJson = new CompdJson()
            {
                DistrMap = distr_map,
                MassMap = mass_map
            };

            object obj = new object();
            Dictionary<double, List<string>> fragments =
                new Dictionary<double, List<string>>();
            List<Dictionary<double, List<string>>> fragmentsContainer =
                new List<Dictionary<double, List<string>>>();
            var map = glycanBuilder.GlycanMaps();

            //Parallel.ForEach(map, new ParallelOptions { MaxDegreeOfParallelism = thread }, pair =>
            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                Dictionary<double, List<string>> temp = new Dictionary<double, List<string>>();
                List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                        .Select(m => Math.Round(m, 2)).ToList();

                foreach (double mass in massList)
                {
                    if (!temp.ContainsKey(mass))
                        temp[mass] = new List<string>();
                    temp[mass].Add(id);
                }

                lock (obj)
                {
                    fragmentsContainer.Add(temp);
                }
                
            });
            foreach (Dictionary<double, List<string>> item in fragmentsContainer)
            {
                foreach (double mass in item.Keys)
                {
                    if (!fragments.ContainsKey(mass))
                        fragments[mass] = new List<string>();
                    fragments[mass].AddRange(item[mass]);
                }
            }

            fragmentsContainer.Clear();
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

            //string fileName = @"C:\Users\iruiz\Downloads\massList.json";
            //string jsonString = JsonSerializer.Serialize(glycanJson);
            //File.WriteAllText(fileName, jsonString);

            //string jsonStringRead = File.ReadAllText(fileName);
            //GlycanJson glycanJsonRead = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            //Assert.AreEqual(map.Count, glycanJsonRead.Fragments.Count);


        }
    }
}
