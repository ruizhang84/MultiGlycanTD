﻿using MultiGlycanTDLibrary.engine.glycan;
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

            List<JsonEntry> entries = new List<JsonEntry>();
            Object obj = new Object();
            int count = 0;
            var composition_map = glycanBuilder.GlycanCompositionMaps();
            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();

            Dictionary<string, List<string>> id_map = new Dictionary<string, List<string>>();
            foreach (var pair in composition_map)
            {
                id_map[pair.Key] = pair.Value.Select(p => p.ID()).ToList();
            }

            CompdJson compdJson = new CompdJson()
            {
                IDMap = id_map,
                DistrMap = distr_map,
                MassMap = mass_map
            };

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
                        JsonEntry g = new JsonEntry()
                        {
                            ID = glycan.ID(),
                            Fragments = massList
                        };
                        entries.Add(g);
                        count++;
                    }
                }
            });
            GlycanJson glycanJson = new GlycanJson()
            {
                Compound = compdJson,
                Entries = entries
            };

            watch.Stop();
            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

            string fileName = @"C:\Users\iruiz\Downloads\massList.json";
            string jsonString = JsonSerializer.Serialize(glycanJson);
            File.WriteAllText(fileName, jsonString);

            string jsonStringRead = File.ReadAllText(fileName);
            GlycanJson glycanJsonRead = JsonSerializer.Deserialize<GlycanJson>(jsonString);
            Assert.AreEqual(count, glycanJsonRead.Entries.Count);


        }
    }
}
