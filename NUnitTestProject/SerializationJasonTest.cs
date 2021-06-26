using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SerializationJasonTest
    {
        public static byte[] SerializeAndCompress(List<double> massList)
        {
            using (MemoryStream ms = new MemoryStream())
            using (GZipStream zs = new GZipStream(ms, CompressionMode.Compress, true))
            {
                BinaryFormatter bf = new BinaryFormatter();
                bf.Serialize(zs, massList);
                return ms.ToArray();
            }
        }

        public static T DecompressAndDeserialize<T>(byte[] data)
        {
            using (MemoryStream ms = new MemoryStream(data))
            using (GZipStream zs = new GZipStream(ms, CompressionMode.Decompress, true))
            {
                BinaryFormatter bf = new BinaryFormatter();
                return (T)bf.Deserialize(zs);
            }
        }

        [Test]
        public void JasonTest()
        {

            GlycanBuilder glycanBuilder =
                new GlycanBuilder(7, 7, 5, 4, 0, true, false, false);
            glycanBuilder.Build();


            //var distr_map = glycanBuilder.GlycanDistribMaps();
            //var mass_map = glycanBuilder.GlycanMassMaps();
            //CompdJson compdJson = new CompdJson()
            //{
            //    DistrMap = distr_map,
            //    MassMap = mass_map
            //};

            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();
            object obj = new object();
            var map = glycanBuilder.GlycanMaps();
            Dictionary<double, List<string>> fragments = new Dictionary<double, List<string>>();
            List<Tuple<string, byte[]>> fragmentsContainer =
                new List<Tuple<string, byte[]>>();

            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                        .Select(m => Math.Round(m, 4)).ToList();
                    lock (obj)
                    {
                        fragmentsContainer.Add(Tuple.Create(id, SerializeAndCompress(massList)));
                    }
                }
            });

            //foreach (Tuple<string, List<double>> item in fragmentsContainer)
            //{
            //    string id = item.Item1;
            //    foreach (double mass in item.Item2)
            //    {
            //        if (!fragments.ContainsKey(mass))
            //            fragments[mass] = new List<string>();
            //        fragments[mass].Add(id);
            //    }
            //}


            watch.Stop();
            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");
            //var composition_map = glycanBuilder.GlycanCompositionMaps();
            //Dictionary<string, List<string>> id_map = new Dictionary<string, List<string>>();
            //foreach (var pair in composition_map)
            //{
            //    id_map[pair.Key] = pair.Value.Select(p => p.ID()).ToList();
            //}

            //GlycanJson glycanJson = new GlycanJson()
            //{
            //    Compound = compdJson,
            //    IDMap = id_map,
            //    FragmentMap = fragments
            //};

            watch.Stop();
            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

            //string fileName = @"C:\Users\iruiz\Downloads\massList.json";
            //string jsonString = JsonSerializer.Serialize(glycanJson);
            //File.WriteAllText(fileName, jsonString);

            //string jsonStringRead = File.ReadAllText(fileName);
            //GlycanJson glycanJsonRead = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            //Assert.AreEqual(map.Count, glycanJsonRead.FragmentMap.Count);


        }
    }
}
