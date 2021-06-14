using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SerializationJasonTestV2
    {
        public List<SortedDictionary<Monosaccharide, int>> ReadFilter()
        {
            string path = @"C:\Users\iruiz\Downloads\MSMS\309 N-glycan search space for HGI.txt";
            Regex HexNAc = new Regex("HexNAc(\\d+)", RegexOptions.Compiled);
            Regex Hex = new Regex("Hex(\\d+)", RegexOptions.Compiled);
            Regex Fuc = new Regex("Fuc(\\d+)", RegexOptions.Compiled);
            Regex NeuAc = new Regex("NeuAc(\\d+)", RegexOptions.Compiled);
            Regex NeuGc = new Regex("NeuGc(\\d+)", RegexOptions.Compiled);

            List<SortedDictionary<Monosaccharide, int>> Filtered = 
                new List<SortedDictionary<Monosaccharide, int>>();
            using (FileStream fileStream = new FileStream(path, FileMode.Open, FileAccess.Read))
            {
                using (StreamReader sr = new StreamReader(fileStream))
                {
                    string line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        if (line.StartsWith("%"))
                        {
                            continue;
                        }
                        SortedDictionary<Monosaccharide, int> temp 
                            = new SortedDictionary<Monosaccharide, int>();
                        if (HexNAc.IsMatch(line))
                        {
                            MatchCollection matches = HexNAc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.HexNAc] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (Hex.IsMatch(line))
                        {
                            MatchCollection matches = HexNAc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.Hex] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (Fuc.IsMatch(line))
                        {
                            MatchCollection matches = HexNAc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.Fuc] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (NeuAc.IsMatch(line))
                        {
                            MatchCollection matches = HexNAc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.NeuAc] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (NeuGc.IsMatch(line))
                        {
                            MatchCollection matches = HexNAc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.NeuGc] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        Filtered.Add(temp);
                    }
                }
            }
            return Filtered;
        }

        [Test]
        public void JasonTestV2()
        {

            GlycanBuilder glycanBuilder =
                new GlycanBuilder(6, 6, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            //Console.WriteLine(map.Count);

            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();
            CompdJson compdJson = new CompdJson()
            {
                DistrMap = distr_map,
                MassMap = mass_map
            };

            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();
            object obj = new object();
            var map = glycanBuilder.GlycanMaps();
            Dictionary<double, List<string>> fragments = new Dictionary<double, List<string>>();
            List<Tuple<string, List<double>>> fragmentsContainer =
                new List<Tuple<string, List<double>>>();

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
                        fragmentsContainer.Add(Tuple.Create(id, massList));
                    }
                }
            });

            foreach (Tuple<string, List<double>> item in fragmentsContainer)
            {
                string id = item.Item1;
                foreach (double mass in item.Item2)
                {
                    if (!fragments.ContainsKey(mass))
                        fragments[mass] = new List<string>();
                    fragments[mass].Add(id);
                }
            }
            watch.Stop();
            Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds} ms");

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

            string fileName = @"C:\Users\iruiz\Downloads\massList.json";
            string jsonString = JsonSerializer.Serialize(glycanJson);
            File.WriteAllText(fileName, jsonString);

            string jsonStringRead = File.ReadAllText(fileName);
            GlycanJson glycanJsonRead = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            //Assert.AreEqual(map.Count, glycanJsonRead.Fragments.Count);


        }
    }
}
