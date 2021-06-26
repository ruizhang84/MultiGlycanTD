using MultiGlycanClassLibrary.util.mass;
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
    public class SerializationJasonTestV3
    {
        public List<SortedDictionary<Monosaccharide, int>> ReadFilter(string path)
        {

            Regex HexNAc = new Regex("HexNAc\\((\\d+)\\)", RegexOptions.Compiled);
            Regex Hex = new Regex("Hex\\((\\d+)\\)", RegexOptions.Compiled);
            Regex Fuc = new Regex("Fuc\\((\\d+)\\)", RegexOptions.Compiled);
            Regex NeuAc = new Regex("NeuAc\\((\\d+)\\)", RegexOptions.Compiled);
            Regex NeuGc = new Regex("NeuGc\\((\\d+)\\)", RegexOptions.Compiled);

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
                            MatchCollection matches = Hex.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.Hex] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (Fuc.IsMatch(line))
                        {
                            MatchCollection matches = Fuc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.Fuc] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (NeuAc.IsMatch(line))
                        {
                            MatchCollection matches = NeuAc.Matches(line);
                            foreach (Match match in matches)
                            {
                                GroupCollection groups = match.Groups;
                                temp[Monosaccharide.NeuAc] = int.Parse(groups[1].Value);
                                break;
                            }
                        }
                        if (NeuGc.IsMatch(line))
                        {
                            MatchCollection matches = NeuGc.Matches(line);
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
            string glycanPath = @"C:\Users\iruiz\Downloads\MSMS\309 N-glycan search space for HGI.txt";
            List<SortedDictionary<Monosaccharide, int>> glycanList
                = ReadFilter(glycanPath);

            GlycanBuilderFiltered glycanBuilder =
                new GlycanBuilderFiltered(glycanList, 7, 7, 5, 4, 0, true, false, false);
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
            Dictionary<double, Dictionary<FragmentTypes, List<string>>> fragments
                = new Dictionary<double, Dictionary<FragmentTypes, List<string>>>();
            List<Tuple<string, FragmentTypes, List<double>>> fragmentsContainer
                = new List<Tuple<string, FragmentTypes, List<double>>>();
            GlycanIonsBuilder.Build.Permethylated = false;
            GlycanIonsBuilder.Build.Derivatization = Glycan.kWater;
            Glycan.To.SetPermethylation(false, false);
            Glycan.To.Derivatization = Glycan.kWater;


            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<Tuple<string, FragmentTypes, List<double>>> fragmentMass
                        = new List<Tuple<string, FragmentTypes, List<double>>>();
                    foreach (FragmentTypes type in GlycanIonsBuilder.Build.Types)
                    {
                        List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan, type)
                                        .Select(m => Math.Round(m, 4)).ToList();
                        fragmentMass.Add(Tuple.Create(id, type, massList));
                    }

                    lock (obj)
                    {
                        fragmentsContainer.AddRange(fragmentMass);
                    }
                }
            });

            foreach (Tuple<string, FragmentTypes, List<double>> item in fragmentsContainer)
            {
                string id = item.Item1;
                FragmentTypes type = item.Item2;

                foreach (double mass in item.Item3)
                {
                    if (!fragments.ContainsKey(mass))
                        fragments[mass] = new Dictionary<FragmentTypes, List<string>>();
                    if (!fragments[mass].ContainsKey(type))
                        fragments[mass][type] = new List<string>();
                    fragments[mass][type].Add(id);
                }
            }

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
