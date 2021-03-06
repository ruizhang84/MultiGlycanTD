using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class GlycanBuilderFilterUnitTest
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
        public void BuildTest()
        {
            GlycanBuilder glycanBuilder =
                new GlycanBuilder(6, 6, 5, 4, 0, true, false, false);
            glycanBuilder.Build();

            string glycanPath = @"C:\Users\iruiz\Downloads\MSMS\309 N-glycan search space for HGI.txt";
            List<SortedDictionary<Monosaccharide, int>> glycanList
                = ReadFilter(glycanPath);

            Assert.AreEqual(287, glycanList.Count);

            foreach (var item in glycanList)
            {
                string output = "";
                foreach(var pair in item)
                {
                    output += pair.Key.ToString() + "-" + pair.Value.ToString() + " ";
                }
                Console.WriteLine(output);
            }


        }
    }
}
