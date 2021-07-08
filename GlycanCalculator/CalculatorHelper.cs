using MultiGlycanTDLibrary.model.glycan;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace GlycanCalculator
{
    public class CalculatorHelper
    {
        public static List<SortedDictionary<Monosaccharide, int>> ReadFilter(string path)
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
    }
}
