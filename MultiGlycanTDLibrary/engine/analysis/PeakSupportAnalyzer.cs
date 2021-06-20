using MultiGlycanTDLibrary.algorithm;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class PeakSupport
    {
        public int Scan { get; set; }
        public IPeak Peak { get; set; }
        public List<string> Isomers { get; set; }
    }

    public class PeakSupportAnalyzer
    {
        private int maxDistance = 10;

        ISearch<PeakSupport> searcher;
        ToleranceBy by = ToleranceBy.Dalton;
        double tol = 0.1;


        protected int[] ConvertID(string id)
        {
            return id.Split(" ").Select(s => int.Parse(s)).ToArray();
        }

        protected int Distance(string id1, string id2)
        {
            int[] table1 = ConvertID(id1);
            int[] table2 = ConvertID(id2);
            // the same glycan type
            if (table1.Length != table2.Length)
                return int.MaxValue;
            //// the same number of branches
            //if (table1.Length == 26) // Complex
            //{
            //    for (int i = 0; i < 4; i++)
            //    {
            //        // attached at least one GlcNAc, meaning a branch exist
            //        if ((table1[i + 6] > 0) != (table2[i + 6] > 0))
            //            return int.MaxValue;
            //    }
            //}
            //else if (table1.Length == 18)  // hybrid
            //{
            //    for (int i = 0; i < 2; i++)
            //    {
            //        // attached at least one Man, meaning a branch exist
            //        if ((table1[i + 6] > 0) != (table2[i + 6] > 0))
            //            return int.MaxValue;
            //        if ((table1[i + 8] > 0) != (table2[i + 8] > 0))
            //            return int.MaxValue;
            //    }
            //}
            // compute the difference of monosaccharide
            int distance = 0;
            for (int i = 0; i < table1.Length; i++)
            {
                distance += Math.Abs(table1[i] - table2[i]);
            }


            return distance;
        }

        public void Init(List<SearchResult> target)
        {
            searcher = new BucketSearch<PeakSupport>(by, tol);
            List<Point<PeakSupport>> points = new List<Point<PeakSupport>>();
            foreach(SearchResult result in target)
            {
                int scan = result.Scan();
                Dictionary<double, IPeak> matched = result.Matches();
                foreach (double mz in matched.Keys)
                {
                    PeakSupport peakSupport = new PeakSupport()
                    {
                        Scan = scan,
                        Peak = matched[mz],
                        Isomers = result.Isomers()
                    };
                    Point<PeakSupport> point = new Point<PeakSupport>(mz, peakSupport);
                    points.Add(point);
                }
            }
            searcher.Init(points);
        }

        public void Analyze(List<SearchResult> results)
        {
            foreach(SearchResult result in results)
            {
                Analyze(result);
            }
        }

        public void Analyze(SearchResult result)
        {
            double bestScore = 0;
            List<string> bestIsomers = new List<string>();

            foreach (string isomer in result.Isomers())
            {
                Dictionary<IPeak, double> peakSupportFrac = new Dictionary<IPeak, double>();
                foreach (double mz in result.Matches().Keys)
                {
                    // count peak support matched mz
                    int count = 0;
                    List<PeakSupport> supports = searcher.SearchContent(mz);
                    if (supports.Count == 0)
                        continue;
                    foreach (PeakSupport supp in supports)
                    {
                        // not the same spectrum
                        if (supp.Scan == result.Scan())
                            continue;
                        // not exceed max distance to consider
                        foreach (string supportIsomer in supp.Isomers)
                        {
                            if (Distance(isomer, supportIsomer) <= maxDistance)
                            {
                                count++;
                                break;
                            }
                        }

                    }

                    peakSupportFrac[result.Matches()[mz]] = count * 1.0 / supports.Count;
                }

                double score = peakSupportFrac.Select(pair =>
                        pair.Key.GetIntensity() * (1.0 - pair.Value)).Sum();
                if (score > bestScore)
                {
                    bestIsomers.Clear();
                    bestScore = score;
                    bestIsomers.Add(isomer);
                }
                else if (score == bestScore)
                {
                    bestIsomers.Add(isomer);
                }
            }
            result.set_isomers(bestIsomers);
            result.set_score(bestScore);
        }
    }
}
