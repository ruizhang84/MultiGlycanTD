using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class GlycanSearch
    {
        ISearch<int> searcher_;

        public GlycanSearch(ISearch<int> searcher)
        {
            searcher_ = searcher;
        }

        void Init(List<IPeak> peaks, int precursorCharge)
        {
            List<Point<int>> points = new List<Point<int>>();
            for(int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                for (int charge = 1; charge <= precursorCharge; charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       util.mass.Spectrum.Proton, charge);
                    points.Add(new Point<int>(mass, i));
                }
            }
            searcher_.Init(points);
        }

        public List<IGlycan> Search(List<IPeak> peaks, int precursorCharge,
            List<IGlycan> candidates)
        {
            Init(peaks, precursorCharge);

            List<IGlycan> results = new List<IGlycan>();

            double bestScore = 0;
            foreach(IGlycan glycan in candidates)
            {
                HashSet<int> matchedIndex = new HashSet<int>();
                List<double> fragments = GlycanIonsBuilder.Build.Fragments(glycan);
                foreach (double mass in fragments)
                {
                    List<int> matched = searcher_.Search(mass, mass);
                    matchedIndex.UnionWith(matched);
                }
                double score = matchedIndex.Select(
                    index => Math.Log(peaks[index].GetIntensity())).Sum();

                // compare score
                if (score > bestScore)
                {
                    bestScore = score;
                    results.Clear();
                    results.Add(glycan);
                }
                else if (score == bestScore)
                {
                    results.Add(glycan);
                }
            }

            return results;
        }

    }
}
