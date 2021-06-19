using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
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
        ISearch<string> searcher_;
        Dictionary<string, List<string>> id_map_;
        Dictionary<double, List<string>> fragments_map_;
        //private readonly int maxCharge = 2;
        private Random random;
        private readonly double lower = 1;
        private readonly double upper = 30;
        public int Seed { get; set; } = 2;
        private int maxCharge = 3; // it is not likely a higher charge for fragments.
        private double tol;
        private ToleranceBy by;

        public GlycanSearch(ISearch<string> searcher, GlycanJson glycanJson)
        {
            searcher_ = searcher;
            id_map_ = glycanJson.IDMap;
            fragments_map_ = glycanJson.Fragments;
            List<Point<string>> points = new List<Point<string>>();
            foreach (double mass in fragments_map_.Keys)
            {
                foreach (string glycan in fragments_map_[mass])
                {
                    points.Add(new Point<string>(mass, glycan));
                }
            }
            searcher_.Init(points);
            random = new Random(Seed);
            tol = searcher.Tolerance();
            by = searcher.ToleranceType();
        }

        public List<SearchResult> Search(List<IPeak> peaks, int precursorCharge,
            List<string> candidates, double ion = 1.0078, bool decoy = false)
        {
            // process composition, id -> compos
            Dictionary<string, string> glycanCandid = new Dictionary<string, string>();
            foreach (string composition in candidates)
            {
                foreach (string glycan in id_map_[composition])
                {
                    glycanCandid[glycan] = composition;
                }
            }

            // search peaks glycan_id -> peak_index -> expect mz
            Dictionary<string, Dictionary<int, double>> matched =
                new Dictionary<string, Dictionary<int, double>>();
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                double randomMass = 0;
                if (decoy)
                    randomMass = random.NextDouble() * (upper - lower) + lower;
                for (int charge = 1; charge <= Math.Min(maxCharge, precursorCharge); charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       ion, charge);
                    if (decoy)
                        mass += randomMass;

                    List<Point<string>> glycans = searcher_.Search(mass);
                    foreach (Point<string> pt in glycans)
                    {
                        string glycan = pt.Content();
                        double expectMZ = util.mass.Spectrum.To.ComputeMZ(pt.Value(), ion, charge);
                        if (!glycanCandid.ContainsKey(glycan))
                            continue;

                        if (!matched.ContainsKey(glycan))
                        {
                            matched[glycan] = new Dictionary<int, double>();
                        }

                        // update matching peaks and expected mass
                        if (!matched[glycan].ContainsKey(i) ||
                            Math.Abs(peak.GetMZ() - expectMZ) < Math.Abs(peak.GetMZ() - matched[glycan][i]))
                        {
                            matched[glycan][i] = expectMZ;
                        }
                    }
                }
            }

            // compute score
            Dictionary<string, SearchResult> results = new Dictionary<string, SearchResult>();
            double bestScore = 0;

            foreach (string isomer in matched.Keys)
            {
                string glycan = glycanCandid[isomer];
                double score = matched[isomer].Select(
                    pair => Math.Log(peaks[pair.Key].GetIntensity())).Sum();
                // compare score
                if (score > bestScore)
                {
                    bestScore = score;
                    results.Clear();
                }
                else if (score < bestScore)
                {
                    continue;
                }
                // build up results
                if (!results.ContainsKey(glycan))
                {
                    results[glycan] = new SearchResult();
                    results[glycan].set_glycan(glycan);
                    results[glycan].set_score(ComputeScore(peaks, matched[isomer]));
                }
                results[glycan].Add(isomer);
            }

            return results.Values.ToList();
        }

        private double Difference(double expect, double obs)
        {
            if (by == ToleranceBy.PPM)
            {
                return util.mass.Spectrum.To.ComputePPM(expect, obs);
            }
            return Math.Abs(expect - obs);
        }

        public double ComputeScore(List<IPeak> peaks, Dictionary<int, double> matches)
        {
            double sum = peaks.Select(p => Math.Sqrt(p.GetIntensity())).Sum();
            double score = matches.Select(
                    pair => 
                    Math.Sqrt(peaks[pair.Key].GetIntensity())
                    * (1 - Math.Pow(Difference(pair.Value, peaks[pair.Key].GetMZ()) / tol, 4))
                    ).Sum();
            return score / sum;
        }

    }
}