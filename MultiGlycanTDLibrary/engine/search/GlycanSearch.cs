using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using SpectrumData;
using SpectrumData.Spectrum;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class PeakMatch
    {
        public IPeak Peak { get; set; }
        // theortical m/z
        public double Expected { get; set; }
        // potential matched to other glycans
        public int Potentials { get; set; }
        // fragment ion types
        public List<FragmentTypes> IonTypes { get; set; }
    }

    public class GlycanSearch
    {
        ISearch<string> searcher_;
        Dictionary<string, List<string>> id_map_;
        Dictionary<double, List<string>> fragments_map_;
        private int maxCharge = 3; // it is not likely a higher charge for fragments.
        private int minMatches = 5; // it is not likely only match a few peaks.
        private double tol;
        private ToleranceBy by;

        public GlycanSearch(ISearch<string> searcher, GlycanJson glycanJson)
        {
            searcher_ = searcher;
            id_map_ = glycanJson.IDMap;
            foreach (double mass in glycanJson.Fragments.Keys)
            {
                fragments_map_[mass] = glycanJson.Fragments[mass][FragmentTypes.Y];
            }
            
            List<Point<string>> points = new List<Point<string>>();
            foreach (double mass in fragments_map_.Keys)
            {
                foreach (string glycan in fragments_map_[mass])
                {
                    points.Add(new Point<string>(mass, glycan));
                }
            }
            searcher_.Init(points);
            tol = searcher.Tolerance();
            by = searcher.ToleranceType();
        }

        //protected void UpdateMatchInfo(Dictionary<string, MatchInfo> matched,
        //    string glycan, int index, double observedMZ, double expectMZ)
        //{
        //    if (!matched.ContainsKey(glycan))
        //    {
        //        matched[glycan] = new MatchInfo();
        //    }

        //    // updat info
        //    matched[glycan].Peaks.Add(index);

        //    // update matching peaks and expected mass
        //    if (!matched[glycan].Expects.ContainsKey(index) ||
        //        Math.Abs(observedMZ - expectMZ) < Math.Abs(observedMZ - matched[glycan].Expects[index]))
        //    {
        //        matched[glycan].Expects[index] = expectMZ;
        //    }
        //}

        //protected Dictionary<string, SearchResult> ComputeResults(
        //    Dictionary<string, MatchInfo> matched,
        //    Dictionary<string, string> glycanCandid,
        //    List<IPeak> peaks)
        //{
        //    Dictionary<string, SearchResult> results = new Dictionary<string, SearchResult>();
        //    double bestScore = 0;
        //    foreach (string isomer in matched.Keys)
        //    {
        //        string glycan = glycanCandid[isomer];
        //        double score = matched[isomer].Peaks.Select(index =>
        //            Math.Log10(peaks[index].GetIntensity())).Sum();

        //        // compare score
        //        if (score > bestScore)
        //        {
        //            bestScore = score;
        //            results.Clear();
        //        }
        //        else if (score < bestScore || score == 0)
        //        {
        //            continue;
        //        }
        //        // check number of matches
        //        if (matched[isomer].Peaks.Count < minMatches)
        //            continue;

        //        // build up results
        //        if (!results.ContainsKey(glycan))
        //        {
        //            results[glycan] = new SearchResult();
        //            results[glycan].set_glycan(glycan);
        //            results[glycan].set_score(ComputeScore(peaks, matched[isomer]));
        //        }
        //        results[glycan].Add(isomer);
        //    }

        //    return results;
        //}

        public List<SearchResult> Search(List<IPeak> peaks, int precursorCharge,
            List<string> candidates, double ion = 1.0078)
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
            //Dictionary<string, MatchInfo> matched =
            //    new Dictionary<string, MatchInfo>();
            //for (int i = 0; i < peaks.Count; i++)
            //{
            //    IPeak peak = peaks[i];
            //    for (int charge = 1; charge <= Math.Min(maxCharge, precursorCharge); charge++)
            //    {
            //        double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
            //           ion, charge);

            //        List<Point<string>> glycans = searcher_.Search(mass);

            //        // make records
            //        foreach (Point<string> pt in glycans)
            //        {
            //            string glycan = pt.Content();
            //            double expectMZ = util.mass.Spectrum.To.ComputeMZ(pt.Value(), ion, charge);

            //            if (!glycanCandid.ContainsKey(glycan))
            //            {
            //                continue;
            //            }


            //            UpdateMatchInfo(matched, glycan, i, peak.GetMZ(), expectMZ);
            //        }
            //    }
            //}

            // compute score

            return new List<SearchResult>();
        }

        private double Difference(double expect, double obs)
        {
            if (by == ToleranceBy.PPM)
            {
                return util.mass.Spectrum.To.ComputePPM(expect, obs);
            }
            return Math.Abs(expect - obs);
        }

        //public double ComputeScore(List<IPeak> peaks, MatchInfo match)
        //{
        //    double sum = peaks.Select(p => Math.Sqrt(p.GetIntensity())).Sum();
        //    double score = match.Peaks.Select(
        //            index =>
        //            Math.Sqrt(peaks[index].GetIntensity())
        //            * (1 - Math.Pow(Difference(match.Expects[index], peaks[index].GetMZ()) / tol, 4))
        //            ).Sum();
        //    return score / sum;
        //}

    }
}