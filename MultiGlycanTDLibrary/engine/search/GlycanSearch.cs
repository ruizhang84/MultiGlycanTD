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
    using GlycanFragments = Dictionary<FragmentTypes, List<string>>;

    public class GlycanSearch
    {
        ISearch<GlycanFragments> searcher_;
        Dictionary<string, List<string>> id_map_;
        readonly double minMZ = 200.0; // it is not useful matched at low m/z
        readonly int maxCharge = 3; // it is not likely a higher charge for fragments.
        readonly int minMatches = 5; // it is not likely only match a few peaks.

        public GlycanSearch(
            ISearch<GlycanFragments> searcher, 
            GlycanJson glycanJson)
        {
            searcher_ = searcher;
            List<Point<GlycanFragments>> points
                = new List<Point<GlycanFragments>>();
            foreach (double mass in glycanJson.FragmentMap.Keys)
            {
                points.Add(new Point<GlycanFragments>
                    (mass, glycanJson.FragmentMap[mass]));
            }
            searcher_.Init(points);
            id_map_ = glycanJson.IDMap;
        }

        protected void UpdateMatch(PeakMatch match, IPeak peaks,
            FragmentTypes type, int potentialMatches, double expectMZ)
        {
            double diff = Math.Abs(expectMZ - peaks.GetMZ());
            double prevDiff = Math.Abs(match.TheoreticMZ - peaks.GetMZ());
            if (prevDiff > diff)
            {
                match.TheoreticMZ = expectMZ;
                match.Potentials = potentialMatches;
                match.IonTypes.Clear();
                match.IonTypes.Add(type);
            }
            else if (prevDiff == diff)
            {
                match.Potentials = Math.Min(match.Potentials, potentialMatches);
                match.IonTypes.Add(type);
            }

        }

        protected List<SearchResult> PickTop(
            Dictionary<string, SearchResult> results)
        {
            List<SearchResult> topResults = new List<SearchResult>();
            foreach (string glycan in results.Keys)
            {
                SearchResult result = results[glycan];
                if (result.Matches.Count < minMatches)
                    continue;
                topResults.Add(result);
            }

            return topResults;
        }

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

            // search peaks glycan_id->peak_index
            Dictionary<string, SearchResult> results
                = new Dictionary<string, SearchResult>();
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                if (peak.GetMZ() < minMZ)
                    continue;

                double tempCharge = precursorCharge;
                if (i < peaks.Count - 1 &&
                    Math.Abs(peaks[i + 1].GetMZ() - 1.0 - peak.GetMZ()) < 0.1)
                    tempCharge = 1.0;
                for (int charge = 1; charge <= Math.Min(maxCharge, tempCharge); charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       ion, charge);
                    List<Point<GlycanFragments>> glycans = searcher_.Search(mass);

                    // make records
                    foreach (Point<GlycanFragments> pt in glycans)
                    {
                        GlycanFragments fragments = pt.Content();
                        double expectMZ = util.mass.Spectrum.To.ComputeMZ(pt.Value(), ion, charge);
                        foreach (FragmentTypes type in fragments.Keys)
                        {
                            foreach (string glycan in fragments[type])
                            {
                                if (!glycanCandid.ContainsKey(glycan))
                                {
                                    continue;
                                }

                                if (!results.ContainsKey(glycan))
                                {
                                    results[glycan] = new SearchResult();
                                    results[glycan].Glycan = glycan;
                                    results[glycan].Composition = glycanCandid[glycan];
                                }

                                if (!results[glycan].Matches.ContainsKey(i))
                                {
                                    results[glycan].Matches[i] = new PeakMatch();
                                    results[glycan].Matches[i].Peak = peaks[i];
                                }

                                UpdateMatch(results[glycan].Matches[i], 
                                    peaks[i], type, fragments[type].Count, expectMZ);
                            }
                        }
                    }
                }
            }

            // pick the top candidates
            return PickTop(results);
        }       

    }
}